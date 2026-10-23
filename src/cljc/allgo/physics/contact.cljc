(ns allgo.physics.contact
  "Where two bodies touch, and how deeply.

  `allgo.geometry.gjk` already answers whether two convex shapes overlap
  and by how much, and that is not enough to stand a wall of bricks up.
  GJK with EPA gives *one* deepest point and a normal. A brick resting on
  another touches it along a whole face, and a solver given a single
  point has to be told about the other three by the next frame's
  correction -- so the brick rocks, settles, rocks again, and a wall of
  them shivers instead of standing.

  What a solver wants is a *manifold*: every point where the two shapes
  are in contact, produced in one go, so that a face resting on a face
  arrives as four points and holds still. That is what this makes.

  ## How a box meets a box

  The separating axis theorem, then clipping.

  Two convex shapes miss each other exactly when some axis exists on
  which their shadows do not overlap. For boxes the only axes worth
  testing are fifteen: three face directions each, and the nine cross
  products of their edge directions. If every one of them overlaps the
  boxes intersect, and the axis with the *least* overlap is the direction
  to push them apart along -- the shortest way out.

  What that axis is tells you the kind of contact:

    a face axis   one box has a face pointing at the other. Take that
                  face as the reference, find the most opposed face on
                  the other box, clip it against the reference's four
                  side planes, and keep whatever lies below the reference
                  plane. Up to four points, which is a resting face.
    an edge axis   two edges crossing in mid-air, like the corner of one
                  brick on the corner of another. That really is a single
                  point, and it is the closest point between the two
                  edge segments.

  Face axes are given a small preference when the overlaps are nearly
  equal. Without it a box resting flat on another flips between a face
  result and an edge result as the numbers jitter in the last bits, and
  four contact points become one and back again -- which the solver sees
  as the support vanishing.

  ## What comes out

  A vector of contacts, each

      {:a :b        body indices, normal points from `a` toward `b`
       :point       where they touch, in world coordinates
       :normal      unit, the direction to separate along
       :depth       signed -- positive is an overlap, negative a gap}

  which is what `allgo.physics.solver` consumes whichever solver it is
  running.

  ## Contacts before the touch

  A contact is offered while the two are still apart, out to `margin`,
  and that margin is as wide as the pair can travel in the step. It is
  what stops a fast thing going through a thin thing: detection runs
  once a step, so a ball crossing a meter and a half between two looks
  is never near a 2cm window and arrives on the far side of the brick
  having been asked about nothing. The broad phase sweeps its boxes
  over the step for the same reason -- a pair has to survive it before
  the margin gets a chance to matter.

  This is the cheap half of continuous collision detection and it is
  not the whole of it. It puts a gap in front of the solver in time to
  be seen; it does not find the moment of impact, so a body is stopped
  at the surface within a step rather than at the instant it arrives,
  and a solver that integrates before it projects can still step over
  one. `allgo.physics.solver` says which does."
  (:require [allgo.array :as a]
            [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.math :as am]
            [allgo.numerics.roots :as roots]
            [allgo.physics.rigid :as rigid]
            [allgo.spatial.sweep :as sweep]
            [clojure.math :as math]))

(def ^:private eps 1e-9)

;; ---------------------------------------------------------------------------
;; Boxes

(defn- half [{:keys [size]}] (mapv #(* 0.5 (double %)) size))

(defn frame
  "A box reduced to what the tests actually ask it: where it is, which
  way its three faces point, and how far it reaches along each.

  Computed once per body rather than per axis. The separating axis test
  asks fifteen questions of a pair and each one used to re-derive both
  boxes' orientations -- ninety quaternion rotations to answer a question
  about two boxes that had not moved in between. Hoisting it took a
  225-body wall from 11.7ms of narrow phase to a fraction of that."
  [{:keys [pos rot] :as b}]
  {:pos pos
   :half (half b)
   :axes [(q/rotate rot [1.0 0.0 0.0])
          (q/rotate rot [0.0 1.0 0.0])
          (q/rotate rot [0.0 0.0 1.0])]})

(def speculative
  "How far apart two bodies can be and still be given a contact.

  A contact that only exists once the bodies already overlap can only
  ever clean up after a collision it was too late to prevent: the solver
  is handed a penetration and has to push it out, and if the stack is
  deep enough that it cannot push it all out in the iterations it has,
  the overlap grows until the separating axis test picks a different axis
  and the two pass through each other. That is how a twenty brick column
  telescopes into the floor.

  A speculative contact is made *before* the touch, carries a negative
  depth -- a gap -- and asks the solver for something weaker than `do not
  overlap`: approach no faster than closes that gap this step. A body
  flying at a wall is then stopped exactly at the surface rather than
  after it. Four times the solver's slop, which is Box2D's figure.

  This is the floor, not the whole margin: what detection actually
  allows is this plus how far the pair can travel before it looks
  again. See `margin`."
  0.02)

(defn travel
  "How far `b` moves this step, as a bound on what detection can miss.

  Linear only. A spinning body's surface moves too, but a contact is
  missed when the *gap closes* faster than detection looks, and closing
  is what translation does -- a brick can spin as fast as it likes about
  its own center without arriving anywhere. Adding the angular term
  would widen every long box's margin for a case that does not tunnel."
  ^double [b ^double dt]
  (if (rigid/inert? b)
    0.0
    (* dt (v/length (:vel b)))))

(defn margin
  "How far apart `a` and `b` can be and still be given a contact.

  The fixed margin plus everything the two of them can do to close it
  before detection runs again. Without the second term the margin is a
  statement about geometry when the thing it has to survive is speed: a
  ball crossing 1.5m in a step steps clean over a 2cm window, is never
  offered a contact at all, and arrives on the far side of the brick.
  With it the gap is seen while it is still a gap, and the solver's
  existing rule for gaps -- approach no faster than closes this one --
  stops the ball at the surface.

  It is bought, not free: a margin meters wide asks the separating axis
  test for the shortest way out between two bodies that are nowhere near
  each other, and that answer is only roughly the direction they will
  actually meet from. See `allgo.physics.solver` on what that costs."
  ^double [a b ^double dt]
  (+ speculative (travel a dt) (travel b dt)))

(def ^:private clip-eps
  "A corner this close to a clip plane counts as inside it.

  Two identical bricks stacked square have all four corners of the
  incident face sitting exactly on the reference face's four side planes.
  Clipped strictly, each of those is a cut rather than a corner, the
  answer comes back with six points instead of four, and which six
  depends on the last bit of the arithmetic -- so the manifold is a
  different manifold every step and nothing can be warm started. Kept as
  corners, the pair reports the same four points, with the same ids,
  for as long as it stands there."
  1e-9)

(defn- clip!
  "Sutherland-Hodgman: the part of a polygon on the inner side of a plane.

  The polygon is the first `cnt` points of `src`, flat, with their ids in
  `sids`; the plane is `[nx ny nz d]`, inside where `(dot n x) <= d`.
  What survives goes into `dst` and `dids`, and the count is returned.

  The ids are the point of it. A corner that survives the clip keeps the
  id it came in with; a point made by cutting an edge gets `[tag id]`,
  after the plane that cut it and the corner the cut started from. So the
  same configuration next step produces the same ids, and a contact can be
  recognized as the one that was there before rather than merely as one
  nearby."
  [^doubles src ^objects sids cnt ^doubles dst ^objects dids ^doubles plane tag]
  (let [cnt (long cnt)
        snx (aget plane 0) sny (aget plane 1) snz (aget plane 2) sd (aget plane 3)]
    (loop [i 0 out 0]
      (if (= i cnt)
        out
        (let [q (rem (inc i) cnt)
              i3 (* 3 i) q3 (* 3 q)
              px (aget src i3) py (aget src (+ i3 1)) pz (aget src (+ i3 2))
              qx (aget src q3) qy (aget src (+ q3 1)) qz (aget src (+ q3 2))
              dp (- (+ (* snx px) (* sny py) (* snz pz)) sd)
              dq (- (+ (* snx qx) (* sny qy) (* snz qz)) sd)
              p-in? (<= dp clip-eps)
              q-in? (<= dq clip-eps)
              o3 (* 3 out)]
          (when (not= p-in? q-in?)
            ;; The cut, where the edge crosses the plane: `v/lerp`.
            (let [t (/ dp (- dp dq))]
              (aset dst o3 (+ px (* (- qx px) t)))
              (aset dst (+ o3 1) (+ py (* (- qy py) t)))
              (aset dst (+ o3 2) (+ pz (* (- qz pz) t)))
              (aset dids out [tag (aget sids i)])))
          (cond
            (and p-in? q-in?)
            (do (aset dst o3 qx) (aset dst (+ o3 1) qy) (aset dst (+ o3 2) qz)
                (aset dids out (aget sids q))
                (recur (inc i) (inc out)))
            p-in? (recur (inc i) (inc out))
            q-in?
            (let [o3 (+ o3 3)]
              (aset dst o3 qx) (aset dst (+ o3 1) qy) (aset dst (+ o3 2) qz)
              (aset dids (inc out) (aget sids q))
              (recur (inc i) (+ out 2)))
            :else (recur (inc i) out)))))))

(defn- face-contacts
  "Contacts from a face of `ref-body` pressing into `inc-body`.

  The reference face's four side planes bound the region the contact can
  lie in; clipping the opposing face to them and keeping what is still
  below the reference plane gives the resting polygon."
  ;; On primitive doubles, the polygon clipped back and forth between two
  ;; small buffers. Written over vec3s, with the polygon a vector of
  ;; `[point id]` pairs, this was most of the narrow phase once the axis
  ;; search was on primitives -- it runs for every touching pair. The
  ;; rewrite kept that version's arithmetic in its order, so the points
  ;; and their ids came out the same to the last bit.
  [ia ib ref-body inc-body ref-axis ref-sign n flip? margin]
  (let [margin (double margin)
        ref-axis (long ref-axis)
        ref-sign (double ref-sign)
        rh (:half ref-body) ras (:axes ref-body) [rpx rpy rpz] (:pos ref-body)
        rpx (double rpx) rpy (double rpy) rpz (double rpz)
        ;; The reference face: its outward normal, and its plane.
        [rax ray raz] (nth ras ref-axis)
        rnx (* (double rax) ref-sign) rny (* (double ray) ref-sign) rnz (* (double raz) ref-sign)
        rha (double (nth rh ref-axis))
        ref-d (+ (* rnx (+ rpx (* rnx rha))) (* rny (+ rpy (* rny rha))) (* rnz (+ rpz (* rnz rha))))
        ias (:axes inc-body) ih (:half inc-body) [ipx ipy ipz] (:pos inc-body)
        ipx (double ipx) ipy (double ipy) ipz (double ipz)
        ;; The incident face: the one of the other box pointing most
        ;; directly back along the reference normal -- the face being
        ;; pressed into. Ties go to the first in order, compared as
        ;; `compare` does, which is how the sort this replaced chose.
        k (loop [k 0 best-k 0 best 0.0]
            (if (= k 6)
              best-k
              (let [[ux uy uz] (nth ias (quot k 2))
                    sign (if (even? k) 1.0 -1.0)
                    dv (+ (* (* (double ux) sign) rnx) (* (* (double uy) sign) rny) (* (* (double uz) sign) rnz))]
                (if (or (zero? k) (neg? (compare dv best)))
                  (recur (inc k) k dv)
                  (recur (inc k) best-k best)))))
        inc-axis (quot (long k) 2)
        inc-sign (if (even? (long k)) 1.0 -1.0)
        ;; Its four corners, counter-clockwise about its normal, ids 0 to
        ;; 3. A quadrilateral clipped by four planes gains at most a corner
        ;; a plane, so eight is room enough.
        ^doubles pa (a/f64 24) ^doubles pb (a/f64 24)
        ^objects ida (object-array 8) ^objects idb (object-array 8)
        _ (let [ii (mod (inc inc-axis) 3)
                jj (mod (inc ii) 3)
                hi (double (nth ih ii)) hj (double (nth ih jj)) hk (double (nth ih inc-axis))
                [nx ny nz] (nth ias inc-axis)
                [ux uy uz] (nth ias ii)
                [wx wy wz] (nth ias jj)
                ux (* (double ux) hi) uy (* (double uy) hi) uz (* (double uz) hi)
                wx (* (double wx) hj) wy (* (double wy) hj) wz (* (double wz) hj)
                cx (+ ipx (* (* (double nx) inc-sign) hk))
                cy (+ ipy (* (* (double ny) inc-sign) hk))
                cz (+ ipz (* (* (double nz) inc-sign) hk))]
            ;; c - u - w, c + u - w, c + u + w, c - u + w. Adding the
            ;; negation is subtracting, exactly.
            (dotimes [c 4]
              (let [su (if (or (= c 0) (= c 3)) -1.0 1.0)
                    sw (if (< c 2) -1.0 1.0)
                    c3 (* 3 c)]
                (aset pa c3 (+ (+ cx (* su ux)) (* sw wx)))
                (aset pa (+ c3 1) (+ (+ cy (* su uy)) (* sw wy)))
                (aset pa (+ c3 2) (+ (+ cz (* su uz)) (* sw wz)))
                ;; Boxed: the ids are objects, and cut points' ids are
                ;; vectors.
                (aset ida c (identity c)))))
        ;; The reference face's four sides, as outward planes tagged 0 to
        ;; 3: each axis but the reference one, outward then inward.
        k0 (if (= ref-axis 0) 1 0)
        k1 (if (= ref-axis 2) 1 2)
        ^doubles plane (a/f64 4)
        [cnt ^doubles pts ^objects ids]
        (loop [t 0 cnt 4 src pa sids ida dst pb dids idb]
          (if (or (= t 4) (zero? (long cnt)))
            [cnt src sids]
            (let [kk (if (< t 2) k0 k1)
                  s (if (even? t) 1.0 -1.0)
                  [ux uy uz] (nth ras kk)
                  ux (double ux) uy (double uy) uz (double uz)
                  sh (* s (double (nth rh kk)))
                  snx (* ux s) sny (* uy s) snz (* uz s)]
              (aset plane 0 snx) (aset plane 1 sny) (aset plane 2 snz)
              (aset plane 3 (+ (* snx (+ rpx (* ux sh))) (* sny (+ rpy (* uy sh))) (* snz (+ rpz (* uz sh)))))
              (recur (inc t) (clip! src sids cnt dst dids plane t) dst dids src sids))))
        cnt (long cnt)
        normal (if flip? (v/negate n) n)
        ^longs kept (long-array 8)]
    ;; A face clipped against a face it exactly coincides with -- two
    ;; identical bricks, square on -- has every corner sitting on a clip
    ;; plane, and Sutherland-Hodgman answers each of those twice. Keeping
    ;; both gives a pair six or eight contact points where it has four,
    ;; and which of them survive changes from step to step. One point per
    ;; place, first id wins, and the ordering is the incident face's own,
    ;; so the choice is the same every step. A point is judged a duplicate
    ;; before its depth is: one that matches a point kept but then found
    ;; too far away to report is still dropped.
    (loop [i 0 nkept 0 out (transient [])]
      (if (= i cnt)
        (persistent! out)
        (let [i3 (* 3 i)
              px (aget pts i3) py (aget pts (+ i3 1)) pz (aget pts (+ i3 2))
              dup? (loop [m 0]
                     (if (= m nkept)
                       false
                       (let [j3 (* 3 (aget kept m))
                             dx (- px (aget pts j3)) dy (- py (aget pts (+ j3 1))) dz (- pz (aget pts (+ j3 2)))]
                         (if (< (math/sqrt (+ (* dx dx) (* dy dy) (* dz dz))) 1e-6)
                           true
                           (recur (inc m))))))]
          (if dup?
            (recur (inc i) nkept out)
            (let [_ (aset kept nkept i)
                  depth (- ref-d (+ (* rnx px) (* rny py) (* rnz pz)))]
              (if (>= depth (- margin))
                (let [hd (* 0.5 depth)]
                  (recur (inc i) (inc nkept)
                         (conj! out
                                {:a ia :b ib
                                   ;; Reported on the reference surface
                                   ;; rather than at the clipped point, so
                                   ;; both bodies agree where the touch is.
                                 :point [(+ px (* rnx hd)) (+ py (* rny hd)) (+ pz (* rnz hd))]
                                 :normal normal
                                   ;; Which corner of which face, pressed
                                   ;; into which face -- the same answer
                                   ;; every step for as long as the two
                                   ;; bodies stay in the same arrangement,
                                   ;; which is what warm starting needs to
                                   ;; recognize it by. `flip?` is in it
                                   ;; because it says which body was the
                                   ;; reference, and that decides what the
                                   ;; rest of the id means.
                                 :id [(boolean flip?) ref-axis ref-sign inc-axis inc-sign (aget ids i)]
                                   ;; Signed: positive is an overlap,
                                   ;; negative a gap the solver is allowed
                                   ;; to see coming.
                                 :depth depth})))
                (recur (inc i) (inc nkept) out)))))))))

(defn- closest-on-segments
  "The midpoint of the shortest link between two segments."
  [p1 q1 p2 q2]
  (let [d1 (v/sub q1 p1)
        d2 (v/sub q2 p2)
        r (v/sub p1 p2)
        a (v/dot d1 d1)
        e (v/dot d2 d2)
        f (v/dot d2 r)
        c (v/dot d1 r)
        b (v/dot d1 d2)
        denom (- (* a e) (* b b))
        s (if (> (abs denom) eps)
            (am/clamp (/ (- (* b f) (* c e)) denom) 0.0 1.0)
            0.0)
        t (am/clamp (/ (+ (* b s) f) (max eps e)) 0.0 1.0)
        s (am/clamp (/ (- (* b t) c) (max eps a)) 0.0 1.0)]
    (v/scale (v/add (v/add-scaled p1 d1 s) (v/add-scaled p2 d2 t)) 0.5)))

(defn- edge-contact
  "The single point where two crossing edges meet."
  [ia ib a b n depth i j]
  (let [pick (fn [body axis-i axis-j]
               ;; The corner furthest along the contact normal decides
               ;; which of the four parallel edges is the one touching.
               (let [h (:half body) as (:axes body)
                     k (first (remove #{axis-i axis-j} (range 3)))
                     sgn (fn [idx dir]
                           (v/scale (nth as idx)
                                    (* (nth h idx)
                                       (if (pos? (* dir (v/dot (nth as idx) n))) 1.0 -1.0))))
                     center (v/add (:pos body) (v/add (sgn axis-i 1.0) (sgn axis-j 1.0)))
                     along (v/scale (nth as k) (nth h k))]
                 [(v/sub center along) (v/add center along)]))
        [p1 q1] (pick a (mod (inc (long i)) 3) (mod (inc (long (mod (inc (long i)) 3))) 3))
        [p2 q2] (pick b (mod (inc (long j)) 3) (mod (inc (long (mod (inc (long j)) 3))) 3))]
    [{:a ia :b ib
      :point (closest-on-segments p1 q1 p2 q2)
      :normal n
      ;; One point, named by the pair of edges that made it.
      :id [:edge i j]
      :depth depth}]))

(def ^:private b-face-bias
  "How much less a face of `b` must overlap than a face of `a` to be taken
  as the reference instead, in meters.

  Two identical bricks one on the other give the two boxes' faces nearly
  the same overlap, and a standing wall tilts each brick by fractions of a
  milliradian -- enough to put either face ahead by a fraction of a
  millimeter, and to swap them from one step to the next. Preferring a's
  only on a tie, as the order of the tests used to, the reference changed
  hands that way all the time, every point on the pair was renamed with
  it, and warm starting found last step's impulse for 54% of the contacts
  in a settled 16x16 wall under XPBD, 90% under TGS. With 5mm it finds
  97% and 98%.

  It decides the reference, not the depth: a face of a that is a few
  millimeters deeper is still a fine one to clip against, and columns and
  walls under all three solvers stood no worse for it, beyond the noise
  between runs. Box2D's box test has the same preference with a tenth of
  its linear slop, half a millimeter; here a millimeter left sequential
  impulse and TGS at 94%, and five took them to 96-98%."
  0.005)

(defn box-box
  "The manifold between two boxes, or an empty vector.

  `a` and `b` are `frame`s, not bodies."
  ;; The fifteen axes on primitive doubles: each box's center, extents and
  ;; axes read out once, and the best axis so far carried in the loop
  ;; rather than in a map rebuilt at every improvement. Written over vec3s
  ;; this asked for a dozen vectors per axis -- the axis, its unit vector,
  ;; the offset between the centers, each projection -- and was most of
  ;; the narrow phase. The rewrite kept that version's arithmetic in its
  ;; order, so the answers came out the same to the last bit.
  [ia ib a b margin]
  (let [margin (double margin)
        as (:axes a)
        bs (:axes b)
        [pax pay paz] (:pos a) [pbx pby pbz] (:pos b)
        [hax hay haz] (:half a) [hbx hby hbz] (:half b)
        [[a0x a0y a0z] [a1x a1y a1z] [a2x a2y a2z]] as
        [[b0x b0y b0z] [b1x b1y b1z] [b2x b2y b2z]] bs
        pax (double pax) pay (double pay) paz (double paz)
        pbx (double pbx) pby (double pby) pbz (double pbz)
        hax (double hax) hay (double hay) haz (double haz)
        hbx (double hbx) hby (double hby) hbz (double hbz)
        a0x (double a0x) a0y (double a0y) a0z (double a0z)
        a1x (double a1x) a1y (double a1y) a1z (double a1z)
        a2x (double a2x) a2y (double a2y) a2z (double a2z)
        b0x (double b0x) b0y (double b0y) b0z (double b0z)
        b1x (double b1x) b1y (double b1y) b1z (double b1z)
        b2x (double b2x) b2y (double b2y) b2z (double b2z)
        dpx (- pbx pax) dpy (- pby pay) dpz (- pbz paz)
        ;; Faces are favored over edge pairs, and a's faces over b's, so
        ;; that near-equal overlaps do not flip the answer on noise: from
        ;; a face to an edge pair, and a resting brick loses three of its
        ;; four contact points for a frame; or from a's face to b's, and
        ;; it keeps four points that are named after the other box and so
        ;; match nothing from last step. An edge pair must beat the best
        ;; face of either box by 2%, and a face of b must beat a's by
        ;; `b-face-bias`.
        ;;
        ;; `t` runs over the axes: 0-2 a's faces, 3-5 b's, then the nine
        ;; edge pairs, a's edge `(quot (- t 6) 3)` against b's `(rem ...)`.
        ;; `kind` is 0 for none yet, 1 for a face of a, 2 of b, 3 an edge
        ;; pair; -1 is a separating axis, and ends the search.
        [kind bi bj bnx bny bnz bdepth]
        (loop [t 0 kind 0 bi 0 bj 0 score 0.0 bnx 0.0 bny 0.0 bnz 0.0 bdepth 0.0 fmin ##Inf]
          (if (= t 15)
            [kind bi bj bnx bny bnz bdepth]
            (let [face? (< t 6)
                  ei (quot (- t 6) 3) ej (rem (- t 6) 3)
                  ;; Component by component rather than as a vector, which
                  ;; would be one more allocation per axis.
                  u (long (if face? t ei))
                  ux (case u 0 a0x 1 a1x 2 a2x 3 b0x 4 b1x b2x)
                  uy (case u 0 a0y 1 a1y 2 a2y 3 b0y 4 b1y b2y)
                  uz (case u 0 a0z 1 a1z 2 a2z 3 b0z 4 b1z b2z)
                  wx (case (long ej) 0 b0x 1 b1x b2x)
                  wy (case (long ej) 0 b0y 1 b1y b2y)
                  wz (case (long ej) 0 b0z 1 b1z b2z)
                  ;; A face axis as it stands; an edge pair's cross product.
                  x (if face? ux (- (* uy wz) (* uz wy)))
                  y (if face? uy (- (* uz wx) (* ux wz)))
                  z (if face? uz (- (* ux wy) (* uy wx)))
                  len (math/sqrt (+ (* x x) (* y y) (* z z)))]
              ;; The cross product of two parallel edges is the zero
              ;; vector, no axis at all, and two boxes sitting square to
              ;; each other have three such pairs. Skipped, and not read
              ;; as separating: that reports every axis-aligned box in the
              ;; world as touching nothing.
              (if (<= len 1e-6)
                (recur (inc t) kind bi bj score bnx bny bnz bdepth fmin)
                (let [inv (/ 1.0 len)
                      nx (* x inv) ny (* y inv) nz (* z inv)
                      d (abs (+ (* dpx nx) (* dpy ny) (* dpz nz)))
                      ra (+ (* hax (abs (+ (* a0x nx) (* a0y ny) (* a0z nz))))
                            (* hay (abs (+ (* a1x nx) (* a1y ny) (* a1z nz))))
                            (* haz (abs (+ (* a2x nx) (* a2y ny) (* a2z nz)))))
                      rb (+ (* hbx (abs (+ (* b0x nx) (* b0y ny) (* b0z nz))))
                            (* hby (abs (+ (* b1x nx) (* b1y ny) (* b1z nz))))
                            (* hbz (abs (+ (* b2x nx) (* b2y ny) (* b2z nz)))))
                      o (- (+ ra rb) d)]
                  ;; Within the margin counts as touching. `o` comes back
                  ;; negative there and stays negative all the way to the
                  ;; solver, which is what tells it this is a gap.
                  (if-not (pos? (+ o margin))
                    [-1 0 0 0.0 0.0 0.0 0.0]
                    (let [;; An edge pair is handicapped against a face: it has
                          ;; to fit two percent better to be taken. Scaled the
                          ;; other way for a gap, whose overlap is negative --
                          ;; times 1.02 there made the edge the better fit, and a
                          ;; box settling flat onto the floor within the margin
                          ;; came back as one point on the floor's far edge.
                          sc (cond face? o (pos? o) (* o 1.02) :else (* o 0.98))
                          fmin (if face? (min fmin o) fmin)
                          ;; What the candidate has to beat. A face of b
                          ;; takes over from a face of a only by the bias. An
                          ;; edge pair takes over from a face only by beating
                          ;; the best face of either box: when a's face is
                          ;; kept over a slightly smaller one of b's, a tilted
                          ;; edge pair can come in between the two, and it
                          ;; is a single point where a resting face is four.
                          better? (cond
                                    (zero? kind) true
                                    (and face? (>= t 3) (= kind 1)) (< (+ sc b-face-bias) score)
                                    (and (not face?) (not= kind 3)) (< sc fmin)
                                    :else (< sc score))]
                      (if better?
                        (recur (inc t)
                               (cond (< t 3) 1 face? 2 :else 3)
                               (if face? (rem t 3) ei) ej
                               sc nx ny nz o fmin)
                        (recur (inc t) kind bi bj score bnx bny bnz bdepth fmin)))))))))
        kind (long kind)]
    (if (<= kind 0)
      []
      (let [bnx (double bnx) bny (double bny) bnz (double bnz)
            ;; Point the axis from a toward b.
            n (if (neg? (+ (* bnx dpx) (* bny dpy) (* bnz dpz)))
                [(- bnx) (- bny) (- bnz)]
                [bnx bny bnz])]
        (case kind
          1 (face-contacts ia ib a b bi
                           (if (pos? (v/dot (nth as bi) n)) 1.0 -1.0) n false margin)
          2 (face-contacts ia ib b a bi
                           (if (pos? (v/dot (nth bs bi) (v/negate n))) 1.0 -1.0)
                           (v/negate n) true margin)
          3 (edge-contact ia ib a b n bdepth bi bj))))))

;; ---------------------------------------------------------------------------
;; Spheres

(defn sphere-sphere [ia ib a b margin]
  (let [margin (double margin)
        d (v/sub (:pos b) (:pos a))
        dist (v/length d)
        r (+ (double (:radius a)) (double (:radius b)))]
    (if (or (>= dist (+ r margin)) (< dist eps))
      []
      (let [n (v/scale d (/ 1.0 dist))]
        [{:a ia :b ib
          :point (v/add-scaled (:pos a) n (* 0.5 (+ dist (- (double (:radius a))
                                                            (double (:radius b))))))
          :normal n
          :id :sphere
          :depth (- r dist)}]))))

(defn sphere-box
  "Sphere `a` against box `b`.

  The nearest point on a box to anything is found by clamping in the
  box's own frame, which is the whole of it -- as long as the center is
  outside. A center that has tunneled inside has no nearest surface
  point in that sense, and is pushed out through whichever face it is
  closest to."
  [ia ib a b margin]
  (let [margin (double margin)
        local (q/rotate (:inv-rot b) (v/sub (:pos a) (:pos b)))
        [hx hy hz] (half b)
        [lx ly lz] local
        clamped [(am/clamp lx (- hx) hx) (am/clamp ly (- hy) hy) (am/clamp lz (- hz) hz)]
        inside? (= clamped (vec local))
        r (double (:radius a))]
    (if inside?
      ;; Out through the nearest face.
      (let [gaps [[(- hx (abs lx)) 0 (if (pos? lx) 1.0 -1.0)]
                  [(- hy (abs ly)) 1 (if (pos? ly) 1.0 -1.0)]
                  [(- hz (abs lz)) 2 (if (pos? lz) 1.0 -1.0)]]
            [gap axis sign] (first (sort-by first gaps))
            n-local (assoc [0.0 0.0 0.0] axis sign)
            n (q/rotate (:rot b) n-local)]
        [{:a ia :b ib
          :point (:pos a)
          :normal (v/negate n)
          :id [:inside axis]
          :depth (+ r gap)}])
      (let [world (v/add (q/rotate (:rot b) clamped) (:pos b))
            d (v/sub world (:pos a))
            dist (v/length d)]
        (if (or (>= dist (+ r margin)) (< dist eps))
          []
          [{:a ia :b ib
            :point world
            :normal (v/scale d (/ 1.0 dist))
            :id :sphere-box
            :depth (- r dist)}])))))

;; Tori

(def ^:private torus-samples
  "How many points round the ring are tried before refining. Enough that
  the deepest one is never missed between two of them on anything the
  size of a tire against anything the size of a floor or a curb."
  24)

(defn- ring-point
  "The point at angle `t` on the circle a torus's tube is swept along."
  [{:keys [pos rot major]} ^double t]
  (let [big (double major)]
    (v/add pos (q/rotate rot [(* big (math/cos t)) (* big (math/sin t)) 0.0]))))

(defn- torus-box-deepest
  [ia ib a b margin]
  (let [r (double (:minor a))
        depth-at (fn ^double [^double t]
                   (let [cs (sphere-box ia ib {:pos (ring-point a t) :radius r} b margin)]
                     (if (seq cs) (double (:depth (first cs))) ##-Inf)))
        step (/ (* 2.0 math/PI) torus-samples)
        best (apply max-key #(depth-at (double %))
                    (map #(* step (double %)) (range torus-samples)))]
    (if (= ##-Inf (depth-at best))
      []
      (let [;; golden section between the best sample's neighbors, 24 narrowings
            t (roots/maximize depth-at (- (double best) step) (+ (double best) step) {:tol 0.0 :max-iter 23})
            t (if (>= (depth-at t) (depth-at best)) t best)]
        (mapv #(assoc % :id :torus-box)
              (sphere-box ia ib {:pos (ring-point a t) :radius r} b margin))))))

(defn torus-box
  "Torus `a` against box `b`.

  A torus is every point within `minor` of a circle, so where it meets a
  box is where a ball of radius `minor`, somewhere on that circle, meets
  it -- and `sphere-box` already answers that. What is left is *where*
  on the circle, which is the one-dimensional question of which angle
  goes deepest: try evenly spaced angles, then narrow in on the best by
  golden section between its neighbors.

  One contact, not a manifold. A tire on the road touches in one patch,
  and that is the case this is for; a torus lying flat on a floor touches
  all the way round and will rock until something else is resting on it."
  [ia ib a b margin]
  (if (empty? (sphere-box ia ib {:pos (:pos a) :radius (+ (double (:major a)) (double (:minor a)))}
                          b margin))
    ;; Not even the ball the whole torus fits inside reaches the box, and
    ;; one test says so where finding the deepest point takes seventy.
    []
    (torus-box-deepest ia ib a b margin)))

;; ---------------------------------------------------------------------------
;; Dispatch

(defn between*
  "`between`, given the boxes' `frame`s so a caller with many pairs can
  build each body's once, and the `margin` out to which a gap still
  counts as a contact."
  [ia ib a b fa fb margin]
  (let [margin (double margin)
        sa (:shape a) sb (:shape b)]
    (cond
      (and (= sa :ball) (= sb :ball)) (sphere-sphere ia ib a b margin)
      (and (= sa :ball) (= sb :box)) (sphere-box ia ib a b margin)
      (and (= sa :box) (= sb :ball))
      (mapv (fn [c] (assoc c :a ia :b ib :normal (v/negate (:normal c))))
            (sphere-box ib ia b a margin))
      (and (= sa :box) (= sb :box)) (box-box ia ib fa fb margin)
      (and (= sa :torus) (= sb :box)) (torus-box ia ib a b margin)
      (and (= sa :box) (= sb :torus))
      (mapv (fn [c] (assoc c :a ia :b ib :normal (v/negate (:normal c))))
            (torus-box ib ia b a margin))
      ;; A torus against a ball or another torus is not asked for yet:
      ;; tires meet the road, and the road is made of boxes.
      :else [])))

(defn between
  "The manifold between two bodies, whatever shapes they are.

  The normal always runs from `a` toward `b`, so a caller never has to
  ask which way round the pair was tested.

  With no `dt` the two are taken as standing still, which is the fixed
  `speculative` margin and the right question to ask of an arrangement
  rather than a moment in a simulation."
  ([ia ib a b] (between ia ib a b 0.0))
  ([ia ib a b dt]
   (between* ia ib a b
             (when (= :box (:shape a)) (frame a))
             (when (= :box (:shape b)) (frame b))
             (margin a b dt))))

(defn- aabb
  "A world-axis box around the body, swept over the step.

  Swept, because the pair has to survive the broad phase before anything
  can ask how far apart it is: a ball that is 1.5m short of a brick at
  the top of the step and 1.5m past it at the bottom has an AABB that
  never once overlaps the brick's, and the narrow phase -- however wide
  its margin -- is never given the pair to look at.

  The sweep is one-sided, the union of where the body is and where it is
  going, rather than a symmetric pad. Same cost, half the false pairs."
  [b ^double dt]
  (let [r (case (:shape b)
            :ball (let [rr (double (:radius b))] [rr rr rr])
            (let [{:keys [half axes]} (frame b)
                  [hx hy hz] half
                  [ax ay az] axes
                  reach (fn [k] (+ (* hx (abs (nth ax k)))
                                   (* hy (abs (nth ay k)))
                                   (* hz (abs (nth az k)))))]
              [(reach 0) (reach 1) (reach 2)]))
        lo (v/sub (:pos b) r)
        hi (v/add (:pos b) r)
        d (if (rigid/inert? b) v/zero (v/scale (:vel b) dt))]
    [(mapv min lo (v/add lo d)) (mapv max hi (v/add hi d))]))

(defn near?
  "Whether two bodies' world boxes come within `gap` of each other.

  A cheap and generous test, for the solver's question of which sleeping
  bodies lean on which: nothing generates contacts between two sleeping
  bodies, so when one wakes this is how it finds the rest of its pile."
  [a b gap]
  (let [[alo ahi] (aabb a 0.0)
        [blo bhi] (aabb b 0.0)
        gap (double gap)]
    (every? (fn [k] (and (<= (double (nth alo k)) (+ (double (nth bhi k)) gap))
                         (<= (double (nth blo k)) (+ (double (nth ahi k)) gap))))
            [0 1 2])))

(defn broad-phase
  "Somewhere for `all` to keep its sorted order between steps.

  Sweep and prune earns its name across frames, not within one. Bodies
  barely move between steps, so the order it sorted them into last time
  is very nearly right this time, and restoring it costs about the number
  of bodies that actually changed places. Thrown away and rebuilt each
  step it would still beat testing every pair, and it would be paying
  `n log n` for something that is usually `n`.

  Pass the same one back every step. It resizes itself when the scene
  grows, so a demo that fires a ball into the wall needs no ceremony."
  ([] (broad-phase 256))
  ([capacity] {:sweep (sweep/sweep (max 8 (long capacity)))
               :capacity (max 8 (long capacity))}))

(defn- grow
  "The broad phase, big enough for `n` bodies."
  [broad ^long n]
  (if (and broad (>= (long (:capacity broad)) n))
    broad
    (broad-phase (max 8 (* 2 n)))))

(defn all
  "Every contact among `bodies`, as a flat vector.

  Two phases, which is the standard shape and the only way this gets
  cheap. The broad phase asks which pairs could possibly touch, by
  sorting world-axis boxes along the axis the scene is most spread along
  and comparing only the ones that overlap on it -- `allgo.spatial.sweep`.
  The narrow phase then runs the separating axis test on the survivors,
  which is where all the cost is: fifteen axes, each a handful of dot and
  cross products.

  It matters more than it looks. A wall of two hundred bricks is twenty
  thousand pairs, of which about six hundred touch; the exact test on all
  of them is most of the frame, and on the survivors it is a fraction of
  it. Two static bodies are never tested at all -- neither can move, so
  nothing they might say matters.

  With no `broad` given this allocates one per call, which is correct and
  slower; `allgo.physics.solver` keeps one in the world and hands it back
  each step.

  `dt` is the step the contacts are for, and it is what makes detection
  survive speed: both the swept boxes of the broad phase and the margin
  the narrow phase allows are how far the two bodies can travel in it.
  With no `dt` nothing moves, which is the arrangement rather than the
  step, and the fixed `speculative` margin."
  ([bodies] (all bodies nil 0.0))
  ([bodies broad] (all bodies broad 0.0))
  ([bodies broad dt]
   (let [bodies (vec bodies)
         dt (double dt)
         n (count bodies)
         broad (grow broad n)
         mins (a/f64 (* 3 n))
         maxs (a/f64 (* 3 n))
         frames (mapv #(when (= :box (:shape %)) (frame %)) bodies)]
     (dotimes [i n]
       (let [[lo hi] (aabb (nth bodies i) dt)]
         (dotimes [k 3]
           (aset ^doubles mins (+ (* 3 i) k) (double (nth lo k)))
           (aset ^doubles maxs (+ (* 3 i) k) (double (nth hi k))))))
     (into []
           (comp (mapcat (fn [[i j]]
                           (let [a (nth bodies i) b (nth bodies j)]
                             ;; Neither of them is going to move this
                             ;; step -- both static, both asleep, or one
                             ;; of each -- so there is nothing a contact
                             ;; between them could be used for. This is
                             ;; where sleeping pays: a settled wall stops
                             ;; being asked about, and the narrow phase,
                             ;; which is nearly all of detection, has
                             ;; almost nothing left to do.
                             (when-not (and (rigid/inert? a) (rigid/inert? b))
                               (between* i j a b (nth frames i) (nth frames j)
                                         (margin a b dt)))))))
           (sweep/overlapping-pairs (:sweep broad) mins maxs n)))))
