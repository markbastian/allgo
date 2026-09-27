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

(defn- radius-on
  "How far the box reaches from its center along `axis`."
  ^double [{:keys [half axes]} axis]
  (let [[hx hy hz] half
        [ax ay az] axes]
    (+ (* hx (abs (v/dot ax axis)))
       (* hy (abs (v/dot ay axis)))
       (* hz (abs (v/dot az axis))))))

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

(defn- overlap-on
  "How much the two boxes overlap along `axis`.

  Three answers, and telling them apart matters. `[o n]` is an overlap of
  `o` along the unit axis `n`; `:separated` is a gap, which proves the
  boxes miss and ends the search; `:degenerate` is an axis that is not an
  axis at all.

  The last is the case that bites. The cross product of two parallel edge
  directions is the zero vector, and two boxes sitting square to each
  other have three such pairs. Reading those as separating -- which is
  what a single nil answer invites -- reports every axis-aligned box in
  the world as touching nothing."
  [a b axis ^double margin]
  (let [len (v/length axis)]
    (if (<= len 1e-6)
      :degenerate
      (let [n (v/scale axis (/ 1.0 len))
            d (abs (v/dot (v/sub (:pos b) (:pos a)) n))
            o (- (+ (radius-on a n) (radius-on b n)) d)]
        ;; Within the margin counts as touching. `o` comes back negative
        ;; there and stays negative all the way to the solver, which is
        ;; what tells it this is a gap, not an overlap.
        (if (pos? (+ o margin)) [o n] :separated)))))

(defn- face-verts
  "The four corners of the box face whose outward normal is `axis-index`
  in direction `sign`, counter-clockwise about that normal."
  [{:keys [pos half axes]} axis-index sign]
  (let [h half
        as axes
        n (v/scale (nth as axis-index) (double sign))
        i (mod (inc (long axis-index)) 3)
        j (mod (inc (long i)) 3)
        u (v/scale (nth as i) (nth h i))
        w (v/scale (nth as j) (nth h j))
        c (v/add pos (v/scale n (nth h axis-index)))]
    [(v/sub (v/sub c u) w)
     (v/sub (v/add c u) w)
     (v/add (v/add c u) w)
     (v/add (v/sub c u) w)]))

(defn- most-opposed-face
  "The face of `b` pointing most directly back along `n`.

  `n` runs from the reference box toward this one, so the face we want
  is the one whose outward normal is closest to `-n` -- the face being
  pressed into."
  [{:keys [axes]} n]
  (let [as axes]
    (first
     (sort-by second
              (for [i (range 3) s [1.0 -1.0]]
                [[i s] (v/dot (v/scale (nth as i) s) n)])))))

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

(defn- clip-against-plane
  "Sutherland-Hodgman: the part of `poly` on the inner side of the plane
  `(dot n x) <= d`.

  `poly` is a sequence of `[point id]`, and the ids are the point of it.
  A corner that survives the clip keeps the id it came in with; a point
  made by cutting an edge gets one built from the plane that cut it and
  the corner the cut started from. So the same configuration next step
  produces the same ids, and a contact can be recognized as the one that
  was there before rather than merely as one nearby."
  [poly n d tag]
  (let [pts (vec poly)
        cnt (count pts)]
    (into []
          (mapcat (fn [i]
                    (let [[p pid] (nth pts i)
                          [q' qid] (nth pts (mod (inc i) cnt))
                          dp (- (v/dot n p) d)
                          dq (- (v/dot n q') d)
                          cut (fn [] [(v/lerp p q' (/ dp (- dp dq))) [tag pid]])]
                      (cond
                        (and (<= dp clip-eps) (<= dq clip-eps)) [[q' qid]]
                        (<= dp clip-eps) [(cut)]
                        (<= dq clip-eps) [(cut) [q' qid]]
                        :else []))))
          (range cnt))))

(defn- distinct-points
  "A transducer over `[point id]` pairs dropping any point already seen.

  Within `tol` counts as already seen. The comparison is against every
  point kept so far, which is fine at four of them and would not be at
  four hundred."
  [tol]
  (let [tol (double tol)]
    (fn [rf]
      (let [seen (volatile! [])]
        (fn
          ([] (rf))
          ([acc] (rf acc))
          ([acc [p _ :as x]]
           (if (some #(< (v/distance p %) tol) @seen)
             acc
             (do (vswap! seen conj p)
                 (rf acc x)))))))))

(defn- face-contacts
  "Contacts from a face of `ref-body` pressing into `inc-body`.

  The reference face's four side planes bound the region the contact can
  lie in; clipping the opposing face to them and keeping what is still
  below the reference plane gives the resting polygon."
  [ia ib ref-body inc-body ref-axis ref-sign n flip? margin]
  (let [margin (double margin)
        h (:half ref-body)
        as (:axes ref-body)
        ref-n (v/scale (nth as ref-axis) (double ref-sign))
        ref-c (v/add (:pos ref-body) (v/scale ref-n (nth h ref-axis)))
        ref-d (v/dot ref-n ref-c)
        [[inc-axis inc-sign] _] (most-opposed-face inc-body ref-n)
        poly (map-indexed (fn [i p] [p i]) (face-verts inc-body inc-axis inc-sign))
        ;; The four sides of the reference face, as planes facing outward,
        ;; each tagged so that a point cut by one can say which.
        sides (map-indexed
               (fn [t [sn sd]] [sn sd t])
               (for [k (range 3) :when (not= k ref-axis) s [1.0 -1.0]]
                 [(v/scale (nth as k) s)
                  (v/dot (v/scale (nth as k) s)
                         (v/add (:pos ref-body) (v/scale (nth as k) (* s (nth h k)))))]))
        clipped (reduce (fn [p [sn sd t]] (if (seq p) (clip-against-plane p sn sd t) p))
                        poly
                        sides)]
    (into []
          ;; A face clipped against a face it exactly coincides with --
          ;; two identical bricks, square on -- has every corner sitting
          ;; on a clip plane, and Sutherland-Hodgman answers each of those
          ;; twice. Keeping both gives a pair six or eight contact points
          ;; where it has four, and which of them survive changes from
          ;; step to step. One point per place, first id wins, and the
          ;; ordering is the incident face's own, so the choice is the
          ;; same every step.
          (comp (distinct-points 1e-6)
                (keep (fn [[p pid]]
                        (let [depth (- ref-d (v/dot ref-n p))]
                          (when (>= depth (- margin))
                      ;; Reported on the reference surface rather than at
                      ;; the clipped point, so both bodies agree where the
                      ;; touch is.
                            {:a ia :b ib
                             :point (v/add-scaled p ref-n (* 0.5 depth))
                             :normal (if flip? (v/negate n) n)
                       ;; Which corner of which face, pressed into which
                       ;; face -- the same answer every step for as long
                       ;; as the two bodies stay in the same arrangement,
                       ;; which is what warm starting needs to recognize
                       ;; it by. `flip?` is in it because it says which
                       ;; body was the reference, and that decides what
                       ;; the rest of the id means.
                             :id [(boolean flip?) ref-axis ref-sign inc-axis inc-sign pid]
                             ;; Signed: positive is an overlap, negative
                             ;; a gap the solver is allowed to see coming.
                             :depth depth})))))
          clipped)))

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

(defn box-box
  "The manifold between two boxes, or an empty vector.

  `a` and `b` are `frame`s, not bodies."
  [ia ib a b margin]
  (let [margin (double margin)
        as (:axes a)
        bs (:axes b)
        face-tests (concat (map-indexed (fn [i ax] [:a i nil ax]) as)
                           (map-indexed (fn [j ax] [:b nil j ax]) bs))
        edge-tests (for [i (range 3) j (range 3)]
                     [:edge i j (v/cross (nth as i) (nth bs j))])
        ;; Faces get a small handicap in their favor, and a's faces one
        ;; over b's. Equal overlaps otherwise flip between a face result
        ;; and an edge result, or between a's face and b's, on
        ;; floating-point noise -- and a resting brick loses three of its
        ;; four contact points for a frame, or keeps four points that are
        ;; named after the other box and so match nothing from last step.
        ;; Two identical bricks stacked square on each other tie exactly,
        ;; which is the commonest case there is.
        best (reduce (fn [best [kind i j ax]]
                       (let [r (overlap-on a b ax margin)]
                         (cond
                           (= r :separated) (reduced nil)
                           (= r :degenerate) best
                           :else
                           (let [[o n] r
                                 score (if (= kind :edge) (* o 1.02) o)]
                             (if (or (nil? best) (< score (:score best)))
                               {:score score :kind kind :i i :j j :normal n :depth o}
                               best)))))
                     nil
                     (concat face-tests edge-tests))]
    (if (nil? best)
      []
      (let [{:keys [kind i j normal depth]} best
            ;; Point the axis from a toward b.
            n (if (neg? (v/dot normal (v/sub (:pos b) (:pos a))))
                (v/negate normal)
                normal)]
        (case kind
          :a (face-contacts ia ib a b i
                            (if (pos? (v/dot (nth as i) n)) 1.0 -1.0) n false margin)
          :b (face-contacts ia ib b a j
                            (if (pos? (v/dot (nth bs j) (v/negate n))) 1.0 -1.0)
                            (v/negate n) true margin)
          :edge (edge-contact ia ib a b n depth i j))))))

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
      (let [phi (/ (- (math/sqrt 5.0) 1.0) 2.0)
            t (loop [lo (- (double best) step) hi (+ (double best) step) k 0]
                (if (= k 24)
                  (* 0.5 (+ lo hi))
                  (let [m1 (- hi (* phi (- hi lo)))
                        m2 (+ lo (* phi (- hi lo)))]
                    (if (> (depth-at m1) (depth-at m2))
                      (recur lo m2 (inc k))
                      (recur m1 hi (inc k))))))
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
