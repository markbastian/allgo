(ns allgo.physics.toi
  "When two moving bodies will first touch, by conservative advancement.

  `allgo.physics.contact` answers where two bodies touch *now*, and the
  speculative margin lets a solver see a touch coming within the step.
  That is enough for the impulse solvers, which are told about a gap and
  refuse to close it faster than it is wide. It is not enough for XPBD,
  which integrates a substep and only then pushes overlaps apart: a gap
  it has not yet reached buys it nothing, so above about 120 meters a
  second it steps over a thin slab and comes out the other side.

  ## The algorithm

  Conservative advancement, after Mirtich. Take the distance between the
  two shapes and the fastest either could be closing it; their quotient
  is a length of time in which nothing can possibly happen, so advance
  by it. The distance is recomputed, it is smaller, the safe step is
  shorter, and the sequence walks down onto the moment of first contact
  without ever stepping past it.

  Two things make the bound safe rather than merely plausible. The
  closing speed is measured along the line between the nearest points,
  which is the direction that matters and the only one; and rotation is
  allowed for by `|omega| * r`, the fastest any point of a body can move
  about its own center, rather than by turning the shapes and asking
  again. That over-states how fast a spinning body closes a gap, which
  is the right way to be wrong -- the answer comes back early, never
  late.

  Distance comes from `allgo.geometry.gjk`, which asks a shape only for
  its support mapping, so a box and a ball are the same problem.

  ## What it is for

  A step that advances to the time of impact rather than through it.
  `allgo.physics.solver`'s XPBD path uses it to choose how long a
  substep may be, which is the difference between stopping a bullet and
  reporting that it was never there."
  (:require [allgo.geometry.gjk :as gjk]
            [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.rigid :as rigid]))

(defn- signs [[x y z]]
  [(if (neg? (double x)) -1.0 1.0)
   (if (neg? (double y)) -1.0 1.0)
   (if (neg? (double z)) -1.0 1.0)])

(defn support
  "A body's support mapping, with its center moved to `at`.

  Only the position is moved and not the orientation, because the
  rotation is accounted for in the closing-speed bound instead. Turning
  the shape for every trial time would cost a support mapping rebuild
  apiece and buy a tighter answer than the bound needs."
  [b at]
  (if (= :ball (:shape b))
    (let [r (double (:radius b))]
      (fn [dir] (v/add at (v/scale (v/normalize dir) r))))
    (let [half (mapv #(* 0.5 (double %)) (:size b))
          rot (:rot b)
          inv (:inv-rot b)]
      (fn [dir]
        (v/add at (q/rotate rot (v/mul half (signs (q/rotate inv dir)))))))))

(defn reach
  "The furthest any point of a body is from its center.

  What `|omega| * r` needs: the most a spin of one radian a second can
  move a part of this body."
  ^double [b]
  (if (= :ball (:shape b))
    (double (:radius b))
    (* 0.5 (v/length (:size b)))))

(def ^:private tolerance
  "How close counts as touching.

  Conservative advancement converges on contact and never quite arrives
  -- each step is a fraction of the remaining distance -- so it has to
  be told when to stop. A quarter of the solver's slop, which is close
  enough that the step which follows has something to solve."
  0.00125)

(def ^:private max-iterations 24)

(defn impact
  "When `a` and `b` first touch within `dt`, and how fast they are
  closing when they do.

  `{:t :closing}`, or nil if they do not touch: either they are moving
  apart or they do not close the distance in time. A `:t` of zero means
  they are already touching, which is the ordinary case and not a
  bullet.

  `:closing` is what a caller needs to ask for a *little past* the
  moment of impact. Arriving exactly on it is not always useful -- a
  position solver has nothing to solve at zero overlap -- and dividing a
  wanted depth by this gives the extra time that buys it."
  [a b ^double dt]
  (let [va (:vel a) vb (:vel b)
        ra (reach a) rb (reach b)
        ;; The most either can add to the closing speed by spinning.
        spin (+ (* (v/length (:omega a)) ra) (* (v/length (:omega b)) rb))]
    (loop [t 0.0 i 0 last-closing 0.0]
      (if (>= i max-iterations)
        {:t t :closing last-closing}
        (let [{:keys [distance direction]}
              (gjk/distance (support a (v/add-scaled (:pos a) va t))
                            (support b (v/add-scaled (:pos b) vb t)))
              d (double distance)]
          (if (<= d tolerance)
            {:t t :closing last-closing}
            (let [n (v/normalize direction)
                  ;; `direction` runs from A's nearest point to B's, so
                  ;; A closes the gap by moving along it and B by moving
                  ;; against it.
                  closing (+ (- (v/dot va n) (v/dot vb n)) spin)]
              (if (<= closing 1e-9)
                nil
                (let [t' (+ t (/ d closing))]
                  (if (> t' dt) nil (recur t' (inc i) closing)))))))))))

(defn time-of-impact
  "The first fraction of `dt` at which `a` and `b` could touch, or nil."
  [a b ^double dt]
  (:t (impact a b dt)))

(defn earliest
  "The soonest any pair of `bodies` could touch within `dt`, or nil.

  Pairs come from a caller that has already narrowed them down --
  `allgo.physics.contact` has the broad phase, and a pair it did not
  offer is a pair whose swept boxes do not meet."
  [bodies pairs ^double dt]
  (reduce (fn [best [i j]]
            (let [a (nth bodies i) b (nth bodies j)]
              (if (and (rigid/inert? a) (rigid/inert? b))
                best
                (if-let [t (time-of-impact a b dt)]
                  (if (or (nil? best) (< (double t) (double best))) t best)
                  best))))
          nil
          pairs))
