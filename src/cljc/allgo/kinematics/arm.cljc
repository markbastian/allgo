(ns allgo.kinematics.arm
  "An anthropomorphic arm: a ball at the shoulder, a hinge at the elbow, a
  ball at the wrist. Seven degrees of freedom for a hand pose that has
  only six.

  That spare degree is what makes a human arm different in kind from the
  six-axis robot of `allgo.kinematics.analytic`, not merely different in
  size. A six-axis arm reaching a given pose has eight ways to do it -- a
  handful of discrete choices. A seven-axis arm has *infinitely* many, a
  continuous family, and you can watch it: hold your hand flat on a table
  and swing your elbow. The hand does not move. The arm is solving the
  same problem a different way, the whole time.

  That motion has a name and a single number. Fix the hand and the elbow
  is confined to a circle -- it must stay one upper arm from the shoulder
  and one forearm from the wrist, and two spheres meet in a circle. Where
  on that circle the elbow sits is the **swivel angle**, and it is the
  redundancy, entire. Pick a swivel and everything else follows exactly.

  Two facts do the work, and both are worth seeing:

  The elbow's *bend* does not depend on the swivel at all. How far apart
  the shoulder and wrist are fixes it through the cosine rule, and
  swinging the elbow round changes nothing about how far apart they are.
  So the elbow bends by one amount and then rotates freely about the line
  between shoulder and wrist.

  And the shoulder is then fully determined -- not free, as three degrees
  of freedom might suggest. The upper arm must point at the elbow, which
  uses two of them, and the elbow hinge must lie across the plane the arm
  makes, which uses the third. A hinge cannot bend out of its own plane,
  and that is what spends the shoulder's last freedom.

  A configuration here is quaternions and an angle rather than a table of
  joint values, because that is what a shoulder is: no ordering of three
  hinges describes a ball joint without inventing an order that anatomy
  does not have, and gimbal lock along with it."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [clojure.math :as math]))

(defn arm
  "An arm from its two segment lengths.

  `:shoulder` places it; everything else is measured from there."
  [{:keys [upper forearm shoulder]
    :or   {upper 0.30 forearm 0.27 shoulder [0.0 0.0 0.0]}}]
  {:upper (double upper)
   :forearm (double forearm)
   :shoulder (vec shoulder)})

(defn span
  "The furthest and nearest the wrist can be from the shoulder."
  [{:keys [upper forearm]}]
  [(abs (- (double upper) (double forearm))) (+ (double upper) (double forearm))])

;; ---------------------------------------------------------------------------
;; Forward

(defn pose
  "Where a configuration puts things.

  Returns `{:shoulder :elbow :wrist :hand}` -- three points and the hand's
  orientation. `:elbow` is the angle the elbow is bent, zero being
  straight, and it bends about the upper arm frame's z axis."
  [{:keys [upper forearm] :as a} {:keys [shoulder elbow wrist]}]
  (let [s (:shoulder a)
        upper-rot shoulder
        e (v/add s (q/rotate upper-rot [(double upper) 0.0 0.0]))
        bend (q/from-axis-angle [0.0 0.0 1.0] (double elbow))
        forearm-rot (q/mul upper-rot bend)
        w (v/add e (q/rotate forearm-rot [(double forearm) 0.0 0.0]))]
    {:shoulder s
     :elbow-position e
     :elbow elbow
     :wrist w
     :hand (q/mul forearm-rot wrist)}))

;; ---------------------------------------------------------------------------
;; The elbow's circle

(defn elbow-bend
  "How far the elbow must bend for the wrist to be `distance` from the
  shoulder, or nil if it cannot be.

  The cosine rule on the triangle the arm makes, and the first of the two
  facts: this does not mention the swivel, because swinging the elbow
  about the shoulder-wrist line does not change how long that line is."
  [{:keys [upper forearm]} distance]
  (let [l1 (double upper) l2 (double forearm) d (double distance)
        c (/ (- (* d d) (* l1 l1) (* l2 l2)) (* 2.0 l1 l2))]
    (when (<= -1.0 c 1.0)
      (math/acos c))))

(defn elbow-circle
  "Where the elbow may be, given where the wrist is: a circle about the
  line from shoulder to wrist.

  Returns `{:center :radius :axis :zero :ninety}` -- `axis` the line the
  elbow swings about, and `zero` and `ninety` two perpendicular directions
  in the circle's plane, so that the elbow at swivel `psi` is
  `center + radius * (cos psi * zero + sin psi * ninety)`.

  `reference` decides where swivel zero points; the default is downward,
  which puts zero at the elbow-down posture a person rests in."
  ([a wrist-position] (elbow-circle a wrist-position [0.0 -1.0 0.0]))
  ([{:keys [upper forearm] :as a} wrist-position reference]
   (let [s (:shoulder a)
         l1 (double upper) l2 (double forearm)
         delta (v/sub wrist-position s)
         d (v/length delta)]
     (when (and (pos? d) (<= (abs (- l1 l2)) d (+ l1 l2)))
       (let [axis (v/scale delta (/ 1.0 d))
             ;; How far along the shoulder-wrist line the circle sits, and
             ;; how wide it is: the intersection of two spheres.
             a' (/ (+ (- (* d d) (* l2 l2)) (* l1 l1)) (* 2.0 d))
             r (math/sqrt (max 0.0 (- (* l1 l1) (* a' a'))))
             ;; Any direction across the axis will do for swivel zero; take
             ;; the reference, flattened into the circle's plane.
             flattened (v/sub reference (v/scale axis (v/dot reference axis)))
             zero (if (< (v/length flattened) 1e-9)
                    ;; The reference lay along the axis and says nothing.
                    (v/normalize (v/cross axis (if (< (abs (nth axis 0)) 0.9)
                                                 [1.0 0.0 0.0] [0.0 1.0 0.0])))
                    (v/normalize flattened))]
         {:center (v/add s (v/scale axis a'))
          :radius r
          :axis axis
          :zero zero
          :ninety (v/cross axis zero)})))))

(defn elbow-at
  "Where the elbow sits at a given swivel angle."
  [circle swivel]
  (let [{:keys [center radius zero ninety]} circle
        psi (double swivel)]
    (v/add center (v/scale (v/add (v/scale zero (math/cos psi))
                                  (v/scale ninety (math/sin psi)))
                           radius))))

;; ---------------------------------------------------------------------------
;; Inverse

(defn solve
  "The configuration that puts the hand at `target` with the elbow swung to
  `swivel`, or nil if the arm cannot reach.

  `target` is `{:pos :rot}` for the wrist. Every swivel gives a different
  configuration and the same hand pose, which is the whole of the
  redundancy and what the demo is for."
  ([a target swivel] (solve a target swivel [0.0 -1.0 0.0]))
  ([a target swivel reference]
   (let [s (:shoulder a)
         w (:pos target)
         circle (elbow-circle a w reference)]
     (when circle
       (let [e (elbow-at circle swivel)
             upper-dir (v/normalize (v/sub e s))
             forearm-vec (v/sub w e)
             ;; The elbow hinge lies across the plane the arm makes. This is
             ;; the second fact: it uses up the shoulder's last freedom, so
             ;; the shoulder is determined rather than free.
             hinge (v/cross upper-dir forearm-vec)]
         (when (> (v/length hinge) 1e-12)
           (let [n (v/normalize hinge)
                 ;; An orthonormal frame: the upper arm, the hinge, and the
                 ;; direction the forearm swings toward.
                 y (v/cross n upper-dir)
                 shoulder-rot (q/from-matrix
                               [[(nth upper-dir 0) (nth y 0) (nth n 0)]
                                [(nth upper-dir 1) (nth y 1) (nth n 1)]
                                [(nth upper-dir 2) (nth y 2) (nth n 2)]])
                 bend (elbow-bend a (v/distance w s))]
             (when bend
               (let [forearm-rot (q/mul shoulder-rot
                                        (q/from-axis-angle [0.0 0.0 1.0] bend))]
                 {:shoulder shoulder-rot
                  :elbow bend
                  ;; Whatever turn the hand still needs, the wrist makes.
                  :wrist (q/mul (q/inverse forearm-rot) (:rot target))})))))))))

(defn swivel-of
  "Which swivel angle a configuration is at.

  The inverse of `solve`'s parameter, for recovering where an arm already
  is before moving it somewhere else."
  ([a configuration] (swivel-of a configuration [0.0 -1.0 0.0]))
  ([a configuration reference]
   (let [{:keys [elbow-position wrist]} (pose a configuration)
         circle (elbow-circle a wrist reference)]
     (when circle
       (let [offset (v/sub elbow-position (:center circle))]
         (math/atan2 (v/dot offset (:ninety circle))
                     (v/dot offset (:zero circle))))))))

(defn swivel-samples
  "`n` configurations evenly spaced around the elbow's circle.

  The family of answers, made countable so it can be drawn. They are not
  different solutions in the sense the six-axis arm has different
  solutions -- there is no gap between them, and any angle in between is
  just as good."
  ([a target n] (swivel-samples a target n [0.0 -1.0 0.0]))
  ([a target n reference]
   (into []
         (keep (fn [i]
                 (let [psi (* 2.0 math/PI (/ (double i) (double n)))]
                   (when-let [c (solve a target psi reference)]
                     (assoc c :swivel psi)))))
         (range n))))
