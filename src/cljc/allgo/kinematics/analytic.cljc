(ns allgo.kinematics.analytic
  "Closed-form inverse kinematics for a six-axis arm with a spherical
  wrist, after Craig chapter 4.

  Inverse kinematics is harder than forward: forward is one matrix
  product, inverse is a system of transcendental equations that in general
  has no closed form at all. Pieper's result is that it does have one
  whenever three consecutive joint axes intersect at a point -- which is
  why almost every six-axis industrial robot is built with its last three
  axes meeting at the wrist. The geometry is chosen to make the algebra
  possible.

  What that buys is *kinematic decoupling*. The wrist centre cannot be
  moved by the wrist, only by the first three joints, so:

    1. find the wrist centre by stepping back from the tool along its own
       approach axis
    2. solve joints 1 to 3 for the position of that point -- pure geometry,
       a triangle
    3. whatever orientation those three left behind, the wrist makes up
       the difference; solve joints 4 to 6 for the rotation that remains

  Each stage has a choice, and the choices multiply: the shoulder can
  reach left or right of the target, the elbow can be up or down, and the
  wrist can arrive flipped. Two by two by two is **eight** configurations
  that put the tool in exactly the same place -- the reason a robot
  program specifies a configuration and not just a pose.

  Fewer than eight come back when the arm cannot do it: a target beyond
  reach gives none, one exactly at the limit of reach gives a pair that
  have collapsed together, and a wrist singularity -- joint 5 straight,
  where joints 4 and 6 turn about the same line -- gives a continuous
  family that is reported as a single representative.

  Every candidate is checked by running it back through forward kinematics
  before it is returned. The algebra below has eight branches and a sign
  error in one of them would otherwise be an arm that goes to the wrong
  place in one configuration out of eight, which is exactly the sort of
  thing that survives a demo. Checking costs one forward pass per
  candidate and makes the function honest: what it returns, it has
  verified."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.kinematics.chain :as k]
            [clojure.math :as math]))

(def ^:private half-pi (/ math/PI 2.0))

(defn- near? [^double a ^double b] (< (abs (- a b)) 1e-9))

(defn spherical-wrist?
  "Whether the last three axes meet at a point, which is what makes a
  closed form possible.

  Read off the table rather than measured: the wrist joints must add no
  offset of their own, and the two twists between them must be right
  angles that cancel."
  [{:keys [links]}]
  (and (= 6 (count links))
       (let [[_ l1 _ _ l4 l5] links]
         (and (near? (:a (nth links 4)) 0.0) (near? (:d (nth links 4)) 0.0)
              (near? (:a (nth links 5)) 0.0) (near? (:d (nth links 5)) 0.0)
              (near? (abs (:alpha l4)) half-pi)
              (near? (abs (:alpha l5)) half-pi)
              (near? (+ (:alpha l4) (:alpha l5)) 0.0)
              (near? (:alpha l1) (- half-pi))))))

(defn- geometry
  "The handful of table entries the algebra needs, under the names Craig
  gives them."
  [{:keys [links]}]
  (let [[l0 l1 l2 l3 l4 _] links]
    {:d1 (:d l0) :d2 (:d l1) :a2 (:a l2) :d3 (:d l2)
     :a3 (:a l3) :d4 (:d l3) :alpha3 (:alpha l3) :alpha4 (:alpha l4)
     :rest (mapv :theta links)}))

(defn- solve-shoulder
  "Joint 1, twice: the arm can pass either side of the target.

  Seen from above the wrist centre sits at a fixed offset `b` from the
  shoulder's turning axis, so the base angle is where the target lies
  minus how far that offset carries it round."
  [px py b]
  (let [r2 (- (+ (* px px) (* py py)) (* b b))]
    (when-not (neg? r2)
      (let [r (math/sqrt r2)]
        ;; The two roots of the offset triangle: reaching forwards, and
        ;; reaching backwards over the shoulder.
        [[(- (math/atan2 py px) (math/atan2 b r)) r]
         [(- (math/atan2 py px) (math/atan2 b (- r))) (- r)]]))))

(defn- solve-elbow
  "Joint 3, twice: elbow up and elbow down.

  The upper arm and the forearm make a triangle with the line from
  shoulder to wrist, and the cosine rule gives the included angle. The
  forearm's own offset means the angle is measured from a line that is not
  quite the forearm, which is what the second arctangent takes out."
  [reach lift {:keys [a2 a3 d4]}]
  (let [rr (math/sqrt (+ (* a3 a3) (* d4 d4)))
        k  (/ (- (+ (* reach reach) (* lift lift)) (* a2 a2) (* a3 a3) (* d4 d4))
              (* 2.0 a2))
        c  (/ k rr)]
    (when (<= (abs c) 1.0)
      (let [phi (math/atan2 d4 a3)
            base (math/acos c)]
        [(- base phi) (- (- base) phi)]))))

(defn- solve-upper-arm
  "Joint 2, once the elbow is known: where the upper arm has to point for
  the forearm to finish at the wrist centre."
  [reach lift theta3 {:keys [a2 a3 d4]}]
  (let [c3 (math/cos theta3) s3 (math/sin theta3)
        wx (+ a2 (* a3 c3) (- (* d4 s3)))
        wy (+ (* a3 s3) (* d4 c3))]
    (- (math/atan2 (- lift) reach) (math/atan2 wy wx))))

(defn- solve-wrist
  "Joints 4, 5 and 6 from the rotation the arm has left to make.

  With the two wrist twists at right angles and cancelling, the three
  rotations compose into a plain Z-Y-Z Euler sequence, and that has the
  textbook two solutions -- one with the middle angle positive, one with
  it negative and the outer two turned half a circle."
  [m alpha4]
  (let [[[_ _ m02] [_ _ m12] [m20 m21 m22]] m
        cb (min 1.0 (max -1.0 m22))
        sb (math/sqrt (max 0.0 (- 1.0 (* cb cb))))
        ;; Which way joint 5's axis points depends on the sign of the
        ;; twist ahead of it.
        flip (- (math/sin alpha4))]
    (if (< sb 1e-7)
      ;; Joint 5 straight: joints 4 and 6 turn about the same line and only
      ;; their sum is determined. One representative, with joint 4 left at
      ;; zero.
      [{:theta4 0.0
        :theta5 (* flip (if (pos? cb) 0.0 math/PI))
        :theta6 (math/atan2 m21 (- m20))
        :singular? true}]
      (for [s [1.0 -1.0]]
        (let [sb (* s sb)
              b  (math/atan2 sb cb)]
          {:theta4 (math/atan2 (/ m12 sb) (/ m02 sb))
           :theta5 (* flip b)
           :theta6 (math/atan2 (/ m21 sb) (/ (- m20) sb))
           :singular? false})))))

(def ^:private pose-tolerance
  "How near a candidate has to land to count as a solution. Loose enough
  for the rounding of a dozen trigonometric calls, tight enough that a
  wrong branch cannot hide."
  1e-6)

(defn- reaches?
  [chain target values]
  (let [p (k/pose chain values)]
    (and (< (v/distance (:pos p) (:pos target)) pose-tolerance)
         ;; A quaternion and its negation are the same orientation.
         (< (min (q/angle (q/between (:rot p) (:rot target)))
                 (q/angle (q/between (:rot p) (mapv - (:rot target)))))
            1e-5))))

(defn solutions
  "Every configuration that puts the tool at `target`, as joint values.

  `target` is `{:rot :pos}` in base coordinates. Up to eight come back,
  ordered shoulder then elbow then wrist. Each has been checked by forward
  kinematics; anything the algebra produced that does not actually reach
  the target is not returned.

  Pass `:limits? true` to keep only configurations the joints can
  actually take."
  ([chain target] (solutions chain target {}))
  ([chain target {:keys [limits?]}]
   (when-not (spherical-wrist? chain)
     (throw (ex-info "closed-form inverse kinematics needs a spherical wrist: six joints, the last three meeting at a point"
                     {:joints (k/joint-count chain)})))
   (let [{:keys [d1 d2 d3 alpha3 alpha4 rest] :as g} (geometry chain)
         ;; Step back from the tool to the last joint frame, whose origin
         ;; is the wrist centre -- the wrist joints add no offset, so they
         ;; cannot move it.
         f6 (k/compose target (k/invert (:tool chain)))
         [px py pz] (:pos f6)
         b (+ d2 d3)
         found
         (for [[theta1 reach] (or (solve-shoulder px py b) [])
               :let [lift (- pz d1)]
               theta3 (or (solve-elbow reach lift g) [])
               :let [theta2 (solve-upper-arm reach lift theta3 g)
                     ;; Whatever rotation the arm has not made, the wrist
                     ;; must. Taking out the twist ahead of joint 4 leaves
                     ;; a plain Euler sequence.
                     r3 (:rot (peek (k/frames chain [(- theta1 (rest 0))
                                                     (- theta2 (rest 1))
                                                     (- theta3 (rest 2))])))
                     m (q/to-matrix (q/mul (q/from-axis-angle [1.0 0.0 0.0] (- alpha3))
                                           (q/mul (q/inverse r3) (:rot f6))))]
               {:keys [theta4 theta5 theta6]} (solve-wrist m alpha4)]
           (mapv k/wrap-angle
                 [(- theta1 (rest 0)) (- theta2 (rest 1)) (- theta3 (rest 2))
                  (- theta4 (rest 3)) (- theta5 (rest 4)) (- theta6 (rest 5))]))]
     (cond->> (vec (distinct (filter #(reaches? chain target %) found)))
       limits? (filterv #(k/within-limits? chain %))))))

(defn nearest
  "The solution that moves the joints least from where they are.

  What a robot controller picks, and for the reason it picks it: the eight
  configurations put the tool in the same place but are not
  interchangeable on the way there. Swapping between two of them means
  swinging the whole arm through a large reconfiguration, which takes time
  and sweeps a volume that may not be empty."
  ([chain target from] (nearest chain target from {}))
  ([chain target from opts]
   (let [candidates (solutions chain target opts)]
     (when (seq candidates)
       (apply min-key
              (fn [s] (reduce + (map (fn [a b] (abs (k/wrap-angle (- a b)))) s from)))
              candidates)))))
