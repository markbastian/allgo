(ns allgo.kinematics.chain
  "Serial kinematic chains, after John J. Craig's Introduction to Robotics.

  A robot arm is a sequence of rigid links joined end to end, and the only
  hard part of describing one is agreeing where each link's frame sits.
  Denavit and Hartenberg's answer is to allow just four numbers per link,
  chosen so that the frames are forced rather than picked:

    alpha  the twist: how far the previous joint axis is rotated about the
           common perpendicular to reach this one
    a      the link length: how far along that perpendicular
    d      the offset: how far along this joint's own axis
    theta  the joint angle: how far around it

  Four numbers cannot describe an arbitrary transform, which takes six --
  and that is the point. The missing two are exactly what is fixed by
  insisting the x axis lie along the common perpendicular. A chain
  described this way has no freedom left to describe it wrongly.

  ## Which convention

  Craig's *modified* parameters, where a link's frame sits at its
  proximal end and the transform is

    Rot_x(alpha) . Trans_x(a) . Rot_z(theta) . Trans_z(d)

  with alpha and a belonging to the *previous* link. The other convention
  in common use -- standard DH, frame at the distal end -- has the same
  four names in a different order and a different transform. The two are
  not interchangeable, and a table read under the wrong one gives an arm
  that bends plausibly and reaches the wrong place. Tables here say which
  they are; there is no way to detect it from the numbers.

  A transform is `{:rot quaternion :pos vector}` rather than a matrix,
  because that is what the rest of this library speaks and what a renderer
  wants. It is checked against Craig's matrix in the tests."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [clojure.math :as math]))

(def identity-transform {:rot q/identity-q :pos v/zero})

(defn compose
  "`a` followed by `b`, read as `a` then applying `b` in `a`'s frame."
  [a b]
  {:rot (q/mul (:rot a) (:rot b))
   :pos (v/add (:pos a) (q/rotate (:rot a) (:pos b)))})

(defn invert [{:keys [rot pos]}]
  (let [r (q/inverse rot)]
    {:rot r :pos (v/negate (q/rotate r pos))}))

(defn link
  "One link. `:joint` is `:revolute` (theta varies) or `:prismatic` (d
  varies); `:limits` is `[lo hi]` on whichever that is, or nil for free."
  [{:keys [alpha a d theta joint limits]
    :or   {alpha 0.0 a 0.0 d 0.0 theta 0.0 joint :revolute}}]
  {:alpha (double alpha) :a (double a) :d (double d) :theta (double theta)
   :joint joint
   :limits (when limits (mapv double limits))})

(defn chain
  "A chain from a sequence of link specifications.

  `:tool` is a fixed transform from the last joint frame to whatever is
  actually doing the work -- a gripper, a fingertip. Kept separate because
  changing the tool must not mean re-deriving the arm."
  ([links] (chain links identity-transform))
  ([links tool]
   {:links (mapv link links) :tool tool}))

(defn joint-count [{:keys [links]}] (count links))

(defn link-transform
  "The transform across one link, given that joint's value.

  A revolute joint adds its value to `theta`, a prismatic one to `d`, so
  the value is an offset from the link's rest pose rather than a
  replacement -- which is what lets a table carry a non-zero home
  position."
  [{:keys [alpha a d theta joint]} value]
  (let [value (double value)
        theta (if (= joint :prismatic) theta (+ theta value))
        d     (if (= joint :prismatic) (+ d value) d)]
    {:rot (q/mul (q/from-axis-angle [1.0 0.0 0.0] alpha)
                 (q/from-axis-angle [0.0 0.0 1.0] theta))
     ;; The translation column of Craig's matrix: along x by the link
     ;; length, then along the twisted z by the offset.
     :pos [a (* (- d) (math/sin alpha)) (* d (math/cos alpha))]}))

(defn frames
  "Every joint frame in base coordinates, starting with the base itself.

  `n + 1` of them for `n` joints, so `frames[i]` is the frame of joint `i`
  and its z axis is that joint's axis -- which is what the Jacobian reads."
  [{:keys [links]} values]
  (reduce (fn [acc [lnk value]]
            (conj acc (compose (peek acc) (link-transform lnk value))))
          [identity-transform]
          (map vector links values)))

(defn pose
  "Where the tool ends up."
  [{:keys [tool] :as c} values]
  (compose (peek (frames c values)) tool))

(defn joint-axis
  "The axis of joint `i`, in base coordinates."
  [frames i]
  (q/rotate (:rot (nth frames i)) [0.0 0.0 1.0]))

(defn jacobian
  "The geometric Jacobian: how the tool moves per unit of each joint.

  Six rows -- three of linear velocity, three of angular -- and a column
  per joint. A revolute joint sweeps the tool around its axis, so its
  linear contribution is the axis crossed with the arm from that joint out
  to the tool, and its angular contribution is the axis itself. A
  prismatic joint slides, so it contributes its axis linearly and nothing
  angularly.

  This is what a numerical solver descends and what says how close the arm
  is to a singularity -- a configuration where some direction of motion
  costs infinite joint speed, which shows up here as the matrix losing
  rank."
  [{:keys [links] :as c} values]
  (let [fs (frames c values)
        tool-pos (:pos (pose c values))]
    (mapv (fn [i]
            (let [axis (joint-axis fs (inc i))
                  origin (:pos (nth fs (inc i)))]
              (if (= :prismatic (:joint (nth links i)))
                (into (vec axis) v/zero)
                (into (vec (v/cross axis (v/sub tool-pos origin))) axis))))
          (range (count links)))))

;; ---------------------------------------------------------------------------
;; Joint limits

(defn limits [{:keys [links]}] (mapv :limits links))

(defn within-limits?
  "Whether every joint is inside its range. Joints with no range are free."
  [c values]
  (every? (fn [[[lo hi] value]]
            (or (nil? lo) (<= (double lo) (double value) (double hi))))
          (map vector (limits c) values)))

(defn clamp
  "The nearest configuration that respects the limits."
  [c values]
  (mapv (fn [lim value]
          (if lim
            (min (double (second lim)) (max (double (first lim)) (double value)))
            value))
        (limits c)
        values))

(defn wrap-angle
  "An angle brought into `(-pi, pi]`.

  Joint solutions come out of inverse trigonometry anywhere on the circle,
  and two that differ by a full turn are the same arm in the same place --
  but only one of them is likely to be inside the joint's range."
  ^double [theta]
  (let [two-pi (* 2.0 math/PI)
        r (rem (double theta) two-pi)]
    (cond (> r math/PI) (- r two-pi)
          (<= r (- math/PI)) (+ r two-pi)
          :else r)))

(defn home
  "Every joint at zero."
  [c] (vec (repeat (joint-count c) 0.0)))
