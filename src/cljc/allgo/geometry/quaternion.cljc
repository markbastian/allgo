(ns allgo.geometry.quaternion
  "Unit quaternions as `[x y z w]`, the vector part first.

  A quaternion is how an orientation is stored when it is going to be
  integrated. The alternatives are worse in specific ways: three Euler
  angles lose a degree of freedom when two axes line up, and a rotation
  matrix drifts out of being a rotation as error accumulates and is
  awkward to pull back -- nine numbers constrained to three. A quaternion
  has one constraint, unit length, and restoring it is a divide.

  The order is `[x y z w]` rather than `[w x y z]` because that is what
  three.js uses, and every one of these ends up in a `THREE.Quaternion`.
  Both conventions are common; mixing them silently gives rotations that
  look almost right."
  (:require [allgo.geometry.vec3 :as v]
            [clojure.math :as math]))

(def identity-q
  "No rotation at all."
  [0.0 0.0 0.0 1.0])

(defn mul
  "`a` then `b` read right to left: the rotation `b` followed by `a`,
  as with matrices."
  [[ax ay az aw] [bx by bz bw]]
  [(+ (* aw bx) (* ax bw) (* ay bz) (- (* az by)))
   (+ (* aw by) (- (* ax bz)) (* ay bw) (* az bx))
   (+ (* aw bz) (* ax by) (- (* ay bx)) (* az bw))
   (- (* aw bw) (* ax bx) (* ay by) (* az bz))])

(defn conjugate [[x y z w]] [(- x) (- y) (- z) w])

(defn dot ^double [[ax ay az aw] [bx by bz bw]]
  (+ (* ax bx) (* ay by) (* az bz) (* aw bw)))

(defn length ^double [q] (math/sqrt (dot q q)))

(defn normalize
  "Back onto the unit sphere, which integration walks off."
  [[x y z w :as q]]
  (let [l (length q)]
    (if (zero? l) identity-q [(/ x l) (/ y l) (/ z l) (/ w l)])))

(defn inverse
  "For a unit quaternion this is the conjugate; the general form divides
  by the squared length, and doing it properly costs one divide and means
  a quaternion that has drifted still inverts to something usable."
  [[x y z w :as q]]
  (let [n (dot q q)]
    (if (zero? n) identity-q [(/ (- x) n) (/ (- y) n) (/ (- z) n) (/ w n)])))

(defn rotate
  "Turns a vector by a quaternion.

  The long way is `q * (v, 0) * q'`, two quaternion products. This is the
  same thing rearranged into two cross products, which is the form worth
  having because it runs on every vertex of everything."
  [[qx qy qz qw] v]
  (let [u [qx qy qz]
        t (v/scale (v/cross u v) 2.0)]
    (v/add v (v/scale t qw) (v/cross u t))))

(defn from-axis-angle
  "A rotation of `angle` radians about `axis`."
  [axis angle]
  (let [[x y z] (v/normalize axis)
        h (* 0.5 (double angle))
        s (math/sin h)]
    [(* x s) (* y s) (* z s) (math/cos h)]))

(defn from-euler
  "From Euler angles in three.js's default XYZ order.

  Which composes as `Rx * Ry * Rz`. Both readings of \"XYZ\" are in use
  and they are opposites, so to be unambiguous: the object turns about its
  own x, then its own y, then its own z; equivalently, acting on a vector,
  the z rotation is applied first. A test pins it, because getting this
  backwards gives rotations that are wrong only when two of the angles are
  nonzero, which is the kind of thing that survives a demo."
  [x y z]
  (let [hx (* 0.5 (double x)) hy (* 0.5 (double y)) hz (* 0.5 (double z))
        c1 (math/cos hx) c2 (math/cos hy) c3 (math/cos hz)
        s1 (math/sin hx) s2 (math/sin hy) s3 (math/sin hz)]
    [(+ (* s1 c2 c3) (* c1 s2 s3))
     (- (* c1 s2 c3) (* s1 c2 s3))
     (+ (* c1 c2 s3) (* s1 s2 c3))
     (- (* c1 c2 c3) (* s1 s2 s3))]))

(defn to-axis-angle
  "`[axis angle]`, the inverse of `from-axis-angle`."
  [q]
  (let [[x y z w] (normalize q)
        s (math/sqrt (max 0.0 (- 1.0 (* w w))))]
    (if (< s 1e-9)
      [[1.0 0.0 0.0] 0.0]
      [[(/ x s) (/ y s) (/ z s)] (* 2.0 (math/acos (min 1.0 (max -1.0 w))))])))

(defn angle
  "How far this quaternion turns, in radians, always the short way round."
  ^double [q]
  (let [[_ _ _ w] (normalize q)]
    (* 2.0 (math/acos (min 1.0 (max -1.0 (abs w)))))))

(defn between
  "The rotation taking `a` to `b`, both unit quaternions."
  [a b]
  (mul b (inverse a)))

(defn slerp
  "Along the shortest arc between two orientations."
  [a b t]
  (let [d (dot a b)
        ;; A quaternion and its negation are the same rotation; taking the
        ;; nearer one is what keeps the path short rather than going the
        ;; long way round.
        [b d] (if (neg? d) [(mapv - b) (- d)] [b d])
        t (double t)]
    (if (> d 0.9995)
      ;; Nearly parallel: lerp, because slerp divides by a sine going to
      ;; zero.
      (normalize (mapv (fn [x y] (+ x (* t (- y x)))) a b))
      (let [theta (math/acos (min 1.0 (max -1.0 d)))
            s (math/sin theta)
            wa (/ (math/sin (* (- 1.0 t) theta)) s)
            wb (/ (math/sin (* t theta)) s)]
        (normalize (mapv (fn [x y] (+ (* wa x) (* wb y))) a b))))))

(defn axis
  "The body's local `axis` expressed in world terms -- the rotation applied
  to a basis vector, which is what a joint frame is made of."
  [q a]
  (rotate q a))

(defn from-vectors
  "The shortest rotation taking one direction onto another.

  Used wherever something has to be pointed at something else: a cylinder
  along a link, a gripper at what it is reaching for. Two directions leave
  the roll about the shared axis undetermined, and this takes the choice
  that turns least.

  Exactly opposite directions have no shortest rotation -- every half turn
  about every perpendicular axis does the job -- so one perpendicular is
  picked, which is as good an answer as any and better than a NaN."
  [from to]
  (let [f (v/normalize from)
        t (v/normalize to)
        d (v/dot f t)]
    (cond
      (> d 0.999999) identity-q
      (< d -0.999999)
      (let [axis (v/cross f (if (< (abs (double (nth f 0))) 0.9) [1.0 0.0 0.0] [0.0 1.0 0.0]))]
        (from-axis-angle axis math/PI))
      :else (from-axis-angle (v/cross f t) (math/acos (min 1.0 (max -1.0 d)))))))

(defn to-matrix
  "The rotation as three rows, `[[r00 r01 r02] [r10 ...] [r20 ...]]`.

  The columns are where the x, y and z axes end up, which is what makes a
  matrix the convenient form when a calculation wants to read one axis out
  -- extracting joint angles from an orientation, say."
  [q]
  (let [[x y z] (rotate q [1.0 0.0 0.0])
        [a b c] (rotate q [0.0 1.0 0.0])
        [d e f] (rotate q [0.0 0.0 1.0])]
    [[x a d]
     [y b e]
     [z c f]]))

(defn from-matrix
  "A quaternion from three rows of a rotation matrix.

  Shepperd's method: four ways to compute the same quaternion, each stable
  where a different component is largest, and taking whichever the trace
  says is safest. The naive single formula divides by a number that goes
  to zero on a half turn."
  [[[r00 r01 r02] [r10 r11 r12] [r20 r21 r22]]]
  (let [trace (+ r00 r11 r22)]
    (normalize
     (cond
       (pos? trace)
       (let [s (* 2.0 (math/sqrt (+ 1.0 trace)))]
         [(/ (- r21 r12) s) (/ (- r02 r20) s) (/ (- r10 r01) s) (* 0.25 s)])

       (and (> r00 r11) (> r00 r22))
       (let [s (* 2.0 (math/sqrt (+ 1.0 r00 (- r11) (- r22))))]
         [(* 0.25 s) (/ (+ r01 r10) s) (/ (+ r02 r20) s) (/ (- r21 r12) s)])

       (> r11 r22)
       (let [s (* 2.0 (math/sqrt (+ 1.0 r11 (- r00) (- r22))))]
         [(/ (+ r01 r10) s) (* 0.25 s) (/ (+ r12 r21) s) (/ (- r02 r20) s)])

       :else
       (let [s (* 2.0 (math/sqrt (+ 1.0 r22 (- r00) (- r11))))]
         [(/ (+ r02 r20) s) (/ (+ r12 r21) s) (* 0.25 s) (/ (- r10 r01) s)])))))
