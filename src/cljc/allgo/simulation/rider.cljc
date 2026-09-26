(ns allgo.simulation.rider
  "Someone on the motorcycle: a jointed figure, pinned to the bike at the
  seat, the grips and the pegs, and holding a riding posture with its own
  muscles.

  ## Two models and five pins

  The rider and the bike cannot be one articulated model. Two hands on
  one handlebar, a seat and two pegs close five loops through the bike's
  frame, and reduced coordinates describe trees: there is no way to say
  that the end of the left forearm and the end of the right both sit on
  the steering head. So the rider is a second model of its own -- a
  ragdoll with a free pelvis -- and the loops are closed by pins
  (`allgo.physics.world`), which are solved as bilateral impulses in the
  same sweep as the tires' contacts. Inside each model the joints are
  exact; between them they are iterated. That is the right split, because
  it puts the iteration where the loop is and nowhere else.

  ## Built sitting down

  The figure is laid out in its riding posture from the start: its rest
  configuration, every joint at zero, *is* the rider on the bike. The
  pelvis sits on the seat; each arm is a two-bone reach from the shoulder
  to its grip and each leg from the hip to its peg, with the elbow and
  knee put where two bones of those lengths have to meet. So nothing has
  to be solved to seat the rider, and a muscle holding a joint at zero is
  holding it in the riding pose.

  ## Muscles

  Each joint is the right kind -- a hinge at the elbows and knees, a ball
  joint at the hips, shoulders, spine and neck -- and each is sprung
  toward its rest angle with a spring and a damper, as a generalized
  force on that joint. Stiff enough that the rider sits up and does not
  fold over the tank; soft enough that the arms follow the bars and the
  body sways as the bike leans and pitches. The torques are a proportional
  and derivative pull on the joint's own coordinate: an angle for a
  hinge, the rotation vector of its quaternion for a ball joint.

  They also lean. Asked for a turn, the rider tips their upper body into
  it (`weight-shift`), which is the other half of how a bike is steered:
  the hands countersteer it over, the body goes with it.

  In a crash the pins let go and the muscles switch off altogether:
  what falls off is a plain ragdoll, held together only by its joint
  limits, and it lands on the bike as well as the road."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]))

;; ---------------------------------------------------------------------------
;; Where the rider touches the bike, in the bike frame's coordinates

(defn anchors
  "The seat, the grips and the pegs, on a bike whose steering head is at
  `head`. The grips are given in the steering head's own frame too,
  because that is what they are fixed to."
  [{:keys [head]}]
  (let [grip-l [0.0 0.12 -0.30]
        grip-r [0.0 0.12 0.30]]
    {:seat [-0.24 0.18 0.0]
     :grip-local {:left grip-l :right grip-r}
     :grips {:left (v/add head grip-l) :right (v/add head grip-r)}
     :pegs {:left [-0.10 -0.14 -0.20] :right [-0.10 -0.14 0.20]}}))

;; ---------------------------------------------------------------------------
;; Geometry

(defn- rotation-from-y
  "The rotation taking +y onto `d`."
  [d]
  (let [d (v/normalize d)
        axis (v/cross [0.0 1.0 0.0] d)
        len (v/length axis)
        c (double (nth d 1))]
    (cond
      (> len 1e-9) (q/from-axis-angle (v/scale axis (/ 1.0 len)) (Math/atan2 len c))
      (pos? c) q/identity-q
      :else (q/from-axis-angle [1.0 0.0 0.0] Math/PI))))

(defn- box-inertia [[sx sy sz] m]
  (let [f (/ (double m) 12.0)
        sx (double sx) sy (double sy) sz (double sz)]
    [[(* f (+ (* sy sy) (* sz sz))) 0.0 0.0]
     [0.0 (* f (+ (* sx sx) (* sz sz))) 0.0]
     [0.0 0.0 (* f (+ (* sx sx) (* sy sy)))]]))

(defn- rotate-inertia
  "`R I R^T`: an inertia given about a shape's own axes, turned into the
  frame the shape is turned by `rot` in."
  [inertia rot]
  (let [cols (mapv #(q/rotate rot %) [[1.0 0.0 0.0] [0.0 1.0 0.0] [0.0 0.0 1.0]])
        r (fn [i k] (double (nth (nth cols k) i)))]
    (vec (for [i (range 3)]
           (vec (for [j (range 3)]
                  (reduce + (for [k (range 3) l (range 3)]
                              (* (r i k) (double (nth (nth inertia k) l)) (r j l))))))))))

(defn- bend
  "Where the joint between two bones of `len1` and `len2` goes, reaching
  from `a` to `c` and bending toward `hint`: an elbow or a knee."
  [a c len1 len2 hint]
  (let [len1 (double len1) len2 (double len2)
        ac (v/sub c a)
        d (min (- (+ len1 len2) 1e-3) (v/length ac))
        u (v/normalize ac)
        along (/ (+ (- (* len1 len1) (* len2 len2)) (* d d)) (* 2.0 d))
        h (Math/sqrt (max 0.0 (- (* len1 len1) (* along along))))
        perp (v/normalize (v/sub hint (v/scale u (v/dot hint u))))]
    (v/add a (v/add (v/scale u along) (v/scale perp h)))))

(defn- angle-between ^double [a b]
  (Math/acos (max -1.0 (min 1.0 (v/dot (v/normalize a) (v/normalize b))))))

(defn- bone
  "A limb from its joint to `bone` -- the vector to its far end -- with
  the mass in the middle and a box of the right thickness along it. Its
  frame is at the joint and unrotated, like every frame in this figure,
  which is what lets every position here be written in one set of
  coordinates."
  [parent joint origin bone mass width depth]
  (let [len (v/length bone)
        rot (rotation-from-y bone)
        size [depth len width]]
    {:parent parent :joint joint :origin {:rot nil :pos origin}
     :mass mass :com (v/scale bone 0.5)
     :inertia (rotate-inertia (box-inertia size mass) rot)
     :shape :box :size size
     :shape-pose {:rot rot :pos (v/scale bone 0.5)}}))

;; ---------------------------------------------------------------------------
;; The figure

(def ^:private upper-arm 0.30)
(def ^:private forearm 0.30)
(def ^:private thigh 0.42)
(def ^:private shin 0.42)

(defn layout
  "Every joint of the seated rider, in the bike frame's coordinates."
  [bike-cfg]
  (let [{:keys [seat grips pegs]} (anchors bike-cfg)
        pelvis (v/add seat [0.02 0.09 0.0])
        spine (v/add pelvis [0.0 0.08 0.0])
        chest (v/add pelvis [0.22 0.43 0.0])
        neck (v/add chest [0.02 0.05 0.0])
        side (fn [s] (if (= s :left) -1.0 1.0))
        shoulder (fn [s] (v/add chest [0.0 0.0 (* 0.19 (side s))]))
        hip (fn [s] (v/add pelvis [0.0 -0.04 (* 0.10 (side s))]))
        elbow (fn [s] (bend (shoulder s) (grips s) upper-arm forearm [0.0 -1.0 (* 0.6 (side s))]))
        knee (fn [s] (bend (hip s) (pegs s) thigh shin [1.0 0.4 (* 0.3 (side s))]))]
    {:seat seat :pelvis pelvis :spine spine :chest chest :neck neck
     :head (v/add neck [0.05 0.22 0.0])
     :shoulder {:left (shoulder :left) :right (shoulder :right)}
     :elbow {:left (elbow :left) :right (elbow :right)}
     :hand grips
     :hip {:left (hip :left) :right (hip :right)}
     :knee {:left (knee :left) :right (knee :right)}
     :foot pegs}))

(def pelvis-link -1)
(def torso 0)
(def head 1)
(def upper-arms {:left 2 :right 4})
(def forearms {:left 3 :right 5})
(def thighs {:left 6 :right 8})
(def shins {:left 7 :right 9})

(defn- hinge
  "A hinge at the join of `upper` and `lower`, about the axis their bend
  turns about, allowed to straighten fully and to fold to `fold` radians.
  Bending is positive: turning about `upper x lower` opens the angle
  between the two bones."
  [link upper lower fold]
  (let [axis (v/normalize (v/cross upper lower))
        bent (angle-between upper lower)]
    (assoc link :axis axis :limit [(- bent) (- (double fold) bent)])))

(defn model
  "The rider as an articulated model, its pelvis the free root, laid out
  on the bike described by `bike-cfg`."
  [bike-cfg]
  (let [{:keys [pelvis spine chest neck head shoulder elbow hand hip knee foot]}
        (layout bike-cfg)
        arm (fn [s]
              (let [up (v/sub (elbow s) (shoulder s))
                    low (v/sub (hand s) (elbow s))]
                [(assoc (bone torso :spherical (v/sub (shoulder s) spine) up 2.2 0.09 0.09)
                        :cone 1.4 :twist 1.0)
                 (hinge (bone (upper-arms s) :revolute up low 1.8 0.08 0.08) up low 2.5)]))
        leg (fn [s]
              (let [up (v/sub (knee s) (hip s))
                    low (v/sub (foot s) (knee s))]
                [(assoc (bone -1 :spherical (v/sub (hip s) pelvis) up 9.0 0.15 0.15)
                        :cone 1.1 :twist 0.6)
                 (hinge (bone (thighs s) :revolute up low 4.5 0.11 0.11) up low 2.6)]))
        [ual fal] (arm :left)
        [uar far] (arm :right)
        [thl shl] (leg :left)
        [thr shr] (leg :right)]
    {:base {:mass 11.0 :com [0.0 0.0 0.0]
            :inertia (box-inertia [0.20 0.16 0.34] 11.0)
            :shape :box :size [0.20 0.16 0.32]
            :shape-pose {:rot q/identity-q :pos [0.0 0.0 0.0]}}
     :links
     [(assoc (bone -1 :spherical (v/sub spine pelvis) (v/sub chest spine) 30.0 0.36 0.20)
             :cone 0.6 :twist 0.5)
      ;; A head is round, and a helmet more so.
      (assoc (bone torso :spherical (v/sub neck spine) (v/sub head neck) 5.0 0.18 0.20)
             :cone 0.8 :twist 1.0
             :shape :ball :radius 0.13
             :shape-pose {:rot q/identity-q :pos (v/scale (v/sub head neck) 0.6)})
      ual fal uar far
      thl shl thr shr]}))

(defn mass
  "What the rider weighs, all of them."
  [rider]
  (+ (double (get-in rider [:base :mass])) (reduce + (map :mass (:links rider)))))

;; ---------------------------------------------------------------------------
;; Pins

(defn pins
  "The five places the rider is fixed to the bike, as world pins between
  model `bike` (its frame, and `steering` for the grips) and model
  `rider`. Each lets go past its `:break-force`: a rider is thrown when
  something asks the seat or the grips to hold more than a person would.
  A landing off a well-built jump stays under them; coming down flat
  from a couple of meters does not."
  [bike-cfg bike rider steering]
  (let [{:keys [seat pelvis elbow hand knee foot]} (layout bike-cfg)
        {:keys [grip-local pegs]} (anchors bike-cfg)]
    (into [{:a [:link rider pelvis-link] :pa (v/sub seat pelvis)
            :b [:link bike -1] :pb seat :break-force 20000.0 :part :seat}]
          (for [s [:left :right]
                pin [{:a [:link rider (forearms s)] :pa (v/sub (hand s) (elbow s))
                      :b [:link bike steering] :pb (grip-local s)
                      :break-force 4500.0 :part [:hand s]}
                     {:a [:link rider (shins s)] :pa (v/sub (foot s) (knee s))
                      :b [:link bike -1] :pb (pegs s)
                      :break-force 8000.0 :part [:foot s]}]]
            pin))))

;; ---------------------------------------------------------------------------
;; Muscles

(def ^:private muscle
  "How hard each joint holds the riding pose: N m per radian and N m s per
  radian. The spine and hips carry the body and are the stiffest; the
  arms are soft enough to follow the bars."
  {torso [900.0 40.0]
   head [60.0 3.0]
   (upper-arms :left) [120.0 4.0] (upper-arms :right) [120.0 4.0]
   (forearms :left) [60.0 2.0] (forearms :right) [60.0 2.0]
   (thighs :left) [500.0 20.0] (thighs :right) [500.0 20.0]
   (shins :left) [300.0 10.0] (shins :right) [300.0 10.0]})

(defn- rotation-vector
  "A ball joint's rotation as an axis scaled by its angle -- the three
  numbers a spring toward rest pulls against."
  [[x y z w]]
  (let [[x y z w] (if (neg? (double w)) [(- x) (- y) (- z) (- w)] [x y z w])
        s (Math/sqrt (+ (* x x) (* y y) (* z z)))]
    (if (< s 1e-9)
      [(* 2.0 x) (* 2.0 y) (* 2.0 z)]
      (let [k (/ (* 2.0 (Math/atan2 s w)) s)]
        [(* k x) (* k y) (* k z)]))))

(defn muscles
  "The torque on every joint of the rider, one per degree of freedom, as
  a spring and damper toward the riding pose. `tone` scales both: 1 is a
  rider riding, zero a ragdoll.

  `targets` moves a ball joint's rest somewhere else -- link index to the
  rotation it should hold instead of none -- which is how the rider
  shifts their weight: the spine is asked to lean, and the same spring
  that held it upright now holds it over."
  ([rider pose tone] (muscles rider pose tone nil))
  ([rider {:keys [q qd]} tone targets]
   (let [tone (double tone)]
     (loop [i 0 offset 0 tau (transient [])]
       (if (= i (count (:links rider)))
         (persistent! tau)
         (let [link (nth (:links rider) i)
               [kp kd] (muscle i)
               kp (* tone (double kp))
               kd (* tone (double kd))]
           (case (:joint link)
             :revolute
             (recur (inc i) (inc offset)
                    (conj! tau (- (- (* kp (double (nth q i))))
                                  (* kd (double (nth qd offset))))))
             :spherical
             (let [rv (rotation-vector (if-let [t (get targets i)]
                                        ;; What is left of the joint's
                                        ;; rotation once the target's is
                                        ;; taken out of it.
                                         (q/mul (q/conjugate t) (nth q i))
                                         (nth q i)))
                   t (fn [k] (- (- (* kp (double (nth rv k))))
                                (* kd (double (nth qd (+ offset k))))))]
               (recur (inc i) (+ offset 3)
                      (-> tau (conj! (t 0)) (conj! (t 1)) (conj! (t 2))))))))))))

(def ^:private max-shift
  "How far a rider leans their upper body into a turn, at most."
  0.35)

(defn weight-shift
  "Where the rider holds their spine and head to help a lean of `lean`
  radians along: the torso tipped toward the inside of the turn by part
  of it -- moving their weight to where the bike is going -- and the head
  tipped part of the way back, to keep the eyes nearer level.

  The figure's frames face +x with +y up, so a lean toward -z, the
  positive way, is a turn of the torso about +x the negative way."
  [lean]
  (let [shift (max (- max-shift) (min max-shift (* 0.6 (double lean))))]
    {torso (q/from-axis-angle [1.0 0.0 0.0] (- shift))
     head (q/from-axis-angle [1.0 0.0 0.0] (* 0.6 shift))}))

;; ---------------------------------------------------------------------------
;; Getting on

(defn start-pose
  "The rider in the saddle of a bike whose frame is at `bike-base`: every
  joint at rest, the pelvis where the layout puts it relative to the
  frame, moving with it."
  [bike-cfg rider bike-base]
  (let [{:keys [pelvis]} (layout bike-cfg)
        {:keys [rot pos vel]} bike-base
        [wx wy wz vx vy vz] (or vel [0.0 0.0 0.0 0.0 0.0 0.0])
        ;; The pelvis is carried by the frame: its velocity is the frame's
        ;; plus the frame's spin times the lever between them, both in the
        ;; frame's own coordinates -- which, the two frames being aligned,
        ;; are the pelvis's too.
        w [wx wy wz]
        lin (v/add [vx vy vz] (v/cross w pelvis))]
    {:q (mapv (fn [l] (if (= :spherical (:joint l)) q/identity-q 0.0)) (:links rider))
     :qd (vec (repeat (reduce + (map #(if (= :spherical (:joint %)) 3 1) (:links rider))) 0.0))
     :base {:rot rot
            :pos (v/add pos (q/rotate rot pelvis))
            :vel (into (vec w) lin)}}))
