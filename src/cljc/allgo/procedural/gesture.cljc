(ns allgo.procedural.gesture
  "Procedural gesture: a figure animated by rhythm and noise, with no
  keyframes.

  Ebert et al., *Texturing & Modeling: A Procedural Approach*, Perlin's
  chapter 12, after his and Goldberg's Improv. A pose is a set of joint
  angles, and an action is not a recorded sequence of poses but a rule
  for each angle: a range it moves in, and a signal that moves it --
  a sine for anything with a beat, noise for anything that should never
  quite repeat. A dancer's elbows pump on the sine; the same figure at
  rest shifts its weight and glances around on noise, slowly, and that
  small, never-repeating motion is most of what makes it look alive
  rather than paused.

  Actions are layered and blended rather than switched. A body action
  sets every joint; an arm or a head action above it sets only the joints
  it names and lets the rest show through, so a figure can wave while it
  dances and nod while it waves. Changing an action fades its layer from
  the old to the new over `:transition` seconds, angle by angle, so
  nothing ever jumps: the angles at any instant are a weighted mix of
  smooth signals with smoothly changing weights.

  The legs are the exception: they are not animated but solved. Each
  foot keeps its place on the floor and the knee is put wherever two
  bones of fixed length allow, so when a dance bobs the hips the knees
  bend to absorb it. Animating them would let the feet slide or sink.

  A pose is `{joint [x y z]}` of Euler angles in radians (`quaternion/
  from-euler`'s order) plus `:root`, the hips' displacement in meters;
  `skeleton` says how the joints hang together, and `figure` turns a
  pose into points in space."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.procedural.noise :as noise]
            [allgo.procedural.shaping :as s]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; The figure

(def skeleton
  "Each joint's parent and its offset from the parent in the parent's
  frame, in meters, for a figure facing +z with its left on +x. The legs
  are listed for their lengths only; `figure` solves them."
  {:pelvis   {:parent nil :offset [0.0 0.98 0.0]}
   :chest    {:parent :pelvis :offset [0.0 0.3 0.0]}
   :neck     {:parent :chest :offset [0.0 0.2 0.0]}
   :head     {:parent :neck :offset [0.0 0.08 0.0]}
   :crown    {:parent :head :offset [0.0 0.22 0.02]}
   :l-shoulder {:parent :chest :offset [0.19 0.15 0.0]}
   :l-elbow  {:parent :l-shoulder :offset [0.0 -0.29 0.0]}
   :l-wrist  {:parent :l-elbow :offset [0.0 -0.26 0.0]}
   :l-hand   {:parent :l-wrist :offset [0.0 -0.08 0.0]}
   :r-shoulder {:parent :chest :offset [-0.19 0.15 0.0]}
   :r-elbow  {:parent :r-shoulder :offset [0.0 -0.29 0.0]}
   :r-wrist  {:parent :r-elbow :offset [0.0 -0.26 0.0]}
   :r-hand   {:parent :r-wrist :offset [0.0 -0.08 0.0]}})

(def ^:private chain-order
  "Parents before children."
  [:pelvis :chest :neck :head :crown
   :l-shoulder :l-elbow :l-wrist :l-hand
   :r-shoulder :r-elbow :r-wrist :r-hand])

(def ^:private hip-width 0.1)
(def ^:private thigh 0.44)
(def ^:private shin 0.43)
(def ^:private ankle-height 0.08)

(def bones
  "What to draw: pairs of points from `figure`."
  [[:pelvis :chest] [:chest :neck] [:neck :head]
   [:neck :l-shoulder] [:l-shoulder :l-elbow] [:l-elbow :l-wrist] [:l-wrist :l-hand]
   [:neck :r-shoulder] [:r-shoulder :r-elbow] [:r-elbow :r-wrist] [:r-wrist :r-hand]
   [:pelvis :l-hip] [:l-hip :l-knee] [:l-knee :l-ankle] [:l-ankle :l-toe]
   [:pelvis :r-hip] [:r-hip :r-knee] [:r-knee :r-ankle] [:r-ankle :r-toe]])

(defn- leg
  "Hip to a planted ankle by two bones: the knee where the two lengths
  meet, bent toward `forward`. A target out of reach is approached along
  the line to it, the leg straight."
  [hip ankle forward]
  (let [d (v/sub ankle hip)
        len (v/length d)
        u (v/scale d (/ 1.0 len))
        reach (* 0.999 (+ thigh shin))
        len (min len reach)
        a (/ (+ (* thigh thigh) (- (* shin shin)) (* len len)) (* 2.0 len))
        h (math/sqrt (max 0.0 (- (* thigh thigh) (* a a))))
        ;; The forward direction, made square to the leg.
        f (v/sub forward (v/scale u (v/dot forward u)))
        fl (v/length f)
        f (if (pos? fl) (v/scale f (/ 1.0 fl)) [0.0 0.0 1.0])
        knee (v/add (v/add hip (v/scale u a)) (v/scale f h))]
    [knee (v/add hip (v/scale u len))]))

(defn figure
  "The points of the figure in `pose`: `{joint [x y z]}` for every joint
  of `skeleton`, plus the solved legs -- hips, knees, ankles, toes --
  with the feet planted at `feet`, `[left right]` ankle positions on the
  floor (by default under the hips at rest)."
  ([pose] (figure pose nil))
  ([pose feet]
   (let [[rx ry rz] (:root pose [0.0 0.0 0.0])
         frames (reduce
                 (fn [acc j]
                   (let [{:keys [parent offset]} (skeleton j)
                         [ax ay az] (get pose j [0.0 0.0 0.0])
                         local (q/from-euler ax ay az)]
                     (if parent
                       (let [[prot ppos] (acc parent)]
                         (assoc acc j [(q/mul prot local) (v/add ppos (q/rotate prot offset))]))
                       (assoc acc j [local (v/add offset [rx ry rz])]))))
                 {} chain-order)
         points (update-vals frames second)
         [prot ppos] (frames :pelvis)
         ;; The legs go where the hips face, turned only about the vertical.
         [fx _ fz] (q/rotate prot [0.0 0.0 1.0])
         fl (math/hypot fx fz)
         forward (if (pos? fl) [(/ fx fl) 0.0 (/ fz fl)] [0.0 0.0 1.0])
         [lfoot rfoot] (or feet [[hip-width ankle-height 0.0] [(- hip-width) ankle-height 0.0]])
         side (fn [points hip-key knee-key ankle-key toe-key x foot]
                (let [hip (v/add ppos (q/rotate prot [x -0.05 0.0]))
                      [knee ankle] (leg hip foot forward)]
                  (assoc points hip-key hip knee-key knee ankle-key ankle
                         toe-key (v/add ankle (v/add (v/scale forward 0.15) [0.0 (- 0.02 ankle-height) 0.0])))))]
     (-> points
         (side :l-hip :l-knee :l-ankle :l-toe hip-width lfoot)
         (side :r-hip :r-knee :r-ankle :r-toe (- hip-width) rfoot)))))

;; ---------------------------------------------------------------------------
;; Signals and actions

(defn signal
  "A smooth signal in [-1, 1] at time `t`: `[:sin phase]`, a sine with
  `phase` in cycles, or `[:noise channel]`, gradient noise along time,
  one independent stream for each `channel` and `seed`. Both run at
  `rate` cycles a second; noise at that rate wanders rather than cycles."
  [[kind arg] ^double t ^double rate basis]
  (case kind
    :sin (math/sin (* 2.0 math/PI (+ (* rate t) (double arg))))
    :noise (max -1.0 (min 1.0 (* 1.8 (double (noise/sample basis (* rate t) (+ 0.5 (* 3.7 (double arg))) 0.31)))))))

(defn- degrees->radians [d] (* (double d) (/ math/PI 180.0)))

(defn- dof-value
  "One degree of freedom of an action at `t`: a number is held, and
  `[lo hi signal]`, optionally `[lo hi signal speed]`, swings between
  `lo` and `hi` -- degrees, or meters for the root -- as the signal goes
  from -1 to 1."
  [dof t rate basis]
  (if (number? dof)
    (double dof)
    (let [[lo hi sig speed] dof
          x (signal sig t (* (double rate) (double (or speed 1.0))) basis)]
      (s/mix (double lo) (double hi) (* 0.5 (+ 1.0 x))))))

(defn evaluate
  "An action's own angles at `t`: `{joint [x y z]}`, radians, for the
  joints it names, and `:root` in meters if it moves the hips."
  [{:keys [rate joints]} t basis]
  (into {}
        (for [[j dofs] joints]
          [j (mapv (fn [dof] (let [x (dof-value dof t rate basis)]
                               (if (= j :root) x (degrees->radians x))))
                   dofs)])))

(def actions
  "Every action, by layer: a body action sets every joint, and the arm
  and head actions over it set only theirs. `:rate` is cycles a second.
  The figure's left is +x, so a positive z turn lifts the left arm out
  to the side and a negative one the right."
  {:body
   {"rest"
    {:rate 0.25
     :joints {:root [0.0 0.0 0.0]
              :pelvis [0.0 0.0 0.0]
              ;; Breathing.
              :chest [[-1.0 2.0 [:sin 0.0]] 0.0 0.0]
              :neck [0.0 0.0 0.0]
              :head [[-2.0 2.0 [:noise 1] 0.6] [-6.0 6.0 [:noise 2] 0.5] 0.0]
              :l-shoulder [[-3.0 3.0 [:sin 0.1]] 0.0 [6.0 8.0 [:sin 0.0]]]
              :l-elbow [-8.0 0.0 0.0]
              :r-shoulder [[-3.0 3.0 [:sin 0.6]] 0.0 [-8.0 -6.0 [:sin 0.0]]]
              :r-elbow [-8.0 0.0 0.0]
              :l-wrist [0.0 0.0 0.0] :r-wrist [0.0 0.0 0.0]}}
    "fidget"
    {:rate 0.3
     :joints {;; The weight goes from one foot to the other, never in step.
              :root [[-0.06 0.06 [:noise 3]] [-0.02 0.0 [:noise 4]] [-0.03 0.03 [:noise 5]]]
              :pelvis [0.0 [-12.0 12.0 [:noise 6]] [-5.0 5.0 [:noise 3]]]
              :chest [[-2.0 6.0 [:noise 7]] [-10.0 10.0 [:noise 8]] [-4.0 4.0 [:noise 9]]]
              :neck [0.0 [-10.0 10.0 [:noise 10]] 0.0]
              :head [[-12.0 10.0 [:noise 11] 1.5] [-30.0 30.0 [:noise 12] 1.2] [-6.0 6.0 [:noise 13]]]
              :l-shoulder [[-15.0 10.0 [:noise 14]] 0.0 [4.0 14.0 [:noise 15]]]
              :l-elbow [[-35.0 -5.0 [:noise 16]] 0.0 0.0]
              :r-shoulder [[-15.0 10.0 [:noise 17]] 0.0 [-14.0 -4.0 [:noise 18]]]
              :r-elbow [[-35.0 -5.0 [:noise 19]] 0.0 0.0]
              :l-wrist [0.0 0.0 0.0] :r-wrist [0.0 0.0 0.0]}}
    "dance"
    {:rate 1.0
     :joints {;; Down on every beat, side to side every two.
              :root [[-0.05 0.05 [:sin 0.0] 0.5] [-0.09 -0.01 [:sin 0.25]] 0.0]
              :pelvis [0.0 [-15.0 15.0 [:sin 0.0] 0.5] [-6.0 6.0 [:sin 0.5] 0.5]]
              :chest [[0.0 8.0 [:sin 0.25]] [-20.0 20.0 [:sin 0.5] 0.5] [-8.0 8.0 [:sin 0.0] 0.5]]
              :neck [0.0 0.0 [-8.0 8.0 [:sin 0.0] 0.5]]
              :head [[-10.0 6.0 [:sin 0.75]] [-10.0 10.0 [:noise 20]] 0.0]
              :l-shoulder [[-60.0 20.0 [:sin 0.0]] 0.0 [10.0 35.0 [:sin 0.5] 0.5]]
              :l-elbow [[-110.0 -60.0 [:sin 0.25]] 0.0 0.0]
              :r-shoulder [[-60.0 20.0 [:sin 0.5]] 0.0 [-35.0 -10.0 [:sin 0.0] 0.5]]
              :r-elbow [[-110.0 -60.0 [:sin 0.75]] 0.0 0.0]
              :l-wrist [[-20.0 20.0 [:sin 0.0]] 0.0 0.0]
              :r-wrist [[-20.0 20.0 [:sin 0.5]] 0.0 0.0]}}
    "tired"
    {:rate 0.18
     :joints {:root [0.0 -0.05 [-0.02 0.02 [:noise 21]]]
              :pelvis [[4.0 8.0 [:noise 22]] 0.0 0.0]
              ;; Slumped, and breathing slow and deep.
              :chest [[14.0 22.0 [:sin 0.0]] [-4.0 4.0 [:noise 23]] 0.0]
              :neck [[15.0 22.0 [:noise 24]] 0.0 0.0]
              :head [[10.0 25.0 [:noise 25] 0.7] [-8.0 8.0 [:noise 26]] [-5.0 5.0 [:noise 27]]]
              :l-shoulder [[2.0 6.0 [:sin 0.2]] 0.0 3.0]
              :l-elbow [-4.0 0.0 0.0]
              :r-shoulder [[2.0 6.0 [:sin 0.2]] 0.0 -3.0]
              :r-elbow [-4.0 0.0 0.0]
              :l-wrist [0.0 0.0 0.0] :r-wrist [0.0 0.0 0.0]}}}
   :arms
   {"none" {:rate 1.0 :joints {}}
    "wave"
    {:rate 1.4
     :joints {:r-shoulder [[-15.0 -5.0 [:noise 30] 0.3] 0.0 [-150.0 -140.0 [:noise 31] 0.3]]
              :r-elbow [[-30.0 -20.0 [:noise 32] 0.3] 0.0 [-40.0 35.0 [:sin 0.0]]]
              :r-wrist [0.0 0.0 [-15.0 15.0 [:sin 0.1]]]}}
    "conduct"
    {:rate 0.9
     :joints {;; One beat in the right hand, and the left shaping the phrase.
              :r-shoulder [[-75.0 -45.0 [:sin 0.0]] 0.0 [-35.0 -20.0 [:sin 0.25]]]
              :r-elbow [[-80.0 -40.0 [:sin 0.1]] 0.0 [-15.0 15.0 [:sin 0.35]]]
              :r-wrist [[-25.0 25.0 [:sin 0.2]] 0.0 0.0]
              :l-shoulder [[-80.0 -30.0 [:noise 33] 0.4] 0.0 [15.0 40.0 [:noise 34] 0.4]]
              :l-elbow [[-90.0 -40.0 [:noise 35] 0.4] 0.0 0.0]
              :l-wrist [[-20.0 20.0 [:noise 36] 0.6] 0.0 0.0]}}
    "shrug"
    {:rate 0.5
     :joints {:l-shoulder [[-30.0 -25.0 [:noise 37]] 0.0 [28.0 34.0 [:noise 38]]]
              ;; Forearms forward, turned out: palms up.
              :l-elbow [[-95.0 -85.0 [:noise 39]] 0.0 [30.0 40.0 [:noise 40]]]
              :r-shoulder [[-30.0 -25.0 [:noise 41]] 0.0 [-34.0 -28.0 [:noise 42]]]
              :r-elbow [[-95.0 -85.0 [:noise 43]] 0.0 [-40.0 -30.0 [:noise 44]]]
              :neck [0.0 0.0 [-8.0 -4.0 [:noise 45]]]}}}
   :head
   {"none" {:rate 1.0 :joints {}}
    "nod"
    {:rate 1.6
     :joints {:head [[-4.0 22.0 [:sin 0.0]] [-3.0 3.0 [:noise 50] 0.3] 0.0]}}
    "shake"
    {:rate 1.8
     :joints {:head [[2.0 6.0 [:noise 51] 0.3] [-28.0 28.0 [:sin 0.0]] 0.0]}}
    "look around"
    {:rate 0.35
     :joints {:neck [[-6.0 6.0 [:noise 52]] [-25.0 25.0 [:noise 53]] 0.0]
              :head [[-18.0 12.0 [:noise 54]] [-45.0 45.0 [:noise 55]] [-8.0 8.0 [:noise 56]]]}}}})

(def layer-order "Bottom to top." [:body :arms :head])

;; ---------------------------------------------------------------------------
;; Blending

(defn start
  "A performance starting at time `t`: each layer in its first action,
  fully in."
  ([t] (start t {}))
  ([t {:keys [body arms head] :or {body "rest" arms "none" head "none"}}]
   {:layers {:body {:from body :to body :at (double t)}
             :arms {:from arms :to arms :at (double t)}
             :head {:from head :to head :at (double t)}}}))

(defn weight
  "How far a layer has come from its old action to its new one at `t`,
  0 to 1, over `transition` seconds, eased at both ends."
  ^double [{:keys [at]} ^double t ^double transition]
  (if (<= transition 0.0) 1.0 (s/smoothstep 0.0 1.0 (/ (- t (double at)) transition))))

(defn change
  "The performance with `layer` turning to `action` from time `t`. A
  change in the middle of another starts from where that one had got
  to: the layer's old action is replaced by a frozen blend of the two,
  so the angles carry on from their current values instead of jumping."
  [{:keys [layers] :as performance} layer action t transition]
  (let [{:keys [from to] :as l} (layers layer)
        w (weight l (double t) (double transition))]
    (assoc-in performance [:layers layer]
              {:from (if (or (>= w 1.0) (= from to)) to {:blend [from to] :weight w})
               :to action
               :at (double t)})))

(defn- mix-pose
  "`a` toward `b` by `w`, joint by joint, over the joints of both; a
  joint missing from one side is taken from `under`."
  [under a b w]
  (into under
        (for [j (distinct (concat (keys a) (keys b)))
              :let [base (get under j [0.0 0.0 0.0])
                    x (get a j base) y (get b j base)]]
          [j (mapv #(s/mix (double %1) (double %2) w) x y)])))

(defn- action-pose
  "What an action -- a name in `layer`, or a frozen blend of two -- puts
  over `under` at `t`."
  [layer action under t basis]
  (if (map? action)
    (let [{[a b] :blend w :weight} action]
      (mix-pose under (action-pose layer a under t basis) (action-pose layer b under t basis) w))
    (merge under (evaluate (get-in actions [layer action]) t basis))))

(defn pose
  "The pose of a performance at time `t`: each layer's blend of its old
  and new action laid over the layers below it, `:transition` seconds
  apart, with `:seed` for the noise. To run faster or slower, run `t`
  so: scaling it here would jump every phase when the speed changed."
  [{:keys [layers]} t {:keys [transition seed] :or {transition 1.0 seed 0}}]
  (let [basis (noise/gradient-basis {:seed seed})]
    (reduce (fn [under layer]
              (let [{:keys [from to] :as l} (layers layer)
                    w (weight l (double t) (double transition))
                    a (action-pose layer from under t basis)
                    b (action-pose layer to under t basis)]
                (mix-pose under a b w)))
            {} layer-order)))
