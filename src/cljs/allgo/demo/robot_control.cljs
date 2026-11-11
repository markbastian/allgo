(ns allgo.demo.robot-control
  "Canvas demo for the rest of chapter 26 of Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition: control,
  uncertainty, people, and reactive robots.

  Computed torque: a two-link arm swinging between two poses three ways
  at once -- open loop from a model a fifth too light (red), PD feedback
  alone (blue) and computed torque from the same wrong model (green) --
  over the reference (faint). MPC: a pendulum swung up by an iLQR plan,
  run open loop on a pendulum heavier than the plan assumed (red), and
  by model predictive control replanning every step with the same wrong
  model (green). Peg in hole: guarded moves under 20 degrees of motion
  error -- straight down, and the two moves that always get in. Hexapod:
  a six-legged walker's tripod gait, a leg lifting higher each time a
  rock blocks its swing. Predicting people: a person crossing a room, the
  robot's belief over where they are going growing from their steps.
  Terrain: a laser's verdicts on the ground just ahead teaching a camera
  what drivable looks like, road first and then grass.

  From LaValle's *Planning Algorithms* (Cambridge University Press,
  2006): the time-optimal speed along a winding road (section 14.6.3),
  as fast as the top speed, the braking and the grip in the bends allow
  -- full throttle or full brakes everywhere it is not at a limit; and
  two tests of a linear system (section 15.1): Lyapunov's, an energy
  that falls along every trajectory, its level sets crossed only
  inward, and Kalman's controllability, the states one bounded push
  after another can reach -- a widening parallelogram when the system is
  controllable, a line when it is not."
  (:require [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as linear-systems]
            [allgo.random :as random]
            [allgo.robotics.control :as control]
            [allgo.robotics.guarded :as guarded]
            [allgo.robotics.humans :as humans]
            [allgo.robotics.reactive :as reactive]
            [allgo.robotics.terrain :as terrain]
            [allgo.robotics.tracking :as track]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:method "computed torque" :seed 1 :speed 2})

;; ---------------------------------------------------------------------------
;; Models and scenes

(def ^:private arm {:m1 1.0 :m2 1.0 :l1 1.0 :l2 1.0 :g 9.81})
(def ^:private light-arm (assoc arm :m1 0.8 :m2 0.8))
(def ^:private poses [[-0.5 0.8] [0.9 -0.4]])

(defn- pendulum [m]
  (fn [[th w] [u]] (let [a (+ (- (* 9.81 (math/sin th))) (* -0.1 w) (/ u m))]
                     [(+ th (* 0.05 w)) (+ w (* 0.05 a))])))

(def ^:private swing {:dynamics (pendulum 1.0) :x0 [0.0 0.0] :goal [math/PI 0.0]
                      :Q [[0.1 0.0] [0.0 0.01]] :R [[0.001]] :Qf [[100.0 0.0] [0.0 10.0]]})

(def ^:private hole {:width 0.4 :depth 0.6})

(def ^:private room-goals {"window" [10.0 5.0] "kitchen" [0.0 5.0] "hallway" [5.0 10.0]})
(def ^:private steps8 (vec (for [a (range 8)] [(* 0.5 (math/cos (* a (/ math/PI 4)))) (* 0.5 (math/sin (* a (/ math/PI 4))))])))

(def ^:private ground-colors {:road [0.5 0.5 0.52] :grass [0.25 0.6 0.2] :rock [0.55 0.4 0.25]})

(defn- ground [rng road-left road-right]
  (let [kind (fn [x y] (cond (and (<= road-left x) (< x road-right)) :road
                             (zero? (mod (+ (* 7 x) (* 13 y)) 23)) :rock
                             :else :grass))]
    {:image {:w 40 :h 30 :px (vec (for [y (range 30) x (range 40)] (mapv #(random/gaussian rng % 0.03) (ground-colors (kind x y)))))}
     :kind kind :flat? (fn [x y] (not= :rock (kind x y)))}))

;; the time-optimal speed's road, a point every few centimetres
(def ^:private road
  (vec (for [i (range 1201) :let [t (* 12.0 (/ i 1200.0))]]
         [t (+ (* 1.6 (math/sin (* 0.7 t))) (* 0.5 (math/sin (* 2.1 t))))])))
(def ^:private limits {:vmax 4.0 :amax 1.5 :alat 2.0})

(def ^:private stable [[0.9 0.2] [-0.1 0.8]])
(def ^:private unstable [[1.02 0.2] [-0.2 1.0]])

(defn- zonotope
  "The states reachable from rest in `k` steps of x' = A x + B u, |u| <= 1,
  B a single column: the Minkowski sum of the segments +-A^i B, as its
  polygon's corners."
  [A B k]
  (let [gens (map (fn [[[x] [y]]] (if (or (neg? y) (and (zero? y) (neg? x))) [(- x) (- y)] [x y]))
                  (take k (iterate (partial lin/mat-mul A) B)))
        sorted (sort-by (fn [[x y]] (math/atan2 y x)) gens)
        start (reduce (fn [[a b] [x y]] [(- a x) (- b y)]) [0.0 0.0] sorted)
        up (reductions (fn [[a b] [x y]] [(+ a (* 2 x)) (+ b (* 2 y))]) start sorted)]
    (vec (concat up (rest (reductions (fn [[a b] [x y]] [(- a (* 2 x)) (- b (* 2 y))]) (last up) sorted))))))

(defn- init [method]
  (let [rng (random/rng (.-seed controls))]
    (case method
      "computed torque"
      (let [traj (track/retime poses 1.0 2.0)
            back (track/retime (vec (reverse poses)) 1.0 2.0)
            ;; there and back, again and again
            looped {:at (fn [t] (let [d (:duration traj) t (mod t (* 2 (+ d 0.5)))]
                                  (if (< t (+ d 0.5)) ((:at traj) t) ((:at back) (- t d 0.5)))))}]
        {:arms (mapv #(vec (take 4000 (track/track-seq arm % looped 0.002)))
                     [(track/feedforward light-arm) (track/pd 100 20) (track/computed-torque light-arm 100 20)])
         :tick 0})
      "MPC"
      (let [plan (:us (nth (control/ilqr-seq swing (repeat 60 [0.0])) 30))
            heavy (pendulum 1.3)]
        {:open (control/rollout heavy [0.0 0.0] (concat plan (repeat 200 [0.0])))
         :mpc (mapv :x (take 260 (control/mpc-seq heavy swing 40 3)))
         :tick 0})
      "peg in hole"
      (let [cone (math/to-radians 20)
            down (guarded/direction (- (/ math/PI 2)))
            left (guarded/direction (- (+ (/ math/PI 2) (math/to-radians 30))))
            right (guarded/direction (math/to-radians -45))
            bottom (fn [p _] (guarded/at-bottom? hole p))]
        {:trials (vec (for [i (range 24)]
                        (if (even? i)
                          {:strategy :direct :path (guarded/guarded-moves hole [0.0 2.0] [[down guarded/on-contact]] cone 0.02 rng)}
                          {:strategy :guarded :path (guarded/guarded-moves hole [0.0 2.0] [[left guarded/on-contact] [right bottom]] cone 0.02 rng)})))
         :tick 0 :trial 0})
      "hexapod"
      (let [rocks (into {} (for [s (range 0 200 5)] [[(nth reactive/legs (mod s 6)) s] (+ 0.12 (* 0.2 (rng)))]))]
        {:states (vec (take 400 (reactive/walk-seq (reactive/walker 0.1) (fn [leg s] (get rocks [leg s] 0.0)) 0.1 0.05)))
         :rocks rocks :tick 0})
      "predicting people"
      (let [walk (vec (reductions (fn [[x y] [dx dy]] [(+ x dx) (+ y dy)]) [5.0 2.0]
                                  (concat (repeat 4 [0.0 0.5]) (repeat 8 [0.35355 0.35355]) (repeat 3 [0.5 0.0]))))
            steps (map (fn [[a b]] [a (mapv - b a)]) (partition 2 1 walk))]
        {:walk walk :beliefs (vec (humans/belief-seq room-goals steps steps8 2.0)) :tick 0})
      "time-optimal speed"
      (assoc (track/speed-profile road limits) :tick 0)
      "controllability and Lyapunov"
      {:tick 0}
      "terrain"
      (let [frames [(ground rng 12 28) (ground rng 0 2)]
            models (vec (map :model (terrain/adapt-seq frames [15 25 25 30] 2 20 2000 rng)))]
        {:frames frames :models models :tick 0}))))

;; ---------------------------------------------------------------------------

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        model (atom nil)]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! [] (reset! model (init (.-method controls))))
            (line! [[ax ay] [bx by] color width]
              (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
              (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx by) (.stroke ctx))
            (dot! [[x y] r color]
              (set! (.-fillStyle ctx) color) (.beginPath ctx) (.arc ctx x y r 0 (* 2 math/PI)) (.fill ctx))
            (text! [s x y color & [size]]
              (set! (.-fillStyle ctx) color) (set! (.-font ctx) (str (or size 12) "px ui-monospace,Menlo,monospace"))
              (.fillText ctx s x y))
            (draw! []
              (let [cw (.-width canvas) ch (.-height canvas) method (.-method controls)]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (swap! model update :tick + (.-speed controls))
                (set! (.-textContent caption)
                      (str method "  |  "
                           (case method
                             "computed torque" (draw-arms! cw ch)
                             "MPC" (draw-pendulums! cw ch)
                             "peg in hole" (draw-peg! cw ch)
                             "hexapod" (draw-hexapod! cw ch)
                             "predicting people" (draw-people! cw ch)
                             "terrain" (draw-terrain! cw ch)
                             "time-optimal speed" (draw-speed! cw ch)
                             "controllability and Lyapunov" (draw-systems! cw ch))))))
            (draw-speed! [cw ch]
              (let [{:keys [s v t tick]} @model
                    total (peek t)
                    now (mod (* 0.01 tick) (+ total 1.0))
                    k (min (dec (count t)) (count (take-while #(< % now) t)))
                    sx (/ (- cw 80) 12.0) mid (* 0.3 ch)
                    px (fn [[x y]] [(+ 40 (* sx x)) (- mid (* sx 0.45 y))])
                    {:keys [vmax alat]} limits
                    fastest (reduce max v)
                    cap (mapv (fn [[a b c]] (let [kk (track/curvature a b c)] (if (pos? kk) (min vmax (math/sqrt (/ alat kk))) vmax)))
                              (partition 3 1 road))
                    [gx gy gw gh] [40 (* 0.58 ch) (- cw 80) (* 0.3 ch)]
                    gpx (fn [si vi] [(+ gx (* gw (/ si (peek s)))) (- (+ gy gh) (* gh (/ vi (+ vmax 0.5))))])]
                ;; the road, colored by speed
                (doseq [i (range 0 (dec (count road)) 4)]
                  (line! (px (road i)) (px (road (min (dec (count road)) (+ i 4))))
                         (str "hsl(" (- 220 (* 220 (/ (v i) fastest))) ",80%,60%)") 5))
                (dot! (px (road k)) 9 "#ffffff")
                ;; speed against distance: the cap from the bends, and the profile under it
                (set! (.-fillStyle ctx) "#161b25") (.fillRect ctx gx gy gw gh)
                (doseq [i (range 0 (count cap) 4)]
                  (dot! (gpx (s (inc i)) (cap i)) 1.2 "rgba(255,143,112,0.7)"))
                (doseq [i (range 0 (- (count v) 4) 4)]
                  (line! (gpx (s i) (v i)) (gpx (s (+ i 4)) (v (+ i 4))) "#ffd479" 2))
                (dot! (gpx (s k) (v k)) 5 "#ffffff")
                (text! "speed against distance: coral the most each bend allows, gold the time-optimal profile" gx (- gy 8) "#9fb4d8" 11)
                (let [slowest (reduce min cap)]
                  (str "across in " (.toFixed total 1) " s at up to " vmax " m/s; driven at the tightest bend's " (.toFixed slowest 2)
                       " m/s throughout it would take " (.toFixed (/ (peek s) slowest) 1) " s"
                       "\nfull throttle to the top speed or the next bend's limit, full brakes into each bend -- a forward pass and a backward pass, the lesser of the two"))))
            (draw-systems! [cw ch]
              (let [{:keys [tick]} @model
                    unstable? (even? (.-seed controls))
                    A (if unstable? unstable stable)
                    P (control/lyapunov A [[1.0 0.0] [0.0 1.0]])
                    ;; left: trajectories and the level sets of V
                    s (min (/ cw 5.5) (/ ch 3.2)) lc [(* 0.27 cw) (* 0.5 ch)]
                    lpx (fn [[x y]] [(+ (lc 0) (* s x)) (- (lc 1) (* s y))])
                    n (min 80 (quot tick 2))
                    starts (for [a (range 8)] [(math/cos (* a (/ math/PI 4))) (math/sin (* a (/ math/PI 4)))])
                    step (fn [[x y]] [(+ (* (get-in A [0 0]) x) (* (get-in A [0 1]) y)) (+ (* (get-in A [1 0]) x) (* (get-in A [1 1]) y))])]
                (when P
                  (doseq [level [0.5 1.0 2.0 4.0 8.0]]
                    ;; x P x = level along each direction
                    (let [pts (for [i (range 73) :let [a (* 2 math/PI (/ i 72)) d [(math/cos a) (math/sin a)]
                                                       q (+ (* (get-in P [0 0]) (d 0) (d 0)) (* 2 (get-in P [0 1]) (d 0) (d 1)) (* (get-in P [1 1]) (d 1) (d 1)))
                                                       r (math/sqrt (/ level q))]]
                                [(* r (d 0)) (* r (d 1))])]
                      (doseq [[a b] (partition 2 1 pts)] (line! (lpx a) (lpx b) "rgba(160,180,210,0.35)" 1)))))
                (.save ctx)
                (.beginPath ctx) (.rect ctx 0 0 (* 0.55 cw) ch) (.clip ctx)
                (doseq [x0 starts]
                  (let [traj (take (inc n) (iterate step x0))]
                    (doseq [[a b] (partition 2 1 traj)
                            :when (< (max (abs (b 0)) (abs (b 1))) 2.6)]
                      (line! (lpx a) (lpx b) "#ffd479" 1.5))))
                (.restore ctx)
                (text! (if P "x' = A x: V = x P x (faint levels) falls along every trajectory" "x' = A x unstable: no P, the trajectories escape")
                       20 (- (lc 1) (* 1.5 s)) "#9fb4d8" 11)
                ;; right: the reachable sets
                (let [k (inc (mod (quot tick 6) 40))
                      cases [["a double integrator, pushed" [[1.0 0.1] [0.0 1.0]] [[0.0] [0.1]]]
                             ["two integrators, one pushed" [[1.0 0.0] [0.0 1.0]] [[0.1] [0.0]]]]]
                  (doseq [[i [label A B]] (map-indexed vector cases)
                          :let [c [(* 0.75 cw) (+ (* 0.3 ch) (* i 0.4 ch))]
                                poly (zonotope A B k)
                                ;; each axis scaled to the set's own extent, so its shape shows as it grows
                                ex (reduce max 1e-9 (map (fn [[x]] (abs x)) poly))
                                ey (reduce max 1e-9 (map (fn [[_ y]] (abs y)) poly))
                                rpx (fn [[x y]] [(+ (c 0) (* (/ (* 0.2 cw) ex) x)) (- (c 1) (* (/ (* 0.15 ch) ey) y))])
                                ok (control/controllable? A B)]]
                    (set! (.-fillStyle ctx) (if ok "rgba(125,255,160,0.25)" "rgba(255,143,112,0.25)"))
                    (set! (.-strokeStyle ctx) (if ok "#7dffa0" "#ff8f70")) (set! (.-lineWidth ctx) 2)
                    (.beginPath ctx)
                    (doseq [[j p] (map-indexed vector poly) :let [[x y] (rpx p)]] (if (zero? j) (.moveTo ctx x y) (.lineTo ctx x y)))
                    (.closePath ctx) (.fill ctx) (.stroke ctx)
                    (dot! c 3 "#ffffff")
                    (text! (str label ": " (if ok "controllable" "not controllable") ", rank "
                                (linear-systems/rank (control/controllability-matrix A B)))
                           (- (c 0) 150) (- (c 1) (* 0.17 ch)) (if ok "#7dffa0" "#ff8f70") 11)))
                (str (if P "stable: A^T P A - P = -I solved by summing (A^T)^k A^k" "unstable: the sum diverges, so there is no Lyapunov function")
                     "  |  right: every state reachable from rest in " (inc (mod (quot tick 6) 40)) " steps of a push of at most 1 (each axis scaled to fit)"
                     "\nseed odd: a stable system, even: an unstable one; a system is controllable when [B, AB, ...] has full rank -- the reachable set then fills the plane")))
            (draw-arms! [cw ch]
              (let [{:keys [arms tick]} @model
                    k (mod (* 5 tick) (count (first arms)))
                    s (min (/ cw 6) (/ ch 4.5)) base [(/ cw 2) (* 0.5 ch)]
                    at (fn [[t1 t2]] (let [e [(+ (base 0) (* s (math/cos t1))) (- (base 1) (* s (math/sin t1)))]]
                                       [e [(+ (e 0) (* s (math/cos (+ t1 t2)))) (- (e 1) (* s (math/sin (+ t1 t2))))]]))
                    draw (fn [q color w] (let [[e h] (at q)] (line! base e color w) (line! e h color w) (dot! h 5 color)))
                    errors (mapv (fn [states] (track/tracking-error (subvec states 0 (inc k)))) arms)]
                (draw (:ref ((first arms) k)) "rgba(255,255,255,0.18)" 14)
                (doseq [[states color] (map vector arms ["#ff8f70" "#7dc4ff" "#7dffa0"])] (draw (:q (states k)) color 3))
                (dot! base 7 "#c9d2e0")
                (str "largest error so far: open loop " (.toFixed (errors 0) 2) ", PD " (.toFixed (errors 1) 2) ", computed torque " (.toFixed (errors 2) 3) " rad"
                     "\nred: feedforward from a model a fifth too light; blue: PD alone; green: computed torque from the same wrong model; faint: where the arm should be")))
            (draw-pendulums! [cw ch]
              (let [{:keys [open mpc tick]} @model
                    k (min (quot tick 2) (dec (min (count open) (count mpc))))
                    s (min (/ cw 6) (/ ch 3))
                    draw (fn [cx th color label]
                           (let [pivot [cx (* 0.5 ch)] bob [(+ cx (* s (math/sin th))) (+ (* 0.5 ch) (* s (math/cos th)))]]
                             (line! pivot bob color 4) (dot! bob 10 color) (dot! pivot 5 "#c9d2e0")
                             (text! label (- cx 60) (- (* 0.5 ch) s 30) color 13)))]
                (draw (* 0.3 cw) (first (open k)) "#ff8f70" "the plan, open loop")
                (draw (* 0.7 cw) (first (mpc k)) "#7dffa0" "MPC, replanning")
                (str "step " k "  |  angles from hanging: open loop " (.toFixed (first (open k)) 2) ", MPC " (.toFixed (first (mpc k)) 2) " (upright is " (.toFixed math/PI 2) ")"
                     "\nboth pendulums 30% heavier than the model the iLQR plan was made for; MPC plans 40 steps ahead from where it actually is, every step, and keeps the first")))
            (draw-peg! [cw ch]
              (let [{:keys [trials tick]} @model
                    s (min (/ cw 3.2) (/ ch 3.2)) ox (/ cw 2) oy (* 0.68 ch)
                    px (fn [[x y]] [(+ ox (* s x)) (- oy (* s y))])
                    w (/ (:width hole) 2) d (:depth hole)
                    k (quot tick 2)
                    shown (count (take-while (fn [i] (> k (reduce + (map #(count (:path (trials %))) (range (inc i)))))) (range (count trials))))
                    done (subvec trials 0 (min shown (count trials)))]
                (set! (.-fillStyle ctx) "#3a4458")
                (.fillRect ctx 0 oy (- ox (* s w)) (- ch oy))
                (.fillRect ctx (+ ox (* s w)) oy cw (- ch oy))
                (.fillRect ctx (- ox (* s w)) (+ oy (* s d)) (* 2 s w) ch)
                (doseq [{:keys [strategy path]} done]
                  (doseq [[a b] (partition 2 1 (take-nth 3 path))]
                    (line! (px a) (px b) (if (= strategy :direct) "rgba(255,143,112,0.35)" "rgba(125,255,160,0.35)") 1.5)))
                (when (< shown (count trials))
                  (let [{:keys [strategy path]} (trials shown)
                        before (reduce + (map #(count (:path (trials %))) (range shown)))
                        upto (min (count path) (max 1 (- k before)))]
                    (doseq [[a b] (partition 2 1 (take upto path))] (line! (px a) (px b) (if (= strategy :direct) "#ff8f70" "#7dffa0") 2.5))
                    (dot! (px (nth path (dec upto))) 6 "#ffffff")))
                (let [score (fn [st] (let [ts (filter #(= st (:strategy %)) done)]
                                       (str (count (filter (fn [t] (guarded/at-bottom? hole (peek (:path t)))) ts)) " of " (count ts))))]
                  (str "reached the bottom: straight down " (score :direct) ", the two guarded moves " (score :guarded)
                       "\nevery motion strays up to 20 degrees; red: straight down, stopping at contact; green: down and left to the floor, then down and right, sliding along it into the hole"))))
            (draw-hexapod! [cw ch]
              (let [{:keys [states rocks tick]} @model
                    k (min (quot tick 3) (dec (count states)))
                    {:keys [legs strides]} (states k)
                    s (min (/ cw 14) (/ ch 9))
                    oy (* 0.62 ch)
                    stride-px (* 0.8 s)
                    leg-x {:rf 3.0 :rm 1.5 :rr 0.0 :lf 3.0 :lm 1.5 :lr 0.0}
                    body-x (+ (* 0.2 cw) (* 0.0 strides))]
                ;; the ground and the rocks coming up, scrolling as it walks
                (line! [0 oy] [cw oy] "#5b6884" 2)
                (doseq [[[leg st] h] rocks
                        :let [x (+ body-x (* s (leg-x leg)) (* stride-px (- st strides)) stride-px)]
                        :when (< -20 x (+ cw 20))]
                  (set! (.-fillStyle ctx) "#6b5a44") (.fillRect ctx (- x 6) (- oy (* 4 s h)) 12 (* 4 s h)))
                (set! (.-fillStyle ctx) "#3a4458") (.fillRect ctx body-x (- oy (* 2.2 s)) (* 3.4 s) (* 0.6 s))
                (doseq [[leg {:keys [state h retries]}] legs
                        :let [x (+ body-x (* s (leg-x leg)) (if (#{:lf :lm :lr} leg) (* 0.2 s) 0.0))
                              lift (case state (:lift :swing) (* 4 s h) :stuck (* 4 s h) 0.0)
                              forward (case state :swing (* 0.5 stride-px) :down (* 0.8 stride-px) 0.0)
                              foot [(+ x forward) (- oy lift)]
                              color (cond (= state :stuck) "#ff8f70" (reactive/grounded? {:state state}) "#9fb4d8" :else "#ffd479")]]
                  (line! [x (- oy (* 1.9 s))] foot color 3)
                  (when (pos? retries) (text! (str retries) (+ x 4) (- oy (* 2.4 s)) "#ff8f70" 11)))
                (str strides " strides, tick " (:ticks (states k)) "  |  standing legs: " (count (filter reactive/grounded? (vals legs)))
                     "\ntripod gait: three legs swing (gold) while three stand (blue); a swing blocked by a rock (red) retracts, lifts 0.05 higher and tries again -- the numbers count the retries")))
            (draw-people! [cw ch]
              (let [{:keys [walk beliefs tick]} @model
                    k (min (quot tick 20) (dec (count beliefs)))
                    s (min (/ (- cw 300) 11) (/ (- ch 160) 11))
                    px (fn [[x y]] [(+ 40 (* s x)) (- (- ch 80) (* s y))])
                    b (beliefs k)]
                (set! (.-strokeStyle ctx) "#5b6884") (set! (.-lineWidth ctx) 1)
                (.strokeRect ctx 40 (- (- ch 80) (* s 10)) (* s 10) (* s 10))
                (doseq [[g p] room-goals] (dot! (px p) 9 "#ffd479") (text! g (+ ((px p) 0) 12) ((px p) 1) "#ffd479" 12))
                (doseq [[a c] (partition 2 1 (take (inc k) walk))] (line! (px a) (px c) "#7dc4ff" 3))
                (let [here (walk k)]
                  (dot! (px here) 7 "#ffffff")
                  (doseq [[[dx dy] p] (humans/predict-step room-goals b here steps8 2.0)]
                    (line! (px here) (px [(+ (here 0) (* 3 dx p)) (+ (here 1) (* 3 dy p))]) "rgba(125,255,160,0.8)" 3)))
                (doseq [[i [g p]] (map-indexed vector (sort-by key b))]
                  (let [x (- cw 260) y (+ 240 (* 44 i))]
                    (text! g x y "#c9d2e0" 13)
                    (set! (.-fillStyle ctx) "#7dffa0") (.fillRect ctx (+ x 80) (- y 12) (* 140 p) 14)
                    (text! (.toFixed p 2) (+ x 84) (+ y 18) "#9fb4d8" 11)))
                (str "step " k "  |  the person is assumed noisily rational: the better a step for a goal, the likelier"
                     "\neach step taken updates the belief over goals (bars); green arrows: the predicted next step, longer the likelier")))
            (draw-terrain! [cw ch]
              (let [{:keys [frames models tick]} @model
                    k (if (< (mod tick 400) 200) 0 1)
                    {:keys [image]} (frames k)
                    mdl (models k)
                    test (frames 0)
                    s (min (/ (- cw 60) 80) (/ (- ch 180) 30))
                    blit (fn [x0 f]
                           (doseq [y (range 30) x (range 40)]
                             (let [[r g b] (f x y)]
                               (set! (.-fillStyle ctx) (str "rgb(" (* 255 r) "," (* 255 g) "," (* 255 b) ")"))
                               (.fillRect ctx (+ x0 (* s x)) (+ 110 (* s y)) (+ s 0.5) (+ s 0.5)))))]
                (blit 20 (fn [x y] ((:px image) (+ (* y 40) x))))
                (set! (.-strokeStyle ctx) "#ffffff") (set! (.-lineWidth ctx) 2)
                (.strokeRect ctx (+ 20 (* s 15)) (+ 110 (* s 25)) (* s 10) (* s 5))
                (text! (if (zero? k) "on the road: the laser sweeps the white box" "off the road: grass in the box now") 20 98 "#9fb4d8" 12)
                (blit (+ 40 (* 40 s)) (fn [x y] (let [c ((:px (:image test)) (+ (* y 40) x))]
                                                  (if (terrain/drivable? mdl c) [0.95 0.45 0.75] (mapv #(* 0.5 %) c)))))
                (text! "the camera's verdict on the road scene: pink drivable" (+ 40 (* 40 s)) 98 "#9fb4d8" 12)
                (str (if (zero? k) "taught by the road just ahead, the camera finds road drivable everywhere -- and not grass"
                         "having driven onto grass, the laser found it flat: now grass is drivable too; rocks never")
                     "\nself-supervised: the short-range laser labels its own training data for a mixture of Gaussians over color, which then sees as far as the camera")))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))]
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "method" #js ["computed torque" "MPC" "peg in hole" "hexapod" "predicting people" "terrain"
                                                    "time-optimal speed" "controllability and Lyapunov"]) restart!)
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.add gui controls "speed" 1 10 1)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "robot-control")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
