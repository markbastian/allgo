(ns allgo.simulation.courses
  "Places to ride the motorcycle, each with a rider who knows the way.

  A course is static boxes and nothing else -- a jump is a box tilted to
  the ramp's slope, a tabletop is a box, a whoop is a box turned onto its
  edge and half buried -- because boxes are what the tires' contact
  already understands, and a torus meeting a tilted box is a tire meeting
  a ramp. Each box carries a `:kind`, which the physics ignores and the
  drawing uses to color it.

  Each course also has an autopilot: the lean and throttle a rider who
  knows the course would ask for. It chooses only what the keys would;
  the balancing and the countersteering are still
  `allgo.simulation.motorcycle`'s rider, which is why a course can be
  ridden by hand exactly as it is ridden by itself. An autopilot is
  `(fn [pose telemetry {:keys [cfg max-lean]}] controls)`.

    :loop        a 35 m circle over a ring of curbs, with a kicker off to
                 the side
    :excitebike  a straight lane between barriers: jumps, a tabletop, a
                 double, whoops and a step-up, ridden flat out for time
    :trials      logs, steps and a rock garden, taken at a walk"
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.physics.rigid :as rigid]
            [allgo.simulation.motorcycle :as m]))

;; ---------------------------------------------------------------------------
;; Shapes

(defn- tag [body kind] (assoc body :kind kind))

(defn slope
  "A ramp whose riding surface runs straight from `[x0 y0]` to `[x1 y1]`
  along +x, `width` across the lane and `thick` deep beneath its surface,
  so that a thick one reads as solid ground rather than a plank on legs."
  ([from to width] (slope from to width 0.4))
  ([[x0 y0] [x1 y1] width thick]
   (let [dx (- (double x1) (double x0))
         dy (- (double y1) (double y0))
         len (Math/sqrt (+ (* dx dx) (* dy dy)))
         a (Math/atan2 dy dx)
         normal [(- (Math/sin a)) (Math/cos a) 0.0]
         mid [(* 0.5 (+ (double x0) (double x1))) (* 0.5 (+ (double y0) (double y1))) 0.0]]
     (tag (rigid/box {:pos (v/add mid (v/scale normal (* -0.5 (double thick))))
                      :size [len thick width]
                      :rot (q/from-axis-angle [0.0 0.0 1.0] a)})
          :ramp))))

(defn block
  "A flat top at height `h` from `x0` to `x1`, solid down to the ground:
  a tabletop, a step, a plateau."
  [x0 x1 h width]
  (let [h (double h)]
    (tag (rigid/box {:pos [(* 0.5 (+ (double x0) (double x1))) (* 0.5 (- h 0.5)) 0.0]
                     :size [(- (double x1) (double x0)) (+ h 0.5) width]})
         :ramp)))

(defn whoop
  "One of a run of small humps, `h` high and `span` long: a gentle slope
  up and one back down. Steep faces are a wall at speed, not a whoop."
  [x h span width]
  (let [x (double x) half (* 0.5 (double span))]
    [(tag (slope [x 0.0] [(+ x half) h] width 0.3) :whoop)
     (tag (slope [(+ x half) h] [(+ x (double span)) 0.0] width 0.3) :whoop)]))

(def ^:private g 9.81)

(defn flight
  "Where a jump's landing belongs, for a bike reaching the foot of a ramp
  of `run` and `h` at `speed`: the gap from the lip to the top of the
  landing slope, and the landing's length.

  The bike leaves the lip at whatever speed the climb has left it -- its
  kinetic energy less `g h` -- along the ramp's angle, and flies a
  parabola. A landing is right when it slopes the way the bike is falling
  as it arrives, so the suspension takes the rest of the speed along the
  slope rather than all of it at once; and the gap is set so the bike
  touches down `touch` of the way down it. That is how a real jump is
  built, and why a jump taken much faster or slower than it was built
  for goes wrong."
  [run h speed touch]
  (let [run (double run) h (double h) touch (double touch)
        angle (Math/atan2 h run)
        v (Math/sqrt (max 1.0 (- (* (double speed) (double speed)) (* 2.0 g h))))
        vx (* v (Math/cos angle))
        vy (* v (Math/sin angle))
        ;; When the arc comes back down to the touchdown height.
        drop (* touch h)
        t (/ (+ vy (Math/sqrt (+ (* vy vy) (* 2.0 g drop)))) g)
        descent (/ (- (* g t) vy) vx)
        ;; A little shallower than the fall, which makes the landing half
        ;; as long again: a bigger place to put the wheels, for a slightly
        ;; firmer arrival. A landing exactly as steep as the arc is only
        ;; right for exactly the right speed.
        landing (/ h (* 0.7 (max 0.15 descent)))]
    ;; The arc is flown from the lip, and it is the center of mass's. The
    ;; center is half a wheelbase past the lip when the rear tire leaves
    ;; it and the front tire half a wheelbase ahead of the center, so the
    ;; tires arrive a wheelbase later than this says -- which the landing
    ;; is long enough to take, touching down short of halfway.
    {:gap (max 0.0 (- (* vx t) (* touch landing)))
     :landing landing}))

(def ^:private touch-down
  "How far down its landing a jump is built to be landed: past halfway, so
  that a bike arriving a little slow still lands on the slope and not on
  its top edge -- which is where a short jump goes wrong, and badly."
  0.4)

(defn jump
  "An up ramp from the ground at `x` rising `h` over `run`, and a landing
  laid out by `flight` for a bike arriving at `speed`."
  [x run h speed width]
  (let [x (double x) run (double run)
        {:keys [gap landing]} (flight run h speed touch-down)
        thick (+ (double h) 1.0)]
    [(slope [x 0.0] [(+ x run) h] width thick)
     (slope [(+ x run gap) h] [(+ x run gap landing) 0.0] width thick)]))

(defn tabletop
  "A jump whose gap is filled in: fly it at `speed` and land on the down
  slope, come up short and land on the top."
  [x run h speed width]
  (let [x (double x) run (double run)
        {:keys [gap landing]} (flight run h speed touch-down)
        thick (+ (double h) 1.0)]
    [(slope [x 0.0] [(+ x run) h] width thick)
     (block (+ x run) (+ x run gap) h width)
     (slope [(+ x run gap) h] [(+ x run gap landing) 0.0] width thick)]))

(defn- lane-steer
  "The lean that keeps a bike heading along +x on the line z = 0.

  Pure pursuit: aim at a point on the line a speed's worth of distance
  ahead, and ask for the lean that the arc through it needs -- the
  steady-turn lean `atan(v^2 k / g)` for that arc's curvature `k`. A
  lean proportional to the heading error does the same job at speed and
  weaves at a walk, where the bike turns too slowly for it; pursuit asks
  for a turn the bike can actually make at the speed it is going."
  [pose speed]
  (let [{:keys [pos rot]} (:base pose)
        v (max 1.5 (double speed))
        ahead (max 4.0 (* 1.2 v))
        fwd (let [f (q/rotate rot [1.0 0.0 0.0])] (v/normalize [(nth f 0) 0.0 (nth f 2)]))
        want (v/normalize [ahead 0.0 (- (double (nth pos 2)))])
        err (Math/atan2 (nth (v/cross fwd want) 1) (v/dot fwd want))
        curvature (/ (* 2.0 (Math/sin err)) ahead)]
    (max -0.35 (min 0.35 (Math/atan (/ (* v v curvature) 9.81))))))

(defn- cruise
  "Throttle to hold `speed`, and the brake if well over it."
  [speed target]
  (let [e (- (double target) (double speed))]
    {:throttle (max 0.0 (min 1.0 (+ 0.1 (* 0.3 e))))
     :brake (if (< e -3.0) (min 1.0 (* 0.2 (- -3.0 e))) 0.0)}))

;; ---------------------------------------------------------------------------
;; The loop

(def ^:private loop-radius 35.0)
(def ^:private loop-center [0.0 0.0 (- loop-radius)])

(defn- on-loop [theta]
  (let [theta (double theta)]
    [[(* loop-radius (Math/sin theta)) 0.0 (- (* loop-radius (Math/cos theta)) loop-radius)]
     theta]))

(defn- loop-obstacles []
  (conj (vec (for [theta [0.55 1.5 2.5 3.4 4.5 5.5]
                   :let [[[x _ z] across] (on-loop theta)]]
               (tag (rigid/box {:pos [x 0.04 z] :size [0.45 0.08 4.0]
                                :rot (q/from-axis-angle [0.0 1.0 0.0] across)})
                    :curb)))
        ;; A kicker off to the right of the start, for riding by hand.
        (tag (rigid/box {:pos [0.0 0.35 22.0] :size [5.0 0.3 4.0]
                         :rot (q/from-axis-angle [0.0 0.0 1.0] (/ (* 9.0 Math/PI) 180.0))})
             :ramp)))

(defn- loop-autopilot
  "Lap the loop: the lean a steady turn of its radius needs at this
  speed, plus whatever turns the bike toward where it should be heading."
  [pose {:keys [speed]} {:keys [max-lean] :or {max-lean 0.7}}]
  (let [{:keys [pos rot]} (:base pose)
        fwd (let [f (q/rotate rot [1.0 0.0 0.0])] (v/normalize [(nth f 0) 0.0 (nth f 2)]))
        r (let [d (v/sub pos loop-center)] [(nth d 0) 0.0 (nth d 2)])
        dist (max 1e-6 (v/length r))
        out (v/scale r (/ 1.0 dist))
        along [(nth out 2) 0.0 (- (double (nth out 0)))]
        want (v/normalize (v/sub along (v/scale out (* 0.08 (- dist loop-radius)))))
        err (Math/atan2 (nth (v/cross fwd want) 1) (v/dot fwd want))
        steady (Math/atan (/ (* (double speed) (double speed)) (* 9.81 loop-radius)))
        ;; Not quite as hard over as the keys allow: leaned right down, a
        ;; rider's inside boot is near enough the ground to catch a curb.
        cap (min 0.5 (double max-lean))]
    (merge (cruise speed 11.0)
           {:lean (max (- cap) (min cap (+ steady (* 1.4 err))))})))

;; ---------------------------------------------------------------------------
;; Excitebike

(def ^:private lane 7.0)

(def excitebike-finish 320.0)

(defn- barriers
  "Hay bales down both sides of the lane, in lengths."
  [length]
  (vec (for [x (range 0.0 length 50.0)
             z [(- (+ (* 0.5 lane) 0.3)) (+ (* 0.5 lane) 0.3)]]
         (tag (rigid/box {:pos [(+ (double x) 25.0) 0.35 z] :size [50.0 0.7 0.5]})
              :barrier))))

(defn- whoops
  "A run of humps from `from` to `to`, each `span` long and nose to tail."
  [from to h span]
  (vec (mapcat #(whoop % h span lane) (range (double from) (double to) (double span)))))

(def excitebike-speed
  "What the lane is built for: every jump's landing is laid out by
  `flight` for a bike arriving at this."
  15.0)

(defn- excitebike-obstacles []
  (let [v excitebike-speed]
    (-> (barriers (+ excitebike-finish 40.0))
        ;; A first jump to get the feel of it.
        (into (jump 40.0 6.0 1.1 v lane))
        (into (whoops 72.0 90.0 0.18 2.6))
        ;; A tabletop: up, across the top, and down.
        (into (tabletop 100.0 6.0 1.5 v lane))
        ;; The big one: a double.
        (into (jump 145.0 7.0 2.0 v lane))
        (into (whoops 185.0 205.0 0.22 2.8))
        ;; A step-up onto a plateau, and a jump off the far end of it.
        (conj (slope [215.0 0.0] [220.0 1.0] lane 2.0)
              (block 220.0 240.0 1.0 lane))
        (conj (slope [240.0 1.0] [244.0 1.8] lane 2.8))
        (into (let [{:keys [gap landing]} (flight 4.0 0.8 (Math/sqrt (- (* v v) (* 2.0 g 1.0))) touch-down)]
                [(slope [(+ 244.0 gap) 1.8] [(+ 244.0 gap landing) 0.0] lane 2.8)]))
        ;; One more to finish on. Not too steep: a short, steep lip kicks
        ;; the rear of the bike up as it leaves, and the bike arrives
        ;; standing on its back wheel, if it arrives at all.
        (into (jump 275.0 6.0 1.4 v lane)))))

(defn airborne?
  "Whether both wheels are off the ground, as a rider feels it: both ends
  of the suspension have run out to full extension, which they only do
  with nothing pressing on them."
  [cfg {:keys [rear-travel front-travel]}]
  (and (> (double rear-travel) (- (double (second (:rear-travel cfg))) 0.02))
       (> (double front-travel) (- (double (second (:front-travel cfg))) 0.01))))

(def ^:private landing-pitch
  "How far nose down a rider wants the bike on the way down: about the
  slope of a landing, so both wheels meet it together."
  -0.2)

(defn- fly
  "Throttle and brake in the air, which is the only way a rider can pitch
  a bike there. Spinning the rear wheel up pushes the bike the other way
  about its axle, nose up; braking it does the reverse. Excitebike made a
  game of this, and it is the same physics."
  [pitch]
  (let [e (- (double pitch) landing-pitch)]
    ;; Gently: a wheel's spin is a small lever on a whole motorcycle, and
    ;; one still being corrected when the bike lands is worse than a bike
    ;; a little off level.
    (cond (> e 0.1) {:throttle 0.0 :brake (min 0.6 (* 2.0 e))}
          (< e -0.1) {:throttle (min 1.0 (* -2.0 e)) :brake 0.0}
          :else {:throttle 0.15 :brake 0.0})))

(defn- excitebike-autopilot
  "Flat out down the middle of the lane, and in the air, leveled for the
  landing with the throttle and brake."
  [pose {:keys [speed pitch] :as t} {:keys [cfg]}]
  (merge (if (airborne? cfg t)
           (fly pitch)
           ;; A touch over what the jumps are built for: the climb up each
           ;; ramp costs more than the energy sum allows for.
           (cruise speed (+ excitebike-speed 0.5)))
         {:lean (lane-steer pose speed)}))

;; ---------------------------------------------------------------------------
;; Trials

(def trials-finish 78.0)

(defn- log
  "A log across the way, `d` thick: a square beam turned on its edge,
  which the tire meets the way it would meet a round one."
  [x d]
  (tag (rigid/box {:pos [(double x) (* 0.5 (double d)) 0.0]
                   :size [(* (double d) (Math/sqrt 0.5)) (* (double d) (Math/sqrt 0.5)) 4.0]
                   :rot (q/from-axis-angle [0.0 0.0 1.0] (/ Math/PI 4.0))})
       :log))

(defn- rocks
  "A rock garden from `x0` to `x1`: small blocks of stone, scattered
  and tilted, the same every time."
  [x0 x1]
  (let [next-r (fn [s] (mod (* 16807 (long s)) 2147483647))
        rs (rest (iterate next-r 42))
        u (fn [i] (/ (double (nth rs i)) 2147483647.0))]
    (vec (for [k (range 22)
               :let [i (* 6 k)
                     x (+ (double x0) (* (u i) (- (double x1) (double x0))))
                     z (- (* 3.0 (u (+ i 1))) 1.5)
                     s (+ 0.2 (* 0.2 (u (+ i 2))))
                     h (+ 0.05 (* 0.08 (u (+ i 3))))]]
           (tag (rigid/box {:pos [x (* 0.4 h) z] :size [s h s]
                            :rot (q/mul (q/from-axis-angle [0.0 1.0 0.0] (* 3.0 (u (+ i 4))))
                                        (q/from-axis-angle [1.0 0.0 0.0] (* 0.3 (- (u (+ i 5)) 0.5))))})
                :rock)))))

(defn- trials-obstacles []
  (-> [(log 8.0 0.22) (log 13.0 0.3)
       ;; Two steps up, a step down, and off the end.
       (block 20.0 26.0 0.18 4.0)
       (block 26.0 33.0 0.34 4.0)
       (block 33.0 39.0 0.16 4.0)]
      (into (rocks 45.0 58.0))
      (conj (log 64.0 0.2) (log 65.8 0.2) (log 71.0 0.26))))

(defn- trials-autopilot
  "At a walk, straight through."
  [pose {:keys [speed]} _]
  (merge (cruise speed 4.5) {:lean (lane-steer pose speed)}))

;; ---------------------------------------------------------------------------

(def courses
  "Every course, in the order they are offered."
  [{:id :loop :name "Loop"
    :obstacles loop-obstacles
    :start {:speed 10.0 :heading 0.0 :pos [0.0 0.0]}
    :autopilot loop-autopilot}
   {:id :excitebike :name "Excitebike"
    :obstacles excitebike-obstacles
    :start {:speed 12.0 :heading 0.0 :pos [0.0 0.0]}
    :autopilot excitebike-autopilot
    :finish excitebike-finish}
   {:id :trials :name "Trials"
    :obstacles trials-obstacles
    :start {:speed 2.0 :heading 0.0 :pos [0.0 0.0]}
    :autopilot trials-autopilot
    :finish trials-finish}])

(defn course [id] (some #(when (= id (:id %)) %) courses))

(defn start-pose
  "Where and how fast the bike begins `course`, overriding its speed with
  `speed` if given."
  ([cfg crs] (start-pose cfg crs nil))
  ([cfg {:keys [start]} speed]
   (m/start-pose cfg (or speed (:speed start)) (:heading start) (:pos start))))

(defn scene
  "`course`'s obstacles on a floor, with the bike at its start."
  [cfg crs opts]
  (m/scene cfg (into [(m/ground 4000.0)] ((:obstacles crs)))
           (start-pose cfg crs (:speed opts))
           opts))

(defn finished?
  "Whether the bike has crossed `course`'s finish line."
  [crs pose]
  (when-let [x (:finish crs)]
    (> (double (first (:pos (:base pose)))) (double x))))
