(ns allgo.simulation.crossbows
  "A game after Lakeside's *Crossbows and Catapults* (1983): two castles
  facing each other across a board, and discs flung at them until one
  side's tower goes over.

  ## The pieces

  Counted off a photograph of a complete set, per side:

    tower        1   the round keep; knock the other one over to win
    wall         5   crenellated sections
    flag         5   pennants on stands
    warrior      4   Vikings on one side, Barbarians on the other
    disc        10   the ammunition, called caroms here
    catapult     1   flings a disc through the air
    crossbow     1   shoots a disc along the ground, like a puck

  The wall sections overlap in the photograph and are the least certain
  count. The rules are not Lakeside's -- the rule sheet was not to hand --
  but a plain version of the idea: take turns, one disc a turn, and the
  first to topple the other side's tower wins. Everything else on the
  board is there to be in the way.

  ## The board

  Built for `allgo.physics.solver`, at about ten times the toy's size so
  that a disc is as big as the solver's contacts are comfortable with.
  Each castle is laid out in its own frame -- `x` to the right of someone
  standing behind it looking at the enemy, `d` back from the middle of
  the board -- so the two are the same layout turned round. The five
  wall pieces make a gate: two either side of a gap a disc can slide
  through, and the fifth laid across the top of the inner two, which is
  why a crossbow shot has a way in at all.

  A disc is a thin square box, not a ball. The engine has no cylinder,
  and a ball rolls where a puck should slide.

  ## A turn

  `fire` puts a disc in the air or on the ground; `advance` steps the
  world and, once everything has stopped moving or eight seconds have
  passed, settles the turn: a tower tipped past forty-five degrees loses,
  and otherwise it is the other side's shot. Both out of discs with both
  towers up is a draw."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.math :as am]
            [allgo.physics.rigid :as rigid]
            [allgo.physics.solver :as solver]
            [clojure.math :as math]))

(def sides [:vikings :barbarians])

(defn other [side] (if (= side :vikings) :barbarians :vikings))

(def pieces
  "How many of each, per side, as counted from the set."
  {:tower 1 :wall 5 :flag 5 :warrior 4 :disc 10 :catapult 1 :crossbow 1})

(def gravity [0.0 -10.0 0.0])

;; ---------------------------------------------------------------------------
;; Frames

(defn- frame
  "A side's axes in the world: `right`, and `back` -- away from the enemy.
  The Vikings hold the +z end of the board, the Barbarians the -z end."
  [side]
  (if (= side :vikings)
    {:right [1.0 0.0 0.0] :back [0.0 0.0 1.0]}
    {:right [-1.0 0.0 0.0] :back [0.0 0.0 -1.0]}))

(defn- place
  "The world point `x` to the right, `y` up and `d` back of the middle of
  the board, for `side`."
  [side x y d]
  (let [{:keys [right back]} (frame side)]
    (v/add (v/add (v/scale right x) [0.0 y 0.0]) (v/scale back d))))

(defn- facing
  "The rotation that turns a piece built facing -z to face `side`'s enemy."
  [side]
  (if (= side :vikings) q/identity-q (q/from-axis-angle [0.0 1.0 0.0] math/PI)))

;; ---------------------------------------------------------------------------
;; The board

(def ^:private tower-size
  "Slimmer than the toy's, which is squat enough that no disc could tip
  it: at 1.4 by 2.6 a perfect catapult shot rocked it fifteen degrees and
  it sat back down. At 1.3 by 2.8 a clean hit high on the face topples
  it, a glancing one rocks it and it sits back, and a crossbow shot at
  its foot barely moves it -- the crossbow is for clearing the way.

  Hits do not add up. A rigid tower that rocks and settles has lost
  nothing, so the game is decided by one good shot, as the toy's is."
  [1.3 2.8 1.3])
(def ^:private wall-size [1.8 1.0 0.5])
(def ^:private flag-size [0.14 1.4 0.14])
(def ^:private warrior-size [0.35 0.7 0.35])
(def disc-size
  "A disc's collider: thin, and square, which a round disc fits inside."
  [0.56 0.14 0.56])

(def layout
  "Where each piece stands, in a side's own frame: `[x d]`, and for the
  wall a height and whether it is the lintel over the gate."
  {:tower [[0.0 7.5]]
   ;; The inner two leave a gate 1.2m wide; the fifth lies across them
   :wall [[-3.4 5.0 0] [-1.5 5.0 0] [1.5 5.0 0] [3.4 5.0 0]
          [0.0 5.0 1]]
   :flag [[-2.4 6.3] [2.4 6.3] [-1.8 8.6] [1.8 8.6] [0.0 9.4]]
   :warrior [[-3.0 4.2] [-0.9 4.2] [0.9 4.2] [3.0 4.2]]
   :catapult [[3.0 9.8]]
   :crossbow [[0.0 2.6]]})

(def launchers
  "Where each weapon releases its disc, in a side's own frame `[x y d]`,
  and how fast it can throw: meters a second, least to most."
  {:catapult {:at [3.0 1.6 9.8] :speed [8.0 20.0] :default-speed 15.0 :default-pitch 35.0}
   :crossbow {:at [0.0 0.08 2.0] :speed [6.0 22.0] :default-speed 15.0 :default-pitch 0.0}})

(defn- piece [side kind size density pos rot]
  (-> (rigid/box {:pos pos :rot rot :size size :density density})
      (assoc :side side :kind kind)))

(defn- castle [side]
  (let [rot (facing side)
        [_ ty _] tower-size
        [_ wy _] wall-size
        [_ fy _] flag-size
        [_ gy _] warrior-size]
    (concat
     (for [[x d] (:tower layout)]
       ;; Light for its size: the toy's tower is a hollow shell.
       (piece side :tower tower-size 0.2 (place side x (* 0.5 ty) d) rot))
     (for [[x d course] (:wall layout)]
       (piece side :wall wall-size 0.6 (place side x (+ (* 0.5 wy) (* course wy)) d) rot))
     (for [[x d] (:flag layout)]
       (piece side :flag flag-size 0.5 (place side x (* 0.5 fy) d) rot))
     (for [[x d] (:warrior layout)]
       (piece side :warrior warrior-size 0.8 (place side x (* 0.5 gy) d) rot)))))

(defn engine-pose
  "Where `side`'s catapult or crossbow stands, and which way it faces.
  They are drawn, not simulated: a disc passes through them, which keeps
  the crossbow from blocking its own gate."
  [side weapon]
  (let [[x d] (first (layout weapon))]
    {:pos (place side x 0.0 d) :rot (facing side)}))

(defn board
  "Every body on the board, the ground first."
  []
  (into [(assoc (rigid/box {:pos [0.0 -0.5 0.0] :size [16.0 1.0 26.0]}) :kind :ground)]
        (mapcat castle sides)))

(def solver-options
  {:solver :tgs :iterations 8 :substeps 4 :friction 0.6 :gravity gravity})

;; ---------------------------------------------------------------------------
;; Shots

(defn launch-point [side weapon]
  (let [[x y d] (:at (launchers weapon))] (place side x y d)))

(defn launch-velocity
  "`speed` along the aim of `side`'s `weapon`: `yaw` degrees to the right
  of straight at the enemy, `pitch` degrees up. A crossbow has no pitch;
  its disc leaves along the ground."
  [side weapon yaw pitch speed]
  (let [{:keys [right back]} (frame side)
        fwd (v/scale back -1.0)
        y (math/to-radians (double yaw))
        p (if (= weapon :crossbow) 0.0 (math/to-radians (double pitch)))
        c (math/cos p)]
    (v/scale (v/add (v/add (v/scale fwd (* c (math/cos y)))
                           (v/scale right (* c (math/sin y))))
                    [0.0 (math/sin p) 0.0])
             (double speed))))

(defn flight-path
  "Where a shot would go if it hit nothing, sampled until it reaches the
  ground or three seconds pass. A crossbow's is a straight line along the
  ground, shortened by friction."
  [side weapon yaw pitch speed]
  (let [p0 (launch-point side weapon)
        v0 (launch-velocity side weapon yaw pitch speed)]
    (if (= weapon :crossbow)
      (let [stop (/ (* (double speed) (double speed)) (* 2.0 0.6 10.0))
            dir (v/normalize v0)]
        (for [i (range 0 31)] (v/add p0 (v/scale dir (* stop (/ i 30.0))))))
      (for [i (range 0 91)
            :let [t (* i (/ 3.0 90.0))
                  p (v/add (v/add p0 (v/scale v0 t)) (v/scale gravity (* 0.5 t t)))]
            :while (>= (nth p 1) 0.0)]
        p))))

(defn- disc [side weapon yaw pitch speed]
  (let [p0 (launch-point side weapon)
        v0 (launch-velocity side weapon yaw pitch speed)
        ;; A catapult throws a disc turning over and over. A crossbow's
        ;; slides flat and square on: spinning, its corners reach 0.79m
        ;; across, and it would catch the sides of the gate.
        spin (if (= weapon :crossbow) v/zero (v/scale (v/normalize (v/cross v0 [0.0 1.0 0.0])) 8.0))]
    (-> (rigid/box {:pos (if (= weapon :crossbow)
                           (assoc p0 1 (+ 0.001 (* 0.5 (nth disc-size 1))))
                           p0)
                    :size disc-size :density 4.5 :vel v0 :omega spin})
        (assoc :side side :kind :disc))))

;; ---------------------------------------------------------------------------
;; Standing and fallen

(defn up-tilt
  "How far a body has turned from upright, in degrees."
  [b]
  (let [[_ uy _] (q/rotate (:rot b) [0.0 1.0 0.0])]
    (math/to-degrees (math/acos (am/clamp (double uy) -1.0 1.0)))))

(defn down?
  "Tipped past forty-five degrees, or off the board."
  [b]
  (or (> (up-tilt b) 45.0) (< (double (nth (:pos b) 1)) -2.0)))

(defn tally
  "For each side, how many of each kind of its own pieces are down."
  [world]
  (reduce (fn [m b]
            (if (and (:side b) (#{:tower :wall :flag :warrior} (:kind b)) (down? b))
              (update-in m [(:side b) (:kind b)] (fnil inc 0))
              m))
          {:vikings {} :barbarians {}}
          (:bodies world)))

;; ---------------------------------------------------------------------------
;; The game

(defn new-game
  "The board settled and asleep, the Vikings to shoot first."
  []
  (let [w (solver/settled (solver/world (board) solver-options))]
    {:world (update w :bodies (fn [bs] (mapv #(if (rigid/static? %) % (rigid/asleep %)) bs)))
     :turn :vikings
     :discs {:vikings (:disc pieces) :barbarians (:disc pieces)}
     :phase :aiming
     :winner nil}))

(defn fire
  "The side to play shoots. Ignored unless it is waiting to aim and has a
  disc left."
  [{:keys [turn phase discs] :as game} weapon yaw pitch speed]
  (if (and (= phase :aiming) (pos? (long (discs turn))))
    (-> game
        (update-in [:world :bodies] conj (disc turn weapon yaw pitch speed))
        (update-in [:discs turn] dec)
        (assoc :phase :flying :clock 0.0 :calm 0))
    game))

(defn- calm? [world]
  (every? (fn [b] (or (rigid/inert? b)
                      (and (< (v/length (:vel b)) 0.05) (< (v/length (:omega b)) 0.2))))
          (:bodies world)))

(defn- resolve-turn [{:keys [world turn discs] :as game}]
  (let [fallen (tally world)
        winner (cond
                 (and (pos? (get-in fallen [:vikings :tower] 0))
                      (pos? (get-in fallen [:barbarians :tower] 0))) :draw
                 (pos? (get-in fallen [(other turn) :tower] 0)) turn
                 (pos? (get-in fallen [turn :tower] 0)) (other turn)
                 (every? zero? (vals discs)) :draw)
        nxt (if (pos? (long (discs (other turn)))) (other turn) turn)]
    (assoc game
           :fallen fallen
           :winner winner
           :phase (if winner :over :aiming)
           :turn (if winner turn nxt))))

(defn advance
  "One step of the world, and the end of the turn once it has come to
  rest -- still for a third of a second, or eight seconds after the shot."
  ([game] (advance game (/ 1.0 60.0)))
  ([{:keys [phase] :as game} dt]
   (let [game (update game :world solver/step dt)]
     (if (not= phase :flying)
       game
       (let [clock (+ (double (:clock game 0.0)) dt)
             calm (if (calm? (:world game)) (inc (long (:calm game 0))) 0)
             game (assoc game :clock clock :calm calm)]
         (if (or (>= calm 20) (>= clock 8.0))
           (resolve-turn game)
           game))))))

;; ---------------------------------------------------------------------------
;; The computer's turn

(defn- enemy-tower [world side]
  (some #(when (and (= (:kind %) :tower) (= (:side %) (other side))) %) (:bodies world)))

(defn aim-at
  "The shot that would land on `target` from `side`'s `weapon`, in a
  vacuum: `{:weapon :yaw :pitch :speed}`.

  The catapult throws at `pitch` degrees and solves for the speed, which
  with the pitch fixed is one square root:
  `v^2 = g d^2 / (2 cos^2 p (d tan p - h))`. The crossbow aims its yaw
  and shoots hard enough to arrive still moving."
  [side weapon target & [pitch]]
  (let [p0 (launch-point side weapon)
        {:keys [right back]} (frame side)
        fwd (v/scale back -1.0)
        dx (v/sub target p0)
        along (v/dot dx fwd)
        across (v/dot dx right)
        yaw (math/to-degrees (math/atan2 across along))
        d (math/hypot along across)
        g 10.0]
    (if (= weapon :crossbow)
      {:weapon :crossbow :yaw yaw :pitch 0.0
       :speed (min 22.0 (+ 4.0 (math/sqrt (* 2.0 0.6 g d))))}
      (let [pitch (double (or pitch 35.0))
            p (math/to-radians pitch)
            h (nth dx 1)
            c (math/cos p)
            v2 (/ (* g d d) (* 2.0 c c (- (* d (math/tan p)) h)))]
        {:weapon :catapult :yaw yaw :pitch pitch :speed (math/sqrt v2)}))))

(defn tower-face
  "The point to aim a catapult at: high on the face of the enemy tower
  that looks at `side`. A disc that comes down on the top of the tower
  sits there; one that hits the face near the top is the one that tips
  it."
  [world side]
  (let [tower (enemy-tower world side)
        {:keys [back]} (frame side)]
    (v/add (:pos tower) (v/add [0.0 0.8 0.0] (v/scale back 0.6)))))

(defn computer-shot
  "The computer's choice for `side`: mostly the catapult at the face of
  the enemy tower, sometimes the crossbow through the gate at its foot, with
  `error` -- a spread in degrees and in fractions of the speed -- drawn
  from `rand`, a function of no arguments returning a number in [0, 1)."
  [{:keys [world]} side error rand]
  (let [tower (enemy-tower world side)
        [x _ z] (:pos tower)
        weapon (if (< (rand) 0.75) :catapult :crossbow)
        target (if (= weapon :catapult) (tower-face world side) [x 0.1 z])
        shot (aim-at side weapon target)
        jitter (fn [] (* 2.0 (- (rand) 0.5)))]
    (-> shot
        (update :yaw + (* (double error) 5.0 (jitter)))
        (update :speed * (+ 1.0 (* (double error) 0.08 (jitter)))))))
