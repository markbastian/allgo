(ns allgo.game.crossbows
  "Crossbows and Catapults, played in the browser.

  The rules and the physics are `allgo.simulation.crossbows`; this is the
  table they are played on. It has a page of its own, `cnc/`, which
  `allgo.app` starts it on.

  Aiming: the arrow keys turn the shot and raise or lower a catapult's
  throw, W and S set the power, Tab changes weapon, and Space shoots.
  The path drawn is only the start of the throw, as far as someone
  sighting a real catapult could judge -- where it lands is found out by
  shooting."
  (:require [allgo.simulation.crossbows :as c]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

;; ---------------------------------------------------------------------------
;; Colors, after the set

(def ^:private palette
  {:vikings {:stone 0xcdb892 :stone-dark 0xa8936c :cloth 0xc8283c :disc 0xd0304a
             :warrior 0x6b4a2e :engine 0x3a2a1e :name "Vikings"}
   :barbarians {:stone 0x8190a8 :stone-dark 0x66738a :cloth 0x2458c8 :disc 0x2a62cc
                :warrior 0xe6e2da :engine 0xc0692e :name "Barbarians"}})

(defn- material [color & [opts]]
  (THREE/MeshStandardMaterial. (clj->js (merge {:color color :roughness 0.75 :metalness 0.05} opts))))

(defn- mesh [^js geometry ^js mat & [x y z]]
  (doto (THREE/Mesh. geometry mat)
    (-> .-position (.set (or x 0) (or y 0) (or z 0)))))

(defn- group [& children]
  (let [g (THREE/Group.)]
    (doseq [ch children] (.add g ch))
    g))

;; ---------------------------------------------------------------------------
;; Pieces. Each is drawn about its body's center, in the body's own frame.

(defn- tower-mesh [side]
  (let [{:keys [stone stone-dark]} (palette side)
        [w h _] (:size (first (filter #(= :tower (:kind %)) (c/board))))
        r (* 0.5 w)
        body (mesh (THREE/CylinderGeometry. r (* 1.04 r) h 24) (material stone))
        top (* 0.5 h)
        merlons (for [i (range 8)
                      :let [a (* i (/ Math/PI 4))]]
                  (doto (mesh (THREE/BoxGeometry. 0.28 0.3 0.22) (material stone)
                              (* 0.88 r (Math/cos a)) (+ top 0.15) (* 0.88 r (Math/sin a)))
                    (-> .-rotation (.set 0 (- a) 0))))
        ;; The door looks at the enemy, which is -z in the body's frame.
        door (mesh (THREE/BoxGeometry. 0.42 0.7 0.1) (material 0x2a1d14) 0 (+ (- top) 0.35) (- (* 0.98 r)))
        band (mesh (THREE/CylinderGeometry. (* 1.02 r) (* 1.02 r) 0.12 24) (material stone-dark) 0 (- top 0.35) 0)]
    (apply group body door band merlons)))

(defn- wall-mesh [side]
  (let [{:keys [stone stone-dark]} (palette side)
        body (mesh (THREE/BoxGeometry. 1.8 1.0 0.5) (material stone))
        merlons (for [x [-0.6 0.0 0.6]]
                  (mesh (THREE/BoxGeometry. 0.34 0.22 0.5) (material stone-dark) x 0.61 0))]
    (apply group body merlons)))

(defn- flag-mesh [side]
  (let [{:keys [cloth]} (palette side)
        pole (mesh (THREE/CylinderGeometry. 0.04 0.04 1.4 8) (material 0x3b3b3b))
        base (mesh (THREE/CylinderGeometry. 0.24 0.28 0.08 16) (material cloth) 0 -0.66 0)
        shape (doto (THREE/Shape.) (.moveTo 0 0) (.lineTo 0.6 -0.15) (.lineTo 0 -0.34) (.lineTo 0 0))
        pennant (doto (mesh (THREE/ShapeGeometry. shape) (material cloth {:side THREE/DoubleSide}) 0.03 0.68 0))]
    (group pole base pennant)))

(defn- warrior-mesh [side]
  (let [{:keys [warrior]} (palette side)
        body (mesh (THREE/CylinderGeometry. 0.13 0.17 0.5 12) (material warrior) 0 -0.1 0)
        head (mesh (THREE/SphereGeometry. 0.12 12 10) (material warrior) 0 0.24 0)
        horns (when (= side :vikings)
                (for [s [-1 1]]
                  (doto (mesh (THREE/ConeGeometry. 0.035 0.16 8) (material 0xeee6d2) (* s 0.12) 0.32 0)
                    (-> .-rotation (.set 0 0 (* s -0.7))))))
        shield (mesh (THREE/CylinderGeometry. 0.13 0.13 0.03 14) (material (:cloth (palette side))) 0 -0.05 -0.17)]
    (.set (.-rotation shield) (/ Math/PI 2) 0 0)
    (apply group body head shield horns)))

(defn- disc-mesh [side]
  (let [{:keys [disc]} (palette side)
        body (mesh (THREE/CylinderGeometry. 0.28 0.28 0.14 24) (material disc {:roughness 0.4}))
        face (mesh (THREE/CylinderGeometry. 0.18 0.18 0.02 24) (material 0xf2e2b8) 0 0.075 0)]
    (group body face)))

(defn- body-mesh [b]
  (case (:kind b)
    :tower (tower-mesh (:side b))
    :wall (wall-mesh (:side b))
    :flag (flag-mesh (:side b))
    :warrior (warrior-mesh (:side b))
    :disc (disc-mesh (:side b))
    ;; The ground is drawn as the board, below; its body needs nothing.
    (THREE/Group.)))

(defn- beam
  "A square timber `thick` across, from point `a` to point `b`."
  [^js mat thick [ax ay az] [bx by bz]]
  (let [a (THREE/Vector3. ax ay az)
        b (THREE/Vector3. bx by bz)
        ^js m (mesh (THREE/BoxGeometry. thick thick (.distanceTo a b)) mat)]
    (.copy (.-position m) (.multiplyScalar (.add (.clone a) b) 0.5))
    (.lookAt m b)
    m))

(defn- catapult-mesh [side]
  ;; Cocked: the arm runs from the axle up to the cup, and the cup is
  ;; where `c/launchers` says the disc leaves from.
  (let [{:keys [engine]} (palette side)
        m (material engine)
        [_ cup-y _] (:at (c/launchers :catapult))
        rails (for [x [-0.45 0.45]] (beam m 0.14 [x 0.07 -0.8] [x 0.07 0.8]))
        uprights (for [x [-0.45 0.45]] (beam m 0.12 [x 0.07 0.35] [x 0.95 0.35]))
        braces (for [x [-0.45 0.45]] (beam m 0.08 [x 0.07 -0.5] [x 0.9 0.3]))
        axle (beam m 0.1 [-0.5 0.95 0.35] [0.5 0.95 0.35])
        arm (beam m 0.12 [0 0.95 0.35] [0 (- cup-y 0.06) 0])
        cup (mesh (THREE/CylinderGeometry. 0.3 0.24 0.12 16) (material (:disc (palette side))) 0 cup-y 0)]
    (apply group axle arm cup (concat rails uprights braces))))

(defn- crossbow-mesh [side]
  (let [{:keys [engine]} (palette side)
        m (material engine)
        stock (mesh (THREE/BoxGeometry. 0.22 0.16 1.4) m 0 0.12 0.2)
        bow (doto (mesh (THREE/TorusGeometry. 0.7 0.05 8 24 Math/PI) m 0 0.14 -0.25)
              (-> .-rotation (.set (/ Math/PI 2) 0 0)))
        legs (for [x [-0.35 0.35]] (mesh (THREE/BoxGeometry. 0.1 0.1 1.0) m x 0.05 0.3))]
    (apply group stock bow legs)))

(defn- board-mesh []
  (let [grass (mesh (THREE/BoxGeometry. 16 1 26) (material 0x55803e) 0 -0.5 0)
        mats (for [side c/sides
                   :let [z (if (= side :vikings) 7.0 -7.0)]]
               (mesh (THREE/BoxGeometry. 12 0.02 10) (material (:stone-dark (palette side)) {:roughness 0.95}) 0 0.011 z))
        river (mesh (THREE/BoxGeometry. 16 0.02 1.2) (material 0x3a78b0 {:roughness 0.3}) 0 0.012 0)]
    (apply group grass river mats)))

;; ---------------------------------------------------------------------------
;; The game page

(defn- aim-for
  "A fresh aim for `side` and `weapon`: pointed at the enemy tower, at the
  weapon's usual pitch and power. The direction is given away; the
  distance is not."
  [game side weapon]
  (let [{:keys [yaw]} (c/aim-at side weapon (c/tower-face (:world game) side))
        {:keys [default-pitch default-speed]} (c/launchers weapon)]
    {:weapon weapon :yaw (Math/round yaw) :pitch default-pitch :speed default-speed}))

(defn- el [tag cls & [text]]
  (doto (js/document.createElement tag)
    (-> .-className (set! cls))
    (-> .-textContent (set! (or text "")))))

(defn init! [^js container]
  (let [scene (THREE/Scene.)
        camera (THREE/PerspectiveCamera. 50 1 0.1 400)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        ^js orbit (OrbitControls. camera (.-domElement renderer))
        ;; HUD
        hud (el "div" "hud")
        banner (el "div" "hud__banner")
        tally (el "div" "hud__tally")
        weapon-line (el "div" "hud__weapon")
        help (el "div" "hud__help" "← → aim · ↑ ↓ raise (catapult) · W S power · Tab weapon · Space shoot · drag to look")
        buttons (el "div" "hud__buttons")
        new-btn (el "button" "hud__button" "New game")
        mode-btn (el "button" "hud__button" "")
        overlay (el "div" "hud__overlay")
        ;; State
        state (atom nil)
        meshes (atom [])
        running? (atom true)
        ;; A tube, not a line: WebGL draws every line one pixel wide.
        path (THREE/Mesh. (THREE/BufferGeometry.)
                          (THREE/MeshBasicMaterial. #js {:color 0xfff27a :transparent true :opacity 0.85}))]
    (.setPixelRatio renderer (min 2 (or js/window.devicePixelRatio 1)))
    (.appendChild container (.-domElement renderer))
    (doto buttons (.appendChild new-btn) (.appendChild mode-btn))
    (doto hud (.appendChild banner) (.appendChild tally) (.appendChild weapon-line)
          (.appendChild buttons) (.appendChild help) (.appendChild overlay))
    (.appendChild container hud)
    (set! (.-background scene) (THREE/Color. 0x9cc4e4))
    (set! (.-fog scene) (THREE/Fog. 0x9cc4e4 40 90))
    (.add scene (THREE/HemisphereLight. 0xdfeaf5 0x4a5530 0.55))
    (let [sun (THREE/DirectionalLight. 0xfff3dc 1.1)]
      (.set (.-position sun) 10 18 6)
      (.add scene sun))
    (.add scene (board-mesh))
    (doseq [side c/sides
            [weapon make] [[:catapult catapult-mesh] [:crossbow crossbow-mesh]]]
      (let [{[x y z] :pos [qx qy qz qw] :rot} (c/engine-pose side weapon)
            ^js m (make side)]
        (.set (.-position m) x y z)
        (.set (.-quaternion m) qx qy qz qw)
        (.add scene m)))
    (.add scene path)
    (set! (.-enableDamping orbit) true)
    (letfn [(human? [side] (or (= :two (:mode @state)) (= side :vikings)))
            (sync-meshes! []
              (let [bodies (:bodies (:world (:game @state)))]
                ;; A disc fired this turn is a body with no mesh yet.
                (while (< (count @meshes) (count bodies))
                  (let [^js m (body-mesh (nth bodies (count @meshes)))]
                    (.add scene m)
                    (swap! meshes conj m)))
                (doseq [[^js m b] (map vector @meshes bodies)]
                  (let [[x y z] (:pos b) [qx qy qz qw] (:rot b)]
                    (.set (.-position m) x y z)
                    (.set (.-quaternion m) qx qy qz qw)))))
            (look-from! [side]
              ;; Behind the shooter's catapult and out to its side, so
              ;; the shooter's own tower is not in the way of the other.
              (let [[x y z] (c/launch-point side :catapult)
                    [tx ty tz] (:pos (some #(when (and (= :tower (:kind %)) (not= side (:side %))) %)
                                           (:bodies (:world (:game @state)))))
                    ;; A side's right and back, along x and z.
                    s (if (= side :vikings) 1 -1)]
                (swap! state assoc :camera-goal {:pos [(+ x (* s 2.0)) (+ y 5.0) (+ z (* s 6.5))]
                                                 :target [(* 0.5 tx) (* 0.4 ty) (* 0.45 tz)]})))
            (draw-aim! []
              (let [{:keys [game aims]} @state
                    {:keys [turn phase]} game
                    {:keys [weapon yaw pitch speed]} (aims turn)
                    pts (vec (c/flight-path turn weapon yaw pitch speed))
                    ;; Only the start of the throw: a third of the arc, or
                    ;; the first few meters of a crossbow's slide.
                    shown (take (max 2 (quot (count pts) 3)) pts)
                    ^js old (.-geometry path)]
                (set! (.-visible path) (and (= phase :aiming) (human? turn)))
                (set! (.-geometry path)
                      (THREE/TubeGeometry.
                       (THREE/CatmullRomCurve3. (clj->js (map (fn [[x y z]] (THREE/Vector3. x y z)) shown)))
                       48 0.05 6 false))
                (.dispose old)))
            (draw-hud! []
              (let [{:keys [game aims mode]} @state
                    {:keys [turn phase winner discs]} game
                    name-of #(:name (palette %))
                    {:keys [weapon pitch speed]} (aims turn)]
                (set! (.-textContent banner)
                      (cond
                        (= winner :draw) "Out of discs -- a draw"
                        winner (str (name-of winner) " win!")
                        (= phase :flying) "…"
                        (human? turn) (str (name-of turn) " to shoot")
                        :else (str (name-of turn) " are aiming…")))
                (set! (.-className banner) (str "hud__banner hud__banner--" (name (or (#{:vikings :barbarians} winner) turn))))
                (set! (.-textContent tally)
                      (str "Vikings " (apply str (repeat (discs :vikings) "●")) (apply str (repeat (- 10 (discs :vikings)) "○"))
                           "   ·   Barbarians " (apply str (repeat (discs :barbarians) "●")) (apply str (repeat (- 10 (discs :barbarians)) "○"))))
                (set! (.-textContent weapon-line)
                      (if (human? turn)
                        (str (if (= weapon :catapult) "Catapult" "Crossbow")
                             (when (= weapon :catapult) (str " · " (Math/round pitch) "° up"))
                             " · power " (.toFixed speed 1))
                        ""))
                (set! (.-textContent mode-btn) (if (= mode :two) "Two players" "Vs computer"))
                (set! (.-textContent overlay)
                      (cond (= winner :draw) "Both towers stand. Press N for a new game."
                            winner (str "The " (name-of (c/other winner)) "' tower has fallen. Press N for a new game.")
                            :else ""))
                (set! (.-hidden overlay) (nil? winner))))
            (new-game! []
              (doseq [^js m @meshes] (.remove scene m))
              (reset! meshes [])
              (let [game (c/new-game)]
                (swap! state assoc
                       :game game
                       :aims {:vikings (aim-for game :vikings :catapult)
                              :barbarians (aim-for game :barbarians :catapult)}
                       :think nil
                       :shown-turn nil))
              (sync-meshes!)
              (look-from! :vikings)
              (draw-aim!)
              (draw-hud!))
            (shoot! [side]
              (let [{:keys [weapon yaw pitch speed]} (get-in @state [:aims side])]
                (swap! state update :game c/fire weapon yaw pitch speed)
                (draw-aim!)
                (draw-hud!)))
            (adjust! [f]
              (let [{:keys [game]} @state
                    side (:turn game)]
                (when (and (= :aiming (:phase game)) (human? side) (nil? (:winner game)))
                  (swap! state update-in [:aims side] f)
                  (draw-aim!)
                  (draw-hud!))))
            (clamp-speed [{:keys [weapon] :as aim}]
              (let [[lo hi] (:speed (c/launchers weapon))]
                (update aim :speed #(max lo (min hi %)))))
            (on-key [^js e]
              (let [step (if (.-shiftKey e) 5.0 1.0)
                    handled
                    (case (.-key e)
                      "ArrowLeft" (adjust! #(update % :yaw - step))
                      "ArrowRight" (adjust! #(update % :yaw + step))
                      "ArrowUp" (adjust! #(update % :pitch (fn [p] (min 75.0 (+ p step)))))
                      "ArrowDown" (adjust! #(update % :pitch (fn [p] (max 10.0 (- p step)))))
                      ("w" "W") (adjust! #(clamp-speed (update % :speed + (* 0.25 step))))
                      ("s" "S") (adjust! #(clamp-speed (update % :speed - (* 0.25 step))))
                      "Tab" (adjust! (fn [{:keys [weapon]}]
                                       (aim-for (:game @state) (:turn (:game @state))
                                                (if (= weapon :catapult) :crossbow :catapult))))
                      " " (let [{:keys [game]} @state]
                            (when (and (= :aiming (:phase game)) (human? (:turn game)) (nil? (:winner game)))
                              (shoot! (:turn game))))
                      ("n" "N") (new-game!)
                      :unhandled)]
                (when-not (= handled :unhandled) (.preventDefault e))))
            (computer! []
              ;; A pause to look like thinking, then the shot, drawn as
              ;; the aim it was.
              (let [{:keys [game think]} @state
                    now (js/performance.now)]
                (cond
                  (nil? think) (swap! state assoc :think now)
                  (> (- now think) 900)
                  (let [side (:turn game)
                        shot (c/computer-shot game side 1.0 js/Math.random)]
                    (swap! state #(-> % (assoc-in [:aims side] shot) (assoc :think nil)))
                    (shoot! side)))))
            (on-resize []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (set! (.-aspect camera) (/ w h))
                  (.updateProjectionMatrix camera)
                  (.setSize renderer w h))))
            (animate []
              (when @running?
                (js/requestAnimationFrame animate)
                (let [before (:game @state)]
                  (swap! state update :game c/advance)
                  (let [{:keys [game shown-turn]} @state]
                    ;; A turn has just been settled: redraw, and walk the
                    ;; camera round to whoever shoots next.
                    (when (or (not= (:phase before) (:phase game)) (not= shown-turn (:turn game)))
                      (when (and (= :aiming (:phase game)) (not= shown-turn (:turn game)))
                        (swap! state assoc :shown-turn (:turn game))
                        (look-from! (:turn game)))
                      (draw-aim!)
                      (draw-hud!))
                    (when (and (= :aiming (:phase game)) (nil? (:winner game)) (not (human? (:turn game))))
                      (computer!))))
                (sync-meshes!)
                (when-let [{[px py pz] :pos [tx ty tz] :target} (:camera-goal @state)]
                  ;; Ease toward the goal rather than jump; let the reader
                  ;; take over with the mouse once it has arrived.
                  (let [^js p (.-position camera) ^js t (.-target orbit)]
                    (.lerp p (THREE/Vector3. px py pz) 0.06)
                    (.lerp t (THREE/Vector3. tx ty tz) 0.06)
                    (when (< (.distanceTo p (THREE/Vector3. px py pz)) 0.05)
                      (swap! state dissoc :camera-goal))))
                (.update orbit)
                (.render renderer scene camera)))]
      (reset! state {:mode :one})
      (.addEventListener new-btn "click" (fn [_] (new-game!)))
      (.addEventListener mode-btn "click" (fn [_]
                                            (swap! state update :mode #(if (= % :two) :one :two))
                                            (new-game!)))
      (js/window.addEventListener "keydown" on-key)
      (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
      (.set (.-position camera) 0 16 26)
      (new-game!)
      (on-resize)
      (animate))))

(defonce ^:private started (atom false))

(defn start!
  "Starts the game in the page's `#crossbows` element, once."
  []
  (when-let [container (and (not @started) (js/document.getElementById "crossbows"))]
    (reset! started true)
    (init! container)))
