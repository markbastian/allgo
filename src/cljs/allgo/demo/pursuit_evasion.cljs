(ns allgo.demo.pursuit-evasion
  "Canvas demo for `allgo.robotics.pursuit`: hide-and-seek played so the
  seeker can prove no one is left hiding, after LaValle, *Planning
  Algorithms*, Cambridge University Press, 2006, section 12.4.

  The pursuer (white) sees what its lines of sight reach (gold); the rest
  is shadow, coral where an evader might still hide and blue where it
  certainly cannot. Plan searches the information space -- a vertex for
  each place and each labeling of the shadow's components -- breadth
  first, then walks the plan it finds, naming the critical events as they
  happen: a component appearing (cleared), disappearing, splitting, or
  merging, the one that can spread contamination back. Drive lets you
  drag the pursuer and try it yourself. With three hooks there is no
  plan: Example 12.7's three hiding places, no two of which can be
  watched from common ground (their V(V(p)) drawn), so the evader can
  always slip between them."
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.robotics.pursuit :as pursuit]
            [allgo.search :as search]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:environment "gallery" :mode "plan" :speed 3})

(defn- gallery
  "A hall `w` long with rooms opening off it, each through a door
  narrower than the room, north at the `north` x's and south at the
  `south`."
  [north south w]
  (vec (concat [[0 0]]
               (mapcat (fn [cx] [[(- cx 1) 0] [(- cx 1) -2] [(- cx 3) -2] [(- cx 3) -5]
                                 [(+ cx 3) -5] [(+ cx 3) -2] [(+ cx 1) -2] [(+ cx 1) 0]])
                       south)
               [[w 0] [w 2]]
               (mapcat (fn [cx] [[(+ cx 1) 2] [(+ cx 1) 4] [(+ cx 3) 4] [(+ cx 3) 7]
                                 [(- cx 3) 7] [(- cx 3) 4] [(- cx 1) 4] [(- cx 1) 2]])
                       (reverse north))
               [[0 2]])))

(def ^:private scenes
  {"gallery" {:polygon (gallery [4 12 20] [8 16] 24) :start [1 1]}
   "two hooks" {:polygon [[-3 -1] [3 -1] [3 -7] [11 -7] [11 -5] [5 -5] [5 1] [-5 1] [-5 -5]
                          [-11 -5] [-11 -7] [-3 -7]]
                :start [0 0]}
   "three hooks" {:polygon [[-3 -1] [3 -1] [3 -7] [11 -7] [11 -5] [5 -5] [5 1] [1 1] [1 3] [7 3]
                            [7 11] [5 11] [5 5] [-1 5] [-1 1] [-5 1] [-5 -5] [-11 -5] [-11 -7] [-3 -7]]
                  :start [0 0]
                  :hides [[10 -6] [-10 -6] [6 10]]}})

(def ^:private hide-colors ["rgba(125,255,160,0.45)" "rgba(127,209,255,0.45)" "rgba(214,160,255,0.45)"])

(defn- event-text [{:keys [contaminated]} [kind & args]]
  (case kind
    :appear "a component appears -- cleared, it was just seen"
    :disappear "a component comes into view and disappears"
    :split "a component splits, each part keeping its label"
    :merge (if (contains? contaminated (last args))
             "two components merge -- contaminated: the work is undone where they met"
             "two cleared components merge, still cleared")))

(defn- watched
  "V(V(p)): the probes seen from somewhere probe `k` sees."
  [{:keys [probes] :as env} k]
  (into #{} (mapcat #(pursuit/visible env (probes %))) (pursuit/visible env (probes k))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        model (atom nil)
        dragging? (atom false)]
    (set! (.-style canvas) "display:block;touch-action:none")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(scene [] (scenes (.-environment controls)))
            (restart!
              ([] (restart! (:start (scene))))
              ([p]
               (let [{:keys [polygon hides]} (scene)
                     env (pursuit/environment polygon 0.5)
                     state (pursuit/start env p)]
                 (reset! model {:env env :state state :trail [p] :event nil
                                :phase (if (= "plan" (.-mode controls)) :searching :driving)
                                :search (pursuit/search-seq env p) :expanded 0 :reached #{}
                                :hides (mapv #(watched env (pursuit/nearest-probe env %)) hides)}))))
            (search! []
              ;; as many search steps as fit in a few milliseconds
              (let [deadline (+ (js/performance.now) 10)]
                (loop [{:keys [search] :as m} @model]
                  (let [[s & more] search]
                    (cond
                      (nil? s) (reset! model (assoc m :phase :no-plan))
                      (:solved? s)
                      (let [path (mapv (comp (:probes (:env m)) first) (search/path s))]
                        (reset! model (assoc m :phase :walking :search nil :path path :t 0.0
                                             :states (vec (pursuit/clearing-seq (:env m) path)))))
                      (> (js/performance.now) deadline) (reset! model m)
                      :else (recur (cond-> (-> m (assoc :search more) (update :expanded inc))
                                     (:current s) (update :reached conj (first (:current s))))))))))
            (walk! []
              (let [{:keys [t path states trail]} @model
                    t' (min (dec (count path)) (+ t (* 0.04 (.-speed controls))))
                    i (long (math/floor t')) j (min (inc i) (dec (count path)))
                    at (polygon/lerp (path i) (path j) (- t' i))
                    events (when (> i (long (math/floor t))) (seq (pursuit/events (:shadow (states (dec i))) (:shadow (states i)))))]
                (swap! model (fn [m] (cond-> (assoc m :t t' :state (assoc (states i) :at at) :trail (conj trail at))
                                       events (assoc :event (event-text (states i) (first events)))
                                       (= t' (dec (count path))) (assoc :phase :done))))))
            (drive! []
              ;; toward the pointer, a third of a probe spacing at a time
              (let [{:keys [env state target trail]} @model]
                (when target
                  (loop [s state k 0 trail trail event nil]
                    (let [d (polygon/distance (:at s) target)]
                      (if (or (< d 1e-6) (= k 3))
                        (swap! model (fn [m] (cond-> (assoc m :state s :trail trail) event (assoc :event event))))
                        (let [p (polygon/lerp (:at s) target (min 1.0 (/ (/ (:spacing env) 3) d)))
                              s' (pursuit/clearing-step env s p)
                              e (first (pursuit/events (:shadow s) (:shadow s')))]
                          (recur s' (inc k) (conj trail p) (if e (event-text s' e) event)))))))))
            (layout []
              (let [{:keys [polygon]} (:env @model)
                    xs (map first polygon) ys (map second polygon)
                    x0 (reduce min xs) x1 (reduce max xs) y0 (reduce min ys) y1 (reduce max ys)
                    cw (.-width canvas) ch (.-height canvas)
                    ;; clear of the title above, the caption below and the controls
                    k (min (/ (- cw 300) (- x1 x0)) (/ (- ch 150) (- y1 y0)))
                    ox (- (* 0.5 (- cw 260)) (* k 0.5 (+ x0 x1)))
                    oy (+ 60 (* 0.5 (- ch 150)) (* k 0.5 (+ y0 y1)))]
                {:k k :px (fn [[x y]] [(+ ox (* k x)) (- oy (* k y))])
                 :uv (fn [sx sy] [(/ (- sx ox) k) (/ (- oy sy) k)])}))
            (path! [pts]
              (.beginPath ctx)
              (doseq [[i [x y]] (map-indexed vector pts)] (if (zero? i) (.moveTo ctx x y) (.lineTo ctx x y))))
            (draw! []
              (let [{:keys [env state trail phase expanded reached hides event]} @model
                    {:keys [probes spacing polygon]} env
                    {:keys [k px]} (layout)
                    cw (.-width canvas) ch (.-height canvas)
                    sq (* k spacing)
                    dirty (pursuit/contaminated-probes state)
                    hidden (set (mapcat identity (:shadow state)))]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (path! (map px polygon)) (.closePath ctx)
                (set! (.-fillStyle ctx) "#161b25") (.fill ctx)
                (path! (map px (pursuit/visibility-polygon env (:at state)))) (.closePath ctx)
                (set! (.-fillStyle ctx) "rgba(255,212,121,0.16)") (.fill ctx)
                (doseq [i hidden :let [[x y] (px (probes i))]]
                  (set! (.-fillStyle ctx) (if (dirty i) "rgba(255,143,112,0.55)" "rgba(127,209,255,0.32)"))
                  (.fillRect ctx (- x (* 0.45 sq)) (- y (* 0.45 sq)) (* 0.9 sq) (* 0.9 sq)))
                (when (= phase :no-plan)
                  (doseq [[c ks] (map vector hide-colors hides) i ks :let [[x y] (px (probes i))]]
                    (set! (.-fillStyle ctx) c) (.fillRect ctx (- x (* 0.5 sq)) (- y (* 0.5 sq)) sq sq)))
                (when (= phase :searching)
                  (set! (.-fillStyle ctx) "rgba(125,255,160,0.8)")
                  (doseq [i reached :let [[x y] (px (probes i))]]
                    (.fillRect ctx (- x 1.5) (- y 1.5) 3 3)))
                (path! (map px polygon)) (.closePath ctx)
                (set! (.-strokeStyle ctx) "#9fb4d8") (set! (.-lineWidth ctx) 2) (.stroke ctx)
                (path! (map px trail))
                (set! (.-strokeStyle ctx) "rgba(255,212,121,0.7)") (set! (.-lineWidth ctx) 1.5) (.stroke ctx)
                (let [[x y] (px (:at state))]
                  (set! (.-fillStyle ctx) "#ffffff") (.beginPath ctx) (.arc ctx x y 6 0 (* 2 math/PI)) (.fill ctx))
                (set! (.-textContent caption)
                      (str (case phase
                             :searching (str "searching the information space: " expanded " vertices expanded, "
                                             (count reached) " places reached")
                             :walking (str "walking the plan: " (count (:contaminated state)) " of "
                                           (count (:shadow state)) " shadow components contaminated")
                             :done (str "cleared: no evader can be anywhere -- found by searching " expanded " vertices")
                             :no-plan (str "no plan from here: all " expanded " reachable vertices searched. "
                                           "Three hiding places whose V(V(p)) (green, blue, violet) share no ground"
                                           " -- switch to drive and try")
                             :driving (str (count (:contaminated state)) " of " (count (:shadow state))
                                           " shadow components contaminated"
                                           (when (pursuit/cleared? state) " -- cleared!")))
                           (when event (str "\n" event))
                           "\ngold: seen  coral: might hide the evader  blue: shadow known clear"
                           (if (= "plan" (.-mode controls)) "  click: start here" "  drag: move the pursuer")))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (case (:phase @model)
                  :searching (search!)
                  :walking (walk!)
                  :driving (drive!)
                  nil)
                (draw!)))
            (pointer [^js e]
              (let [r (.getBoundingClientRect canvas)]
                ((:uv (layout)) (- (.-clientX e) (.-left r)) (- (.-clientY e) (.-top r)))))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))]
      (.addEventListener canvas "pointerdown"
                         (fn [^js e]
                           (let [p (pointer e) {:keys [env state]} @model]
                             (when (polygon/contains-point? (:polygon env) p)
                               (if (= "plan" (.-mode controls))
                                 (restart! p)
                                 (do (reset! dragging? true)
                                     (when (pursuit/sees? env (:at state) p) (swap! model assoc :target p))))))))
      (.addEventListener canvas "pointermove"
                         (fn [^js e]
                           (when @dragging?
                             (let [p (pointer e) {:keys [env state]} @model]
                               ;; only where the pursuer can go straight, without passing through a wall
                               (when (and (polygon/contains-point? (:polygon env) p) (pursuit/sees? env (:at state) p))
                                 (swap! model assoc :target p))))))
      (.addEventListener js/window "pointerup" (fn [_] (reset! dragging? false)))
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "environment" (clj->js (keys scenes))) #(restart!))
        (.onChange (.add gui controls "mode" #js ["plan" "drive"]) #(restart!))
        (.add gui controls "speed" 1 10 1)
        (.add gui #js {:restart #(restart!)} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "pursuit-evasion")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
