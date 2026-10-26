(ns allgo.demo.city-plan
  "Canvas demo for `allgo.procedural.city.roads` and `.lots`: Parish and
  Mueller's city, grown as a plan.

  The roads grow a proposal at a time, in the order their L-system's
  delays put them: the highways first, reaching for the centers of
  population (the red of the background), then the streets filling in
  behind them in each district's pattern. Watch a street meet another --
  it stops at the crossing, or is pulled to a junction just ahead of it,
  or stretched to a road just short of it -- and the tree becomes a net.
  When the growth stops, the blocks are cut out of the net and divided
  into lots, colored by zone: downtown, commercial, residential. The gray
  left in the middle of a block is land with no street to face, which the
  paper throws away.

  `pattern` is the street-pattern map: `mixed` gives a radial old town
  (Paris) around the main center, grids (New York) in the busy districts
  around it, and the paper's basic rule in the thin outskirts; the others
  apply one rule everywhere, `terrain` being San Francisco's, which
  follows the hills. Scroll to zoom, drag to pan."
  (:require [allgo.procedural.city.land :as land]
            [allgo.procedural.city.lots :as lots]
            [allgo.procedural.city.roads :as roads]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls
  #js {:seed     2
       :pattern  "mixed"
       :speed    6
       :show     "lots"
       :backdrop "population"})

(def ^:private zone-colors
  {:downtown "#d9675e" :commercial "#e3a94f" :residential "#7fae6e"})

(defn- backdrop
  "The land, painted once to an offscreen canvas `n` pixels square: water,
  and the population or the elevation."
  [l n what]
  (let [c (js/document.createElement "canvas")
        _ (set! (.-width c) n)
        _ (set! (.-height c) n)
        ctx (.getContext c "2d")
        img (.createImageData ctx n n)
        data (.-data img)
        size (:size l)]
    (dotimes [j n]
      (dotimes [i n]
        (let [x (- (* (/ (+ i 0.5) n) 2 size) size)
              y (- size (* (/ (+ j 0.5) n) 2 size))
              o (* 4 (+ i (* n j)))
              e (min 1.0 (/ ((:elevation l) x y) 70.0))
              base (+ 0.30 (* 0.28 e))
              [r g b] (cond
                        ((:water? l) x y) [0.22 0.38 0.56]
                        (= what "population") (let [p ((:population l) x y)]
                                                [(+ base (* 0.35 p)) (+ base (* 0.02 p)) (- (* 0.85 base) (* 0.05 p))])
                        :else [(* 0.9 base) (+ base 0.05) (* 0.75 base)])]
          (aset data o (* 255 r)) (aset data (+ o 1) (* 255 g))
          (aset data (+ o 2) (* 255 b)) (aset data (+ o 3) 255))))
    (.putImageData ctx img 0 0)
    c))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        last-frame (atom nil)
        state (atom nil)
        ;; The view: meters per pixel and the point at the middle.
        view (atom {:scale nil :cx 0.0 :cy 0.0})]
    (set! (.-style canvas) "display:block;cursor:grab;touch-action:none")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;padding:3px 8px;border-radius:4px;background:rgba(10,14,22,0.7);"
               "font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;pointer-events:none"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(regenerate! []
              (let [l (land/land {:seed (.-seed controls) :pattern (keyword (.-pattern controls))})]
                (reset! state {:land l
                               :back (backdrop l 520 (.-backdrop controls))
                               :roads (roads/start l {:seed (.-seed controls)})
                               :blocks nil :lots nil})
                (swap! view assoc :scale nil)))
            (fit! []
              (let [{:keys [land]} @state
                    w (.-width canvas) h (.-height canvas)]
                (when (and land (nil? (:scale @view)))
                  (reset! view {:scale (/ (* 2.1 (:size land)) (min w h)) :cx 0.0 :cy 0.0}))))
            (to-screen [[x y]]
              (let [{:keys [scale cx cy]} @view]
                [(+ (* 0.5 (.-width canvas)) (/ (- x cx) scale))
                 (- (* 0.5 (.-height canvas)) (/ (- y cy) scale))]))
            (poly-path! [pts]
              (.beginPath ctx)
              (doseq [[k p] (map-indexed vector pts)]
                (let [[sx sy] (to-screen p)]
                  (if (zero? k) (.moveTo ctx sx sy) (.lineTo ctx sx sy))))
              (.closePath ctx))
            (finish! []
              ;; Blocks and lots, once the roads are done: a moment's work,
              ;; done between frames so the caption can say so first.
              (swap! state assoc :blocks :pending)
              (set! (.-textContent caption) "cutting blocks and lots...")
              (js/setTimeout
               (fn []
                 (let [{:keys [land roads] :as s} @state
                       t0 (js/performance.now)
                       blocks (lots/blocks (roads/graph roads))
                       ls (lots/lots blocks (:zone land) {:seed (.-seed controls)})]
                   (when (identical? s @state)
                     (swap! state assoc :blocks blocks :lots ls :lot-ms (- (js/performance.now) t0)))))
               30))
            (draw! []
              (let [{:keys [land back roads blocks lots lot-ms]} @state
                    w (.-width canvas) h (.-height canvas)
                    {:keys [scale]} @view
                    size (:size land)
                    [x0 y0] (to-screen [(- size) size])
                    [x1 y1] (to-screen [size (- size)])
                    show (.-show controls)
                    segs (roads/segments roads)]
                (set! (.-fillStyle ctx) "#1a1f28")
                (.fillRect ctx 0 0 w h)
                (set! (.-imageSmoothingEnabled ctx) true)
                (.drawImage ctx back x0 y0 (- x1 x0) (- y1 y0))
                (when (and (vector? blocks) (not= show "roads"))
                  (set! (.-fillStyle ctx) "rgba(70,72,78,0.9)")
                  (doseq [{:keys [points]} blocks] (poly-path! points) (.fill ctx)))
                (when (and (vector? lots) (= show "lots"))
                  (set! (.-strokeStyle ctx) "rgba(30,30,30,0.55)")
                  (set! (.-lineWidth ctx) 0.6)
                  (doseq [{:keys [points zone]} lots]
                    (poly-path! points)
                    (set! (.-fillStyle ctx) (zone-colors zone))
                    (.fill ctx)
                    (when (< scale 2.0) (.stroke ctx))))
                (set! (.-lineCap ctx) "round")
                (doseq [[kind color width] [[:street "#f1efe8" 9.0] [:highway "#f6c94c" 18.0]]]
                  (set! (.-strokeStyle ctx) color)
                  (set! (.-lineWidth ctx) (max (if (= kind :highway) 2.2 1.0) (/ width scale)))
                  (.beginPath ctx)
                  (doseq [{:keys [a b] k :kind} segs :when (= k kind)]
                    (let [[ax ay] (to-screen a) [bx by] (to-screen b)]
                      (.moveTo ctx ax ay) (.lineTo ctx bx by)))
                  (.stroke ctx))
                ;; Bridges, over the top.
                (set! (.-strokeStyle ctx) "#fff7c2")
                (doseq [{:keys [a b bridge?]} segs :when bridge?]
                  (let [[ax ay] (to-screen a) [bx by] (to-screen b)]
                    (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx by) (.stroke ctx)))
                ;; Where the queue's next proposals are waiting.
                (when-not (roads/done? roads)
                  (set! (.-fillStyle ctx) "#ff5a7a")
                  (doseq [[_ {:keys [a]}] (take 60 (:queue roads))]
                    (let [[sx sy] (to-screen a)]
                      (.beginPath ctx) (.arc ctx sx sy 2.5 0 (* 2 math/PI)) (.fill ctx))))
                (set! (.-textContent caption)
                      (str (count segs) " segments (" (count (filter #(= :highway (:kind %)) segs)) " highway), step "
                           (:steps roads) ", " (count (:queue roads)) " waiting"
                           (cond
                             (vector? lots) (str "  |  " (count blocks) " blocks, " (count lots) " lots in "
                                                 (.toFixed lot-ms 0) " ms")
                             (= blocks :pending) "  |  cutting blocks and lots..."
                             :else "")))))
            (frame [now]
              (when @running?
                (js/requestAnimationFrame frame)
                (fit!)
                ;; `speed` steps every sixtieth of a second, however long
                ;; the frames are taking.
                (let [{:keys [land roads blocks]} @state
                      dt (if-let [t0 @last-frame] (min 250.0 (- now t0)) 16.7)
                      steps (max 1 (js/Math.round (* (.-speed controls) (/ dt 16.7))))]
                  (reset! last-frame now)
                  (if (roads/done? roads)
                    (when (nil? blocks) (finish!))
                    (swap! state update :roads
                           (fn [r] (loop [r r k steps]
                                     (if (or (zero? k) (roads/done? r)) r (recur (roads/step land r) (dec k))))))))
                (draw!)))
            (resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h)
                           (or (not= w (.-width canvas)) (not= h (.-height canvas))))
                  (set! (.-width canvas) w)
                  (set! (.-height canvas) h)
                  (swap! view assoc :scale nil))))]
      (let [gui (GUI. #js {:container container})]
        (-> (.add gui controls "seed" 1 99 1) (.onFinishChange regenerate!))
        (-> (.add gui controls "pattern" #js ["mixed" "grid" "radial" "organic" "terrain"]) (.onChange regenerate!))
        (.add gui controls "speed" 1 400 1)
        (.add gui controls "show" #js ["roads" "blocks" "lots"])
        (-> (.add gui controls "backdrop" #js ["population" "elevation"])
            (.onChange (fn [] (swap! state assoc :back (backdrop (:land @state) 520 (.-backdrop controls))))))
        (.add gui #js {:finish (fn [] (swap! state update :roads
                                             (fn [r] (loop [r r] (if (roads/done? r) r (recur (roads/step (:land @state) r)))))))}
              "finish")
        (.add gui #js {:regrow regenerate!} "regrow"))
      ;; Drag to pan, scroll to zoom about the pointer.
      (let [drag (atom nil)]
        (.addEventListener canvas "pointerdown"
                           (fn [^js e]
                             (try (.setPointerCapture canvas (.-pointerId e)) (catch :default _ nil))
                             (reset! drag [(.-clientX e) (.-clientY e)])))
        (.addEventListener canvas "pointermove"
                           (fn [^js e]
                             (when-let [[x y] @drag]
                               (let [s (:scale @view)]
                                 (swap! view #(-> % (update :cx - (* s (- (.-clientX e) x)))
                                                  (update :cy + (* s (- (.-clientY e) y))))))
                               (reset! drag [(.-clientX e) (.-clientY e)]))))
        (doseq [ev ["pointerup" "pointercancel"]]
          (.addEventListener canvas ev (fn [_] (reset! drag nil))))
        (.addEventListener canvas "wheel"
                           (fn [^js e]
                             (.preventDefault e)
                             (let [r (.getBoundingClientRect canvas)
                                   px (- (.-clientX e) (.-left r)) py (- (.-clientY e) (.-top r))
                                   {:keys [scale cx cy]} @view
                                   ;; The point under the pointer stays put.
                                   wx (+ cx (* scale (- px (* 0.5 (.-width canvas)))))
                                   wy (- cy (* scale (- py (* 0.5 (.-height canvas)))))
                                   k (math/exp (* 0.0015 (.-deltaY e)))
                                   s2 (max 0.15 (min 12.0 (* scale k)))]
                               (reset! view {:scale s2
                                             :cx (- wx (* s2 (- px (* 0.5 (.-width canvas)))))
                                             :cy (+ wy (* s2 (- py (* 0.5 (.-height canvas)))))})))
                           #js {:passive false}))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (regenerate!)
      {:start (fn [] (when-not @running? (reset! running? true) (reset! last-frame nil) (resize!)
                               (js/requestAnimationFrame frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "city-plan")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
