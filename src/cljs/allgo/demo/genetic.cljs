(ns allgo.demo.genetic
  "Canvas demo for `allgo.procedural.genetic`: breed textures by choosing
  the ones you like.

  Twelve random textures to start. Click the ones worth keeping and
  `breed`: one chosen gives twelve variations on it, several give their
  children. The chosen ones are kept, first in the next screen, so
  nothing liked is lost; `back` returns to the previous screen. Double
  click a texture to see it large, and again to return. The expression
  of the last one clicked is printed underneath -- that is the whole of
  its genome, with its palette.

  Every pixel is the expression evaluated, on the CPU, so a screen is
  drawn coarse first and sharpened over a few seconds."
  (:require [allgo.demo.progressive :as progressive]
            [allgo.procedural.genetic :as genetic]
            [allgo.random :as random]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private n 12)
(def ^:private gap 10)
(def ^:private budget-ms 24)
(def ^:private background [0.06 0.07 0.09])

(def ^:private ^js controls
  #js {:mutation 0.2
       :seed     1})

(defn- layout
  "Tiles for `n` over a `w` by `h` canvas, leaving `top` pixels for the
  title and `bottom` for the caption: as many columns as make the tiles
  largest."
  [w h top bottom]
  (let [h (- h top bottom)
        best (apply max-key :size
                    (for [cols (range 1 (inc n))
                          :let [rows (long (math/ceil (/ n cols)))
                                size (math/floor (min (/ (- w (* gap (inc cols))) cols)
                                                      (/ (- h (* gap (inc rows))) rows)))]]
                      {:cols cols :rows rows :size size}))
        {:keys [cols rows size]} best
        x0 (/ (- w (+ (* cols size) (* gap (dec cols)))) 2)
        y0 (+ top (/ (- h (+ (* rows size) (* gap (dec rows)))) 2))]
    (vec (for [i (range n)]
           {:x (math/floor (+ x0 (* (mod i cols) (+ size gap))))
            :y (math/floor (+ y0 (* (quot i cols) (+ size gap))))
            :size size}))))

(defn- tile-at [tiles px py]
  (first (keep-indexed (fn [i {:keys [x y size]}]
                         (when (and (<= x px) (< px (+ x size)) (<= y py) (< py (+ y size))) i))
                       tiles)))

(defn- shade
  "The texture at a pixel of a tile: the tile spans [-2, 2] on a side, on
  the plane z = 0.3."
  [tex {:keys [x y size]} px py]
  (let [u (- (* 4.0 (/ (- px x) size)) 2.0)
        v (- 2.0 (* 4.0 (/ (- py y) size)))]
    (tex u v 0.3)))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        overlay (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        octx (.getContext overlay "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        state (atom nil)
        job (atom nil)]
    (set! (.-style canvas) "display:block")
    (set! (.-style overlay) "position:absolute;left:0;top:0;cursor:pointer")
    (set! (.-style caption)
          (str "position:absolute;left:12px;right:12px;bottom:8px;max-height:56px;overflow:hidden;"
               "font:11px ui-monospace,Menlo,monospace;color:#9aa6bb;pointer-events:none;white-space:normal"))
    (doto container (.appendChild canvas) (.appendChild overlay) (.appendChild caption))
    (letfn [(textures [population] (mapv genetic/texture population))
            (start-generation! [population gen]
              (swap! state (fn [s]
                             (-> s
                                 (update :history (fnil conj []) (select-keys s [:population :gen]))
                                 (assoc :population population :textures (textures population)
                                        :gen gen :selected #{} :focus nil :shown nil))))
              (redraw!))
            (rng-for [gen] (random/rng (.-seed controls) gen 0 0))
            (fresh! []
              (reset! state nil)
              (start-generation! (genetic/population (rng-for 0) n) 0))
            (breed! []
              (let [{:keys [population selected gen]} @state
                    gen (inc gen)
                    parents (map population (sort selected))]
                (start-generation! (genetic/breed (rng-for gen) parents n (.-mutation controls)) gen)))
            (back! []
              (when-let [prev (peek (:history @state))]
                (when (:population prev)
                  (swap! state (fn [s] (assoc s :population (:population prev)
                                              :textures (textures (:population prev))
                                              :gen (:gen prev) :selected #{} :focus nil :shown nil
                                              :history (pop (:history s)))))
                  (redraw!))))
            (sampler [w h]
              (let [{:keys [textures focus tiles]} @state]
                (if focus
                  (let [size (min w (- h 100))
                        tile {:x (/ (- w size) 2) :y (+ 36 (/ (- h 100 size) 2)) :size size}
                        tex (nth textures focus)]
                    (fn [px py]
                      (if (tile-at [tile] px py) (shade tex tile px py) background)))
                  (fn [px py]
                    (if-let [i (tile-at tiles px py)]
                      (shade (nth textures i) (nth tiles i) px py)
                      background)))))
            (draw-overlay! []
              (let [w (.-width overlay) h (.-height overlay)
                    {:keys [tiles selected focus]} @state]
                (.clearRect octx 0 0 w h)
                (when-not focus
                  (set! (.-lineWidth octx) 3)
                  (set! (.-strokeStyle octx) "#f2c14e")
                  (doseq [i selected
                          :let [{:keys [x y size]} (nth tiles i)]]
                    (.strokeRect octx (- x 2.5) (- y 2.5) (+ size 5) (+ size 5))))))
            (write-caption! []
              (let [{:keys [gen selected population shown focus]} @state]
                (set! (.-textContent caption)
                      (str "generation " gen " · " (count selected) " chosen · "
                           (if focus "double-click to return" "click to choose, double-click to enlarge")
                           (when-let [i (or focus shown)]
                             (str "  —  " (genetic/expr->string (:expr (nth population i)))))))))
            (redraw! []
              (let [w (.-width canvas) h (.-height canvas)]
                (when (and (pos? w) (pos? h))
                  (swap! state assoc :tiles (layout w h 36 64))
                  (reset! job (progressive/job ctx w h (sampler w h) [8 4 2 1]))
                  (draw-overlay!)
                  (write-caption!))))
            (resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h)
                           (or (not= w (.-width canvas)) (not= h (.-height canvas))))
                  (doseq [^js c [canvas overlay]]
                    (set! (.-width c) w)
                    (set! (.-height c) h))
                  (redraw!))))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (when-let [{:keys [^js img] :as j} @job]
                  (when-not (progressive/done? j)
                    (progressive/work! j budget-ms)
                    (.putImageData ctx img 0 0)))))]
      (.addEventListener overlay "click"
                         (fn [^js e]
                           (let [{:keys [tiles focus]} @state]
                             (when-not focus
                               (when-let [i (tile-at tiles (.-offsetX e) (.-offsetY e))]
                                 (swap! state (fn [s]
                                                (-> s
                                                    (update :selected #(if (% i) (disj % i) (conj % i)))
                                                    (assoc :shown i))))
                                 (draw-overlay!)
                                 (write-caption!))))))
      (.addEventListener overlay "dblclick"
                         (fn [^js e]
                           (let [{:keys [tiles focus]} @state]
                             (if focus
                               (swap! state assoc :focus nil)
                               (when-let [i (tile-at tiles (.-offsetX e) (.-offsetY e))]
                                 ;; The two clicks of a double click toggled
                                 ;; the choice twice: it is as it was.
                                 (swap! state assoc :focus i :shown i)))
                             (redraw!))))
      (let [gui (GUI. #js {:container container})]
        (.add gui #js {:breed breed!} "breed")
        (.add gui #js {:back back!} "back")
        (.add gui #js {:newPopulation fresh!} "newPopulation")
        (.add gui controls "mutation" 0.05 0.6 0.05)
        (-> (.add gui controls "seed" 1 99 1) (.onFinishChange fresh!)))
      (.observe (js/ResizeObserver. resize!) container)
      (fresh!)
      (resize!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "genetic")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
