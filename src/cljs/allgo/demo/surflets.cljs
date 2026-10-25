(ns allgo.demo.surflets
  "Canvas demo for `allgo.procedural.surflet`: a plane of gradient noise,
  shown as what it is a sum of.

  Every lattice point carries an arrow, its gradient, and owns one
  surflet: a ramp up the arrow, faded to nothing a cell away. Hover a
  point to see its surflet's reach outlined, and `show` it alone; drag
  from a point to turn its arrow, and only the cells around it change --
  noise edited by hand, which works because a surflet is local. With
  `kernel` at `radial` the hump is round instead of square, and the
  caption counts the price: a round hump has to reach further to cover
  the middles of the cells, so more surflets add up at every point.

  The picture is computed here from `surflet/hump` and the lattice's
  gradients, over flat arrays, rather than by calling `surflet/basis` at
  every pixel: the sum is the same (the tests hold the basis to it), and
  this way a drag redraws while the pointer moves."
  (:require [allgo.procedural.noise :as noise]
            [allgo.procedural.surflet :as surflet]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls
  #js {:kernel "separable"
       :radius 1.5
       :show   "sum"
       :seed   3
       :arrows true})

(def ^:private dark [0.15 0.33 0.62])
(def ^:private mid [0.10 0.11 0.14])
(def ^:private light [0.98 0.70 0.30])

(defn- color [v]
  (let [t (math/tanh (* 1.4 v))
        [a b] (if (neg? t) [mid dark] [mid light])
        t (abs t)]
    (mapv #(+ %1 (* t (- %2 %1))) a b)))

(defn- opts []
  {:seed (.-seed controls)
   :kernel (keyword (.-kernel controls))
   :radius (.-radius controls)
   :dimensions 2})

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        overrides (atom {})
        hover (atom nil)
        drag (atom nil)
        dirty? (atom true)
        running? (atom false)]
    (set! (.-style canvas) "display:block;cursor:crosshair;touch-action:none")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;padding:3px 8px;border-radius:4px;background:rgba(10,14,22,0.7);"
               "font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;pointer-events:none"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(layout []
              ;; Lattice units to pixels: about seven cells across the
              ;; shorter side, lattice point (0, 0) at the middle.
              (let [w (.-width canvas) h (.-height canvas)
                    cell (/ (min w h) 6.5)]
                {:w w :h h :cell cell :ox (* 0.5 w) :oy (* 0.5 h)}))
            (to-lattice [{:keys [cell ox oy]} px py]
              [(/ (- px ox) cell) (/ (- oy py) cell)])
            (to-screen [{:keys [cell ox oy]} x y]
              [(+ ox (* x cell)) (- oy (* y cell))])
            (gradients [{:keys [w h] :as l} n]
              ;; Every lattice point that can reach the screen, and its
              ;; gradient's x and y: the plane is z = 0, so z drops out.
              (let [[x0 y1] (to-lattice l 0 0) [x1 y0] (to-lattice l w h)
                    i0 (- (long (math/floor x0)) n) i1 (+ (long (math/ceil x1)) n)
                    j0 (- (long (math/floor y0)) n) j1 (+ (long (math/ceil y1)) n)
                    p (noise/permutation (.-seed controls))
                    nx (inc (- i1 i0)) ny (inc (- j1 j0))
                    gx (js/Float64Array. (* nx ny)) gy (js/Float64Array. (* nx ny))]
                (doseq [i (range i0 (inc i1)) j (range j0 (inc j1))]
                  (let [[x y] (surflet/gradient-at p @overrides [i j 0])
                        k (+ (- i i0) (* nx (- j j0)))]
                    (aset gx k x) (aset gy k y)))
                {:i0 i0 :j0 j0 :nx nx :ny ny :gx gx :gy gy}))
            (render! []
              (let [{:keys [w h cell] :as l} (layout)
                    radial? (= "radial" (.-kernel controls))
                    radius (.-radius controls)
                    n (surflet/reach (opts))
                    gain (if radial? 0.45 1.0)
                    {:keys [i0 j0 nx gx gy]} (gradients l n)
                    only (when (= "hovered" (.-show controls)) @hover)
                    step 2
                    img (.createImageData ctx w h)
                    data (.-data img)]
                (doseq [py (range 0 h step) px (range 0 w step)]
                  (let [[x y] (to-lattice l (+ px 0.5) (+ py 0.5))
                        fi (long (math/floor x)) fj (long (math/floor y))
                        v (loop [di (- 1 n) dj (- 1 n) acc 0.0]
                            (cond
                              (> dj n) acc
                              (> di n) (recur (- 1 n) (inc dj) acc)
                              :else
                              (let [i (+ fi di) j (+ fj dj)
                                    dx (- x i) dy (- y j)
                                    wt (if radial?
                                         (surflet/hump (/ (math/sqrt (+ (* dx dx) (* dy dy))) radius))
                                         (* (surflet/hump dx) (surflet/hump dy)))
                                    k (+ (- i i0) (* nx (- j j0)))
                                    on? (or (nil? only) (= only [i j]))]
                                (recur (inc di) dj
                                       (if (and on? (pos? wt))
                                         (+ acc (* wt (+ (* (aget gx k) dx) (* (aget gy k) dy))))
                                         acc)))))
                        [r g b] (color (* gain v))]
                    (doseq [sy (range py (min h (+ py step))) sx (range px (min w (+ px step)))]
                      (let [o (* 4 (+ sx (* sy w)))]
                        (aset data o (* 255 r)) (aset data (+ o 1) (* 255 g))
                        (aset data (+ o 2) (* 255 b)) (aset data (+ o 3) 255)))))
                (.putImageData ctx img 0 0)
                ;; The lattice: a dot at every point, an arrow for its
                ;; gradient, and the hovered one's reach.
                (let [[x0 y1] (to-lattice l 0 0) [x1 y0] (to-lattice l w h)]
                  (doseq [i (range (long (math/floor x0)) (inc (long (math/ceil x1))))
                          j (range (long (math/floor y0)) (inc (long (math/ceil y1))))]
                    (let [[sx sy] (to-screen l i j)
                          k (+ (- i i0) (* nx (- j j0)))
                          ax (aget gx k) ay (aget gy k)
                          edited? (contains? @overrides [i j 0])
                          hot? (= @hover [i j])]
                      (set! (.-fillStyle ctx) (cond hot? "#ffffff" edited? "#f2c14e" :else "rgba(230,236,245,0.75)"))
                      (.beginPath ctx) (.arc ctx sx sy (if hot? 4 2.5) 0 (* 2 math/PI)) (.fill ctx)
                      (when (.-arrows controls)
                        (let [len (* 0.28 cell) ex (+ sx (* len ax)) ey (- sy (* len ay))
                              a (math/atan2 (- ey sy) (- ex sx))]
                          (set! (.-strokeStyle ctx) (if edited? "#f2c14e" "rgba(230,236,245,0.75)"))
                          (set! (.-lineWidth ctx) 1.5)
                          (.beginPath ctx) (.moveTo ctx sx sy) (.lineTo ctx ex ey)
                          (.moveTo ctx (- ex (* 6 (math/cos (- a 0.5)))) (- ey (* 6 (math/sin (- a 0.5)))))
                          (.lineTo ctx ex ey)
                          (.lineTo ctx (- ex (* 6 (math/cos (+ a 0.5)))) (- ey (* 6 (math/sin (+ a 0.5)))))
                          (.stroke ctx)))
                      (when hot?
                        (set! (.-strokeStyle ctx) "rgba(255,255,255,0.8)")
                        (set! (.-lineWidth ctx) 1)
                        (.setLineDash ctx #js [4 4])
                        (.beginPath ctx)
                        (if radial?
                          (.arc ctx sx sy (* radius cell) 0 (* 2 math/PI))
                          (.rect ctx (- sx cell) (- sy cell) (* 2 cell) (* 2 cell)))
                        (.stroke ctx)
                        (.setLineDash ctx #js [])))))
                (set! (.-textContent caption)
                      (str (.-kernel controls) ": " (* 4 n n) " surflets visited a point here, "
                           (* 8 n n n) " in 3D   |   "
                           (count @overrides) " edited   |   drag from a lattice point to turn its gradient"))))
            (nearest [px py]
              (let [l (layout)
                    [x y] (to-lattice l px py)
                    i (math/round x) j (math/round y)
                    [sx sy] (to-screen l i j)]
                (when (< (math/hypot (- px sx) (- py sy)) (max 12 (* 0.25 (:cell l))))
                  [i j])))
            (resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h)
                           (or (not= w (.-width canvas)) (not= h (.-height canvas))))
                  (set! (.-width canvas) w)
                  (set! (.-height canvas) h)
                  (reset! dirty? true))))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (when @dirty?
                  (reset! dirty? false)
                  (render!))))
            (local [^js e]
              (let [r (.getBoundingClientRect canvas)]
                [(- (.-clientX e) (.-left r)) (- (.-clientY e) (.-top r))]))]
      (let [gui (GUI. #js {:container container})
            redraw! #(reset! dirty? true)]
        (-> (.add gui controls "kernel" #js ["separable" "radial"]) (.onChange redraw!))
        (-> (.add gui controls "radius" 1 2 0.05) (.onChange redraw!))
        (-> (.add gui controls "show" #js ["sum" "hovered"]) (.onChange redraw!))
        (-> (.add gui controls "seed" 0 99 1) (.onChange (fn [] (reset! overrides {}) (redraw!))))
        (-> (.add gui controls "arrows") (.onChange redraw!))
        (.add gui #js {:resetEdits (fn [] (reset! overrides {}) (redraw!))} "resetEdits"))
      (.addEventListener canvas "pointerdown"
                         (fn [^js e]
                           (let [[px py] (local e)]
                             (when-let [ij (nearest px py)]
                               (try (.setPointerCapture canvas (.-pointerId e)) (catch :default _ nil))
                               (reset! drag ij)))))
      (.addEventListener canvas "pointermove"
                         (fn [^js e]
                           (let [[px py] (local e)]
                             (if-let [[i j] @drag]
                               (let [l (layout)
                                     [sx sy] (to-screen l i j)
                                     dx (- px sx) dy (- sy py)
                                     len (math/hypot dx dy)]
                                 (when (> len 4)
                                   ;; As long as the lattice's own: the twelve
                                   ;; gradients are all root 2.
                                   (let [s (/ (math/sqrt 2.0) len)]
                                     (swap! overrides assoc [i j 0] [(* s dx) (* s dy) 0.0])
                                     (reset! dirty? true))))
                               (let [h (nearest px py)]
                                 (when (not= h @hover)
                                   (reset! hover h)
                                   (reset! dirty? true)))))))
      (doseq [ev ["pointerup" "pointercancel"]]
        (.addEventListener canvas ev (fn [_] (reset! drag nil))))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      {:start (fn [] (when-not @running? (reset! running? true) (reset! dirty? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "surflets")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
