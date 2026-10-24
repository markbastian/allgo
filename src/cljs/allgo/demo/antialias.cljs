(ns allgo.demo.antialias
  "Canvas demo for `allgo.procedural.antialias`: a textured plane running
  off to the horizon, drawn two ways side by side.

  The plane is where aliasing is at its worst, because one picture holds
  every scale at once: a square of the checkerboard covers a hundred
  pixels at the bottom of the screen and a hundred squares fall into one
  pixel at the horizon. Point sampling (`point`) draws the near squares
  perfectly and the far ones as moiré -- false patterns at the wrong
  frequency -- that crawl as the camera moves. The other methods each
  average the texture over the pixel's footprint:

  - `filtered` integrates it, in closed form for the checker and stripes
    and by leaving out too-fine octaves for the noise: one sample, and
    the far field goes smoothly to gray.
  - `grid` and `jittered` take `samples` x `samples` points. The grid's
    moiré moves further out but stays; jitter turns it into noise.
  - `adaptive` takes five, and more only where those disagree --
    `threshold` is how much they may.

  The caption counts what each side spent per pixel. Each method is
  given the same square footprint, its longer side, so all of them blur
  a little along the plane: that is the overblur anisotropic filtering
  exists to remove."
  (:require [allgo.demo.progressive :as progressive]
            [allgo.procedural.antialias :as aa]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private method-names ["point" "filtered" "grid" "jittered" "adaptive"])

(def ^:private ^js controls
  #js {:left      "point"
       :right     "filtered"
       :pattern   "checker"
       :samples   4
       :threshold 0.1
       :scale     1.0
       :tilt      12})

(def ^:private budget-ms 24)
(def ^:private eye-height 1.0)
(def ^:private fov 0.9)

(def ^:private dark [0.13 0.17 0.26])
(def ^:private light [0.93 0.89 0.80])

(defn- tone [v]
  (let [v (min 1.0 (max 0.0 v))]
    (mapv #(+ %1 (* v (- %2 %1))) dark light)))

(defn- sky [t]
  [(+ 0.62 (* -0.25 t)) (+ 0.74 (* -0.2 t)) (+ 0.86 (* -0.08 t))])

(defn- projector
  "`(fn [px py] -> [x z] or nil)`: where the pixel's ray meets the plane,
  in texture units, looking along -z tilted down by `tilt` degrees."
  [w h]
  (let [tilt (math/to-radians (.-tilt controls))
        k (/ (math/tan (* 0.5 fov)) (* 0.5 h))
        scale (.-scale controls)
        c (math/cos tilt) s (math/sin tilt)]
    (fn [px py]
      (let [u (* k (- px (* 0.5 w)))
            v (* k (- (* 0.5 h) py))
            ;; Forward (0, -s, -c), up (0, c, -s), right (1, 0, 0).
            dy (- (* v c) s)
            dz (- (+ c (* v s)))]
        (when (neg? dy)
          (let [t (/ eye-height (- dy))]
            [(* scale t u) (* scale t dz)]))))))

(defn- sampler [w h ^js stats]
  (let [project (projector w h)
        {:keys [point filtered]} (aa/patterns (.-pattern controls))
        n (.-samples controls)
        threshold (.-threshold controls)
        method-of (fn [px] (if (< px (* 0.5 w)) (.-left controls) (.-right controls)))]
    (fn [px py]
      (if-let [[x z] (project px py)]
        (let [;; The footprint: how far the plane point moves to the next
              ;; pixel across and down, the longer of the two.
              [x1 z1] (or (project (inc px) py) [x z])
              [x2 z2] (or (project px (inc py)) [##Inf ##Inf])
              fw (min 1e6 (max (math/hypot (- x1 x) (- z1 z)) (math/hypot (- x2 x) (- z2 z))))
              side (if (< px (* 0.5 w)) "l" "r")
              [v cost] (case (method-of px)
                         "point" [(point x z) 1]
                         "filtered" [(filtered x z fw) 1]
                         "grid" [(aa/supersample point x z fw n :grid [px py]) (* n n)]
                         "jittered" [(aa/supersample point x z fw n :jittered [px py]) (* n n)]
                         "adaptive" (aa/adaptive point x z fw threshold 3))]
          (aset stats side (+ (aget stats side) cost))
          (aset stats (str side "n") (inc (aget stats (str side "n"))))
          (tone v))
        (sky (min 1.0 (/ (- (* 0.5 h) py) (* 0.5 h))))))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        job (atom nil)
        stats #js {:l 0 :ln 0 :r 0 :rn 0}]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;padding:3px 8px;border-radius:4px;background:rgba(10,14,22,0.7);"
               "font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;pointer-events:none"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (let [w (.-width canvas) h (.-height canvas)]
                (when (and (pos? w) (pos? h))
                  (doseq [k ["l" "ln" "r" "rn"]] (aset stats k 0))
                  (reset! job (progressive/job ctx w h (sampler w h stats) [8 4 2 1])))))
            (per-pixel [k] (let [n (aget stats (str k "n"))]
                             (if (pos? n) (.toFixed (/ (aget stats k) n) 1) "-")))
            (resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h)
                           (or (not= w (.-width canvas)) (not= h (.-height canvas))))
                  (set! (.-width canvas) w)
                  (set! (.-height canvas) h)
                  (restart!))))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (when-let [{:keys [^js img] :as j} @job]
                  (when-not (progressive/done? j)
                    (progressive/work! j budget-ms)
                    (.putImageData ctx img 0 0)
                    ;; The divide, drawn over the picture.
                    (let [x (math/floor (* 0.5 (.-width canvas)))]
                      (set! (.-fillStyle ctx) "#f2c14e")
                      (.fillRect ctx (dec x) 0 2 (.-height canvas)))
                    (set! (.-textContent caption)
                          (str (.-left controls) ": " (per-pixel "l") " samples a pixel   |   "
                               (.-right controls) ": " (per-pixel "r") " samples a pixel"))))))]
      (let [gui (GUI. #js {:container container})]
        (-> (.add gui controls "left" (clj->js method-names)) (.onChange restart!))
        (-> (.add gui controls "right" (clj->js method-names)) (.onChange restart!))
        (-> (.add gui controls "pattern" (clj->js (sort (keys aa/patterns)))) (.onChange restart!))
        (-> (.add gui controls "samples" 2 8 1) (.onFinishChange restart!))
        (-> (.add gui controls "threshold" 0.01 0.5 0.01) (.onFinishChange restart!))
        (-> (.add gui controls "scale" 0.25 4 0.25) (.onFinishChange restart!))
        (-> (.add gui controls "tilt" 2 40 1) (.onFinishChange restart!)))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "antialias")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
