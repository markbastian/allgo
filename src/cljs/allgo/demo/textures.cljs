(ns allgo.demo.textures
  "Canvas demo for `allgo.procedural.texture`: brick, marble, wood, stone
  and metal, each a function of a point in space.

  `gallery` shows every texture at once; click one to open it. `view`
  lays it on a plane, or carves a ball out of it -- a solid texture is
  defined everywhere, so the ball needs no wrapping and has no seam.

  `bump` is how far the lighting believes the texture's height: at zero
  every material is a flat picture of itself, which is the clearest way to
  see how much of brick, cobbles or hammered copper is the light across
  their bumps rather than their color. `light` swings the lamp around.

  Nothing is stored and nothing runs on the GPU: every pixel asks the
  texture what is there, three times over for the bumps. That is slow
  enough to watch, so the picture is drawn coarse first and sharpened over
  a few frames."
  (:require [allgo.demo.fps :as fps]
            [allgo.demo.progressive :as progressive]
            [allgo.procedural.texture :as tx]
            ["lil-gui" :default GUI]))

(def ^:private names (mapv :name tx/catalog))

(def ^:private ^js controls
  #js {:texture "gallery"
       :view    "plane"
       :zoom    1.0
       :seed    0
       :bump    0.04
       :light   135})

(def ^:private budget-ms 24)
(def ^:private blocks [8 4 2 1])

(defn- light []
  ;; The lamp's direction: swung about the view axis by `light` degrees,
  ;; 135 being up and to the left, and always partly toward the viewer.
  (let [a (* (.-light controls) (/ js/Math.PI 180))
        x (* 0.78 (js/Math.cos a)) y (* 0.78 (js/Math.sin a)) z 0.62
        l (js/Math.hypot x y z)]
    [(/ x l) (/ y l) (/ z l)]))

(defn- layout
  "Tiles for the gallery: the largest square that fits them all."
  [w h n]
  (let [best (reduce (fn [best cols]
                       (let [rows (js/Math.ceil (/ n cols))
                             s (js/Math.floor (min (/ w cols) (/ h rows)))]
                         (if (> s (:size best)) {:size s :cols cols} best)))
                     {:size 0 :cols 1}
                     (range 1 (inc n)))]
    (assoc best :ox (js/Math.floor (/ (- w (* (:cols best) (:size best))) 2)))))

(def ^:private tile
  "The size of a gallery tile a texture's `:scale` is meant for."
  360)

(def ^:private background [0.02 0.027 0.05])

(defn- sampler
  "`(fn [px py] -> [r g b])` for the current controls, over a `w` by `h`
  canvas."
  [w h]
  (let [seed (.-seed controls)
        opts {:bump (.-bump controls) :light (light)}
        zoom (.-zoom controls)
        sphere? (= "sphere" (.-view controls))
        ;; `span`: how many tiles' worth the square shows, so that a
        ;; texture opened full window keeps the size it had in the
        ;; gallery and shows more of itself, rather than being blown up.
        one (fn [{:keys [make scale]} span]
              (let [tex (make {:seed seed})
                    scale (* scale zoom span)]
                (fn [u v size]
                  ;; u, v in [0, 1) across a square `size` pixels wide.
                  (if sphere?
                    (or (tx/sphere-pixel tex (- (* 2 u) 1) (- 1 (* 2 v)) (* 0.5 scale)
                                         (/ 2 size) opts)
                        background)
                    (tx/plane-pixel tex (* scale (- u 0.5)) (* scale (- 0.5 v)) (/ scale size) opts)))))]
    (if (= "gallery" (.-texture controls))
      (let [{:keys [size cols ox]} (layout w h (count tx/catalog))
            fs (mapv #(one % 1) tx/catalog)]
        (fn [px py]
          (let [c (js/Math.floor (/ (- px ox) size)) r (js/Math.floor (/ py size))
                k (+ c (* r cols))]
            (if (and (<= 0 c) (< c cols) (< k (count fs)))
              ((nth fs k) (/ (- px ox (* c size)) size) (/ (- py (* r size)) size) size)
              background))))
      ;; A plane covers the window; a ball fits inside it.
      (let [size (if sphere? (* 0.94 (min w h)) (max w h))
            f (one (first (filter #(= (.-texture controls) (:name %)) tx/catalog))
                   (max 1 (/ size tile)))]
        (fn [px py] (f (/ (+ px (/ (- size w) 2)) size) (/ (+ py (/ (- size h) 2)) size) size))))))

(defn- draw-labels! [^js ctx w h]
  (when (= "gallery" (.-texture controls))
    (let [{:keys [size cols ox]} (layout w h (count tx/catalog))]
      (set! (.-font ctx) "12px system-ui, sans-serif")
      (set! (.-textBaseline ctx) "bottom")
      (doseq [[k {:keys [name]}] (map-indexed vector tx/catalog)]
        (let [x (+ ox (* size (mod k cols)) 6)
              y (- (* size (inc (quot k cols))) 6)]
          (set! (.-fillStyle ctx) "rgba(0,0,0,0.55)")
          (.fillRect ctx (- x 3) (- y 16) (+ 6 (.-width (.measureText ctx name))) 18)
          (set! (.-fillStyle ctx) "#e8ecf4")
          (.fillText ctx name x y))))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        running? (atom false)
        tick-fps! (fps/meter! container)
        state (atom nil)]
    (set! (.-style canvas) "display:block;cursor:pointer")
    (.appendChild container canvas)
    (letfn [(restart! []
              (let [w (.-width canvas) h (.-height canvas)]
                (when (and (pos? w) (pos? h))
                  (reset! state (progressive/job ctx w h (sampler w h) blocks)))))
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
                (when-let [{:keys [^js img w h] :as job} @state]
                  (when-not (progressive/done? job)
                    (let [t0 (js/performance.now)]
                      (progressive/work! job budget-ms)
                      (.putImageData ctx img 0 0)
                      (draw-labels! ctx w h)
                      (tick-fps! (- (js/performance.now) t0)))))))]
      (let [gui (GUI. #js {:container container})
            picker (.add gui controls "texture" (clj->js (into ["gallery"] names)))]
        (.onChange picker restart!)
        (doto gui
          (-> (.add controls "view" #js ["plane" "sphere"]) (.onChange restart!))
          (-> (.add controls "zoom" 0.25 4 0.05) (.onFinishChange restart!))
          (-> (.add controls "seed" 0 99 1) (.onFinishChange restart!))
          (-> (.add controls "bump" 0 0.15 0.005) (.onFinishChange restart!))
          (-> (.add controls "light" 0 360 5) (.onFinishChange restart!)))
        ;; A click on a gallery tile opens that texture; a click on an open
        ;; one goes back to the gallery.
        (.addEventListener canvas "click"
                           (fn [^js e]
                             (if (= "gallery" (.-texture controls))
                               (let [rect (.getBoundingClientRect canvas)
                                     x (- (.-clientX e) (.-left rect)) y (- (.-clientY e) (.-top rect))
                                     {:keys [size cols ox]} (layout (.-width canvas) (.-height canvas)
                                                                    (count tx/catalog))
                                     c (js/Math.floor (/ (- x ox) size)) r (js/Math.floor (/ y size))
                                     k (+ c (* r cols))]
                                 (when (and (<= 0 c) (< c cols) (< k (count names)))
                                   (.setValue picker (nth names k))))
                               (.setValue picker "gallery")))))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "textures")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
