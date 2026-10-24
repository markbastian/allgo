(ns allgo.demo.volumes
  "Canvas demo for `allgo.procedural.volume`: hypertexture, fire, steam
  and cloud, each a density in space, seen by marching rays through it.

  Drag to walk around it. `steps` is how many samples a ray takes through
  the volume and `shadows` how many it takes toward the sun from each:
  at 0 shadows a cloud loses all of its shape, which is the clearest
  demonstration there is that a cloud is drawn by the light that does
  not get through it. `time` moves the noise, not the gas: steam rises
  and flames climb because the turbulence is read further along each
  instant, Ebert's way of animating a solid space.

  Every pixel is a ray of a hundred-odd density lookups, on the CPU, so
  the picture is drawn coarse first and sharpened over a few seconds;
  `detail` at `fast` stops a pass short."
  (:require [allgo.demo.fps :as fps]
            [allgo.demo.progressive :as progressive]
            [allgo.procedural.volume :as vol]
            ["lil-gui" :default GUI]))

(def ^:private names (mapv :name vol/catalog))

(def ^:private ^js controls
  #js {:volume  "cumulus"
       :seed    0
       :time    0.0
       :steps   64
       :shadows 6
       :light   40
       :detail  "fast"
       :azimuth 30
       :elevation 15})

(def ^:private budget-ms 24)

(defn- light []
  ;; Swung about the vertical by `light` degrees from the viewer's right,
  ;; and always well above the horizon.
  (let [a (* (+ (.-azimuth controls) (.-light controls)) (/ js/Math.PI 180))
        x (* 0.7 (js/Math.sin a)) z (* 0.7 (js/Math.cos a)) y 0.7
        l (js/Math.hypot x y z)]
    [(/ x l) (/ y l) (/ z l)]))

(defn- sky
  "The background: a sky going pale toward the horizon."
  [v]
  (let [t (max 0.0 (min 1.0 (* 0.5 (+ 1 v))))]
    [(+ 0.55 (* -0.30 t)) (+ 0.68 (* -0.22 t)) (+ 0.82 (* -0.10 t))]))

(defn- sampler [w h]
  (let [{:keys [make]} (first (filter #(= (.-volume controls) (:name %)) vol/catalog))
        volume (make {:seed (.-seed controls) :time (.-time controls)})
        size (min w h)
        az (* (.-azimuth controls) (/ js/Math.PI 180))
        el (* (.-elevation controls) (/ js/Math.PI 180))
        opts {:steps (.-steps controls) :shadow-steps (.-shadows controls) :light (light)}]
    (fn [px py]
      (let [u (/ (- (* 2 px) w) size) v (/ (- h (* 2 py)) size)
            bg (sky v)]
        (vol/tone (vol/march volume (vol/camera-ray az el 3.0 0.8 u v)
                             (assoc opts :background bg :jitter (vol/pixel-jitter px py))))))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        running? (atom false)
        tick-fps! (fps/meter! container)
        state (atom nil)]
    (set! (.-style canvas) "display:block;cursor:grab")
    (.appendChild container canvas)
    (letfn [(restart! []
              (let [w (.-width canvas) h (.-height canvas)]
                (when (and (pos? w) (pos? h))
                  (reset! state (progressive/job ctx w h (sampler w h)
                                                 (if (= "fast" (.-detail controls))
                                                   [16 8 4 2]
                                                   [16 8 4 2 1]))))))
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
                (when-let [{:keys [^js img] :as job} @state]
                  (when-not (progressive/done? job)
                    (let [t0 (js/performance.now)]
                      (progressive/work! job budget-ms)
                      (.putImageData ctx img 0 0)
                      (tick-fps! (- (js/performance.now) t0)))))))]
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "volume" (clj->js names)) (.onChange restart!))
          (-> (.add controls "seed" 0 99 1) (.onFinishChange restart!))
          (-> (.add controls "time" 0 4 0.05) (.onFinishChange restart!))
          (-> (.add controls "steps" 16 160 8) (.onFinishChange restart!))
          (-> (.add controls "shadows" 0 16 1) (.onFinishChange restart!))
          (-> (.add controls "light" -180 180 5) (.onFinishChange restart!))
          (-> (.add controls "detail" #js ["fast" "full"]) (.onChange restart!))
          (-> (.add controls "azimuth" -180 180 1) (.listen) (.onFinishChange restart!))
          (-> (.add controls "elevation" -60 80 1) (.listen) (.onFinishChange restart!))))
      ;; Drag to orbit. Pointer events, so a finger works as a mouse does.
      (let [drag (atom nil)]
        (set! (.. canvas -style -touchAction) "none")
        (.addEventListener canvas "pointerdown"
                           (fn [^js e]
                             (try (.setPointerCapture canvas (.-pointerId e)) (catch :default _ nil))
                             (reset! drag [(.-clientX e) (.-clientY e)])))
        (.addEventListener canvas "pointermove"
                           (fn [^js e]
                             (when-let [[x y] @drag]
                               (set! (.-azimuth controls)
                                     (-> (.-azimuth controls) (- (* 0.4 (- (.-clientX e) x)))
                                         (+ 540) (mod 360) (- 180)))
                               (set! (.-elevation controls)
                                     (-> (.-elevation controls) (+ (* 0.4 (- (.-clientY e) y)))
                                         (max -60) (min 80)))
                               (reset! drag [(.-clientX e) (.-clientY e)])
                               (restart!))))
        (doseq [ev ["pointerup" "pointercancel"]]
          (.addEventListener canvas ev (fn [_] (reset! drag nil)))))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "volumes")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
