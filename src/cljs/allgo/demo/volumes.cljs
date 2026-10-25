(ns allgo.demo.volumes
  "Canvas demo for `allgo.procedural.volume`: hypertexture, fire, steam,
  cumulus, cirrus and stratus, each a density in space, seen by marching
  rays through it.

  Drag to walk around it. `steps` is how many samples a ray takes through
  the volume and `shadows` how many it takes toward the sun from each:
  at 0 shadows a cloud loses all of its shape, which is the clearest
  demonstration there is that a cloud is drawn by the light that does
  not get through it. `time` moves the noise, not the gas: steam rises
  and flames climb because the turbulence is read further along each
  instant, Ebert's way of animating a solid space.

  Every pixel is a ray of a hundred-odd density lookups, on the CPU, so
  the picture is drawn coarse first and sharpened over a few seconds;
  `detail` at `fast` stops a pass short.

  `play` animates it, which is far too slow to do live, so it is done
  the way film is: a loop of `loop` seconds is rendered once, small, a
  frame at a time, and then played back at full speed. Only eight frames
  a second are rendered, and playback fades each into the next, which
  gas -- having no edges to judder -- forgives. Noise in time has no
  period, so the loop's end would jump back to its start; instead a few
  extra frames are rendered past the end and faded into the first ones,
  so the steam keeps rising through the seam. `rate` is how much `time`
  a second of the loop covers, and `frame size` the long side of a frame
  in pixels: a frame of 240 is a second or two of rays."
  (:require [allgo.demo.fps :as fps]
            [allgo.demo.progressive :as progressive]
            [allgo.procedural.volume :as vol]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

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
       :elevation 15
       :loop    2.5
       :rate    0.4
       :frame   240})

(def ^:private budget-ms 24)

(def ^:private loop-fps
  "Frames a second of the rendered loop. Few, because each is seconds of
  work, and enough for gas once neighbors are faded into each other."
  8)

(defn- light []
  ;; Swung about the vertical by `light` degrees from the viewer's right,
  ;; and always well above the horizon.
  (let [a (* (+ (.-azimuth controls) (.-light controls)) (/ math/PI 180))
        x (* 0.7 (math/sin a)) z (* 0.7 (math/cos a)) y 0.7
        l (math/sqrt (+ (* x x) (* y y) (* z z)))]
    [(/ x l) (/ y l) (/ z l)]))

(defn- sky
  "The background: a sky going pale toward the horizon."
  [v]
  (let [t (max 0.0 (min 1.0 (* 0.5 (+ 1 v))))]
    [(+ 0.55 (* -0.30 t)) (+ 0.68 (* -0.22 t)) (+ 0.82 (* -0.10 t))]))

(defn- entry []
  (first (filter #(= (.-volume controls) (:name %)) vol/catalog)))

(defn- sampler [w h time]
  (let [volume ((:make (entry)) {:seed (.-seed controls) :time time})
        size (min w h)
        az (* (.-azimuth controls) (/ math/PI 180))
        el (* (.-elevation controls) (/ math/PI 180))
        opts {:steps (.-steps controls) :shadow-steps (.-shadows controls) :light (light)}]
    (fn [px py]
      (let [u (/ (- (* 2 px) w) size) v (/ (- h (* 2 py)) size)
            bg (sky v)]
        (vol/tone (vol/march volume (vol/camera-ray az el 3.0 0.8 u v)
                             (assoc opts :background bg :jitter (vol/pixel-jitter px py))))))))

(defn- frame-canvas [w h]
  (let [c (js/document.createElement "canvas")]
    (set! (.-width c) w)
    (set! (.-height c) h)
    c))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        tick-fps! (fps/meter! container)
        state (atom nil)
        ;; The loop, while it is being rendered and once it is: the frames
        ;; done, the job on the next, and when playback began.
        anim (atom nil)
        play-button (atom nil)]
    (set! (.-style canvas) "display:block;cursor:grab")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;padding:3px 8px;border-radius:4px;background:rgba(10,14,22,0.7);"
               "font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;pointer-events:none;display:none"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(say! [text]
              (set! (.. caption -style -display) (if text "block" "none"))
              (set! (.-textContent caption) (or text "")))
            (stop-play! []
              (say! nil)
              (when @anim
                (reset! anim nil)
                (when-let [^js b @play-button] (.name b "play"))))
            (restart! []
              (stop-play!)
              (let [w (.-width canvas) h (.-height canvas)]
                (when (and (pos? w) (pos? h))
                  (reset! state (progressive/job ctx w h (sampler w h (.-time controls))
                                                 (if (= "fast" (.-detail controls))
                                                   [16 8 4 2]
                                                   [16 8 4 2 1]))))))
            (next-job [{:keys [fw fh t0 dt frames] :as a}]
              ;; One pass, no coarse preview: the pixels not yet done are
              ;; transparent, so the frame wipes in over the last one.
              (let [c (frame-canvas fw fh)]
                (assoc a :canvas c
                       :job (progressive/job (.getContext c "2d") fw fh
                                             (sampler fw fh (+ t0 (* dt (count frames))))
                                             [1]))))
            (seam [{:keys [n fade frames fw fh]}]
              ;; The loop proper: over its first `fade` frames, the ones
              ;; rendered past its end fade out into them.
              (vec (for [i (range n)]
                     (if (>= i fade)
                       (nth frames i)
                       (let [c (frame-canvas fw fh) cx (.getContext c "2d")]
                         (.drawImage cx (nth frames i) 0 0)
                         (set! (.-globalAlpha cx) (- 1.0 (/ i fade)))
                         (.drawImage cx (nth frames (+ n i)) 0 0)
                         c)))))
            (start-play! []
              (if-not (:moves? (entry))
                (say! (str (.-volume controls) " is a solid, and time does nothing to it"))
                (let [w (.-width canvas) h (.-height canvas)
                      k (/ (.-frame controls) (max w h))
                      n (max 2 (math/round (* loop-fps (.-loop controls))))]
                  (when-let [^js b @play-button] (.name b "stop"))
                  (reset! anim (next-job {:fw (max 1 (math/round (* k w))) :fh (max 1 (math/round (* k h)))
                                          :n n :fade (max 1 (quot n 4))
                                          :t0 (.-time controls)
                                          :dt (/ (.-rate controls) loop-fps)
                                          :frames []})))))
            (show! [^js img]
              (set! (.-imageSmoothingEnabled ctx) true)
              (set! (.-imageSmoothingQuality ctx) "high")
              (.drawImage ctx img 0 0 (.-width canvas) (.-height canvas)))
            (resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h)
                           (or (not= w (.-width canvas)) (not= h (.-height canvas))))
                  (set! (.-width canvas) w)
                  (set! (.-height canvas) h)
                  (restart!))))
            (bake! [{:keys [n fade frames ^js canvas job] :as a}]
              ;; Render the next frame a slice at a time, over the last
              ;; one, and move on when it is done.
              (let [t0 (js/performance.now)
                    done? (progressive/work! job budget-ms)]
                (.putImageData (.getContext canvas "2d") (:img job) 0 0)
                (when-let [prev (peek frames)] (show! prev))
                (show! canvas)
                (tick-fps! (- (js/performance.now) t0))
                (let [frames (cond-> frames done? (conj canvas))]
                  (say! (str "rendering frame " (min (inc (count frames)) (+ n fade)) " of " (+ n fade)
                             " -- it plays once they are all done"))
                  (reset! anim (cond
                                 (not done?) a
                                 (< (count frames) (+ n fade)) (next-job (assoc a :frames frames))
                                 :else (assoc a :frames frames :job nil
                                              :loop (seam (assoc a :frames frames))
                                              :since (js/performance.now)))))))
            (play! [{:keys [n t0 dt since] frames :loop}]
              ;; Between frames i and i + 1, faded from one to the next.
              (let [x (mod (* loop-fps (/ (- (js/performance.now) since) 1000.0)) n)
                    i (long (math/floor x))
                    u (- x i)]
                (show! (nth frames i))
                (set! (.-globalAlpha ctx) u)
                (show! (nth frames (mod (inc i) n)))
                (set! (.-globalAlpha ctx) 1.0)
                (set! (.-time controls) (+ t0 (* dt x)))
                (say! (str "loop of " n " frames   |   time " (.toFixed (+ t0 (* dt x)) 2)))))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (if-let [a @anim]
                  (if (:job a) (bake! a) (play! a))
                  (when-let [{:keys [^js img] :as job} @state]
                    (when-not (progressive/done? job)
                      (let [t0 (js/performance.now)]
                        (progressive/work! job budget-ms)
                        (.putImageData ctx img 0 0)
                        (tick-fps! (- (js/performance.now) t0))))))))]
      (let [gui (GUI. #js {:container container})
            actions #js {:play (fn [] (if @anim (restart!) (start-play!)))}]
        (doto gui
          (-> (.add controls "volume" (clj->js names)) (.onChange restart!))
          (-> (.add controls "seed" 0 99 1) (.onFinishChange restart!))
          (-> (.add controls "time" 0 4 0.05) (.listen) (.onFinishChange restart!))
          (-> (.add controls "steps" 16 160 8) (.onFinishChange restart!))
          (-> (.add controls "shadows" 0 16 1) (.onFinishChange restart!))
          (-> (.add controls "light" -180 180 5) (.onFinishChange restart!))
          (-> (.add controls "detail" #js ["fast" "full"]) (.onChange restart!))
          (-> (.add controls "azimuth" -180 180 1) (.listen) (.onFinishChange restart!))
          (-> (.add controls "elevation" -60 80 1) (.listen) (.onFinishChange restart!)))
        (let [f (.addFolder gui "animation")]
          (reset! play-button (.add f actions "play"))
          (-> (.add f controls "loop" 1 6 0.5) (.name "loop (s)") (.onFinishChange restart!))
          (-> (.add f controls "rate" 0.1 1.5 0.05) (.name "time per second") (.onFinishChange restart!))
          (-> (.add f controls "frame" #js [160 240 360]) (.name "frame size") (.onChange restart!))))
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
