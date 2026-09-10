(ns allgo.demo.lorenz
  "Real-time integration of the Lorenz system, drawn as it is computed.

  The point is to make the integrator visible rather than the attractor.
  Lorenz is chaotic at the classical parameters, so any error is amplified
  rather than averaged away: drop the order or stretch the step and the
  trajectory does not merely blur, it leaves for a different part of the
  attractor. Every method from `allgo.numerics` is on the menu, and
  the readout reports what the integrator is doing while it does it."
  (:require [allgo.demo.fps :as fps]
            [allgo.numerics :as num]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private max-points 120000)
(def ^:private by-name (into {} (map (juxt :name identity)) num/first-order))

(def ^:private controls
  #js {:method         "DOPRI5(4)"
       :dt             0.004
       :stepsPerFrame  6
       :sigma          10.0
       :rho            28.0
       :beta           (/ 8.0 3.0)
       :x0             1.0
       :y0             1.0
       :z0             1.0
       :running        true})

(defn- make-integrator []
  (num/integrator (by-name (.-method controls))
                  (num/lorenz {:sigma (.-sigma controls)
                               :rho   (.-rho controls)
                               :beta  (.-beta controls)})
                  0.0
                  [(.-x0 controls) (.-y0 controls) (.-z0 controls)]
                  (.-dt controls)))

;; The attractor lives around z = 25 and spans some 50 units; shift and
;; shrink it into a comfortable view volume.
(def ^:private view-scale 1.1)
(def ^:private view-shift 25.0)

(defn- line-geometry []
  (doto (THREE/BufferGeometry.)
    (.setAttribute "position" (THREE/BufferAttribute. (js/Float32Array. (* max-points 3)) 3))
    (.setAttribute "color" (THREE/BufferAttribute. (js/Float32Array. (* max-points 3)) 3))))

(defn- compact!
  "Drop the oldest half of the trail when the buffer fills, so the run can
  continue rather than freezing or wrapping into a stray segment across the
  attractor."
  [^js geo n]
  (let [half (quot n 2)]
    (doseq [attr ["position" "color"]]
      (let [^js arr (.-array (.getAttribute geo attr))]
        (.copyWithin arr 0 (* half 3) (* n 3))))
    half))

(defn- push-point! [^js geo i [x y z] speed]
  (let [^js pos (.-array (.getAttribute geo "position"))
        ^js col (.-array (.getAttribute geo "color"))
        b       (* i 3)
        ;; hue by speed: the fast sweeps through the lobes read differently
        ;; from the slow turns at their outer edge
        hue     (/ (min 1.0 (/ speed 120.0)) 1.6)
        c       (.setHSL (THREE/Color.) (- 0.66 hue) 0.85 0.6)]
    (aset pos b (* view-scale x))
    (aset pos (+ b 1) (* view-scale (- z view-shift)))
    (aset pos (+ b 2) (* view-scale y))
    (aset col b (.-r c))
    (aset col (+ b 1) (.-g c))
    (aset col (+ b 2) (.-b c))))

(defn init! [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 55 (/ (.-clientWidth container) (.-clientHeight container)) 0.1 5000)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        geo      (line-geometry)
        trail    (THREE/Line. geo (THREE/LineBasicMaterial. #js {:vertexColors true}))
        head     (THREE/Mesh. (THREE/SphereGeometry. 0.6 16 12)
                              (THREE/MeshBasicMaterial. #js {:color 0xffffff}))
        readout  (js/document.createElement "div")
        running? (atom false)
        tick-fps! (fps/meter! container)
        state    (atom {:integ (make-integrator) :n 0 :steps 0})]
    (set! (.-className readout) "numeric-readout")
    (.appendChild container readout)
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (set! (.-frustumCulled trail) false)
    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
    (.setSize renderer (.-clientWidth container) (.-clientHeight container))
    (set! (.-borderRadius (.-style (.-domElement renderer))) "8px")
    (.appendChild container (.-domElement renderer))
    (.add scene trail)
    (.add scene head)
    (.set (.-position camera) 0 0 95)
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 0 0)
      (letfn [(restart! []
                (swap! state assoc :integ (make-integrator) :n 0)
                (.setDrawRange geo 0 0))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h))))
              (integrate! []
                (let [{:keys [integ n steps]} @state
                      f (:f integ)]
                  (loop [s integ i 0 n n]
                    (if (= i (.-stepsPerFrame controls))
                      (swap! state assoc :integ s :n n
                             :steps (+ steps (.-stepsPerFrame controls)))
                      (let [s'  (num/step s)
                            n'  (if (>= n max-points) (compact! geo n) n)]
                        (if (:accepted s')
                          (let [y (:y s')
                                d (f (:t s') y)
                                speed (Math/sqrt (reduce + (map * d d)))]
                            (push-point! geo n' y speed)
                            (recur s' (inc i) (inc n')))
                          (recur s' (inc i) n')))))))
              (publish! []
                (let [{:keys [integ n]} @state]
                  (.setDrawRange geo 0 n)
                  (set! (.-needsUpdate (.getAttribute geo "position")) true)
                  (set! (.-needsUpdate (.getAttribute geo "color")) true)
                  (when (pos? n)
                    (let [[x y z] (:y integ)]
                      (.set (.-position head) (* view-scale x) (* view-scale (- z view-shift)) (* view-scale y))))
                  (set! (.-textContent readout)
                        (let [stages (:stages (:method integ))
                              evals  (* stages (:steps @state))]
                          (str (.-method controls)
                               "   t=" (.toFixed (:t integ) 2)
                               "   h=" (.toExponential (:h integ) 2)
                               (when-let [e (:error integ)] (str "   err=" (.toExponential e 2)))
                               "\n" stages " f-evals/step"
                               "   " (.toLocaleString evals) " total")))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0 (js/performance.now)]
                    (when (.-running controls) (integrate!))
                    (publish!)
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (let [gui (GUI. #js {:container container})]
          (-> (.add gui controls "method" (clj->js (mapv :name num/first-order)))
              (.onChange restart!))
          (-> (.add gui controls "dt" 0.0005 0.02 0.0005) (.onChange restart!))
          (.add gui controls "stepsPerFrame" 1 40 1)
          (.add gui controls "running")
          (let [system (.addFolder gui "system")]
            (-> (.add system controls "sigma" 1 20 0.1) (.onChange restart!))
            (-> (.add system controls "rho" 1 60 0.1) (.onChange restart!))
            (-> (.add system controls "beta" 0.5 5 0.01) (.onChange restart!))
            (.close system))
          (let [start (.addFolder gui "initial state")]
            (-> (.add start controls "x0" -25 25 0.001) (.onChange restart!))
            (-> (.add start controls "y0" -25 25 0.001) (.onChange restart!))
            (-> (.add start controls "z0" 0 60 0.001) (.onChange restart!))
            (.close start))
          (.add gui #js {:restart restart!} "restart"))
        {:start (fn [] (when-not @running? (reset! running? true) (on-resize) (animate)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "lorenz")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
