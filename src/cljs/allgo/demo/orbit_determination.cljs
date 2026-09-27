(ns allgo.demo.orbit-determination
  "Recovering an orbit from noisy tracking data (Montenbruck & Gill ch. 8),
  and with it every other chapter at once.

  What is animated here is not time but *convergence*. The true orbit is
  drawn in green and never moves. The estimate starts kilometers away in
  blue, knowing only a pile of ranges corrupted by noise, and each
  Gauss-Newton iteration snaps it closer until the two are indistinguishable
  and the residuals have fallen to the noise floor.

  The chart is the honest part. Position error is only knowable because this
  is a simulation with a truth to compare against; in real work all you ever
  see is the residual, which is why the two are plotted together.

  They do not always fall in step, and that is the lesson worth staying for.
  Shorten the arc to three hours and only two stations ever see the
  satellite: eleven ranges for six unknowns. The residuals still drop to the
  noise floor -- the fit is excellent -- while the true error climbs to
  nine kilometers. A good fit is not a good orbit. The covariance is what
  tells you which you have, reporting a formal uncertainty of some seven
  kilometers, and the readout says so."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.estimation :as est]
            [allgo.astro.frames :as fr]
            [allgo.astro.geodesy :as gd]
            [allgo.astro.kepler :as kep]
            [allgo.astro.time :as atime]
            [allgo.astro.variational :as va]
            [allgo.demo.fps :as fps]
            [allgo.numerics :as num]
            [allgo.numerics.linear :as lin]
            [allgo.random :as random]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private scale (/ 20.0 c/R-earth))
(def ^:private mu c/GM-earth)
(def ^:private elevation-mask (* 10.0 c/degrees))

(def ^:private ^js controls
  #js {:noiseMeters   10.0
       :arcHours      6.0
       :sampleSeconds 150.0
       :startErrorKm  2.7
       :autoIterate   true})

(def ^:private site-list
  [[35.0 -117.0] [-25.0 28.0] [40.0 140.0] [-33.0 151.0] [51.0 0.0]])

(def ^:private stations
  (mapv (fn [[la lo]] (gd/geodetic->cartesian (* la c/degrees) (* lo c/degrees) 0.0)) site-list))

(defn- station-eci [s secs]
  (let [mjd (+ c/mjd-J2000 (/ secs 86400.0))]
    (lin/mat-vec (fr/terrestrial->celestial (atime/utc->tt mjd) mjd) s)))

(defn- accel [_ r _]
  (let [d (Math/sqrt (reduce + (map * r r)))]
    (mapv #(* (- (/ mu (* d d d))) %) r)))

(def ^:private rk4 (first (filter #(= "RK4" (:name %)) num/first-order)))

(defn- arc
  "Propagate once, sampling state and transition matrix at each time. Done
  in one pass rather than restarting from the epoch per observation, which
  would make the cost quadratic in the length of the arc."
  [x0 times]
  (loop [integ (num/integrator rk4 (va/rhs accel) 0.0
                               (va/initial (subvec (vec x0) 0 3) (subvec (vec x0) 3 6))
                               20.0 {:adaptive? false})
         ts times
         out []]
    (if (empty? ts)
      out
      (let [s (num/step-until integ (first ts))]
        (recur s (rest ts) (conj out (va/unpack (:y s))))))))

(defn- setup
  "Build a truth orbit, simulate observations from it, and start the estimate
  somewhere wrong."
  []
  (let [sigma  (/ (.-noiseMeters controls) 1000.0)
        span   (* 3600.0 (.-arcHours controls))
        times  (vec (range 60.0 span (.-sampleSeconds controls)))
        truth  (let [[r v] (kep/elements->state mu {:a 7500.0 :e 0.02 :i 0.95
                                                    :raan 1.1 :argp 2.0 :nu 0.4})]
                 (vec (concat r v)))
        t-arc  (arc truth times)
        obs    (vec (for [[i secs] (map-indexed vector times)
                          [idx s]  (map-indexed vector stations)
                          :let [rs (:r (nth t-arc i))
                                se (station-eci s secs)]
                          :when (> (:elevation (gd/look-angles se rs)) elevation-mask)]
                      {:i i :t secs :station idx
                       :measured (+ (Math/sqrt (reduce + (map * (mapv - rs se) (mapv - rs se))))
                                    (* sigma (random/gaussian rand 0.0 1.0)))}))
        e      (.-startErrorKm controls)
        guess  (mapv + truth [(* e 0.74) (* e -0.56) (* e 0.37)
                              (* e 7.4e-4) (* e 3.7e-4) (* e -5.6e-4)])]
    {:times times :truth truth :truth-arc t-arc :obs obs :sigma sigma
     :x guess :iter 0 :history [] :arc (arc guess times)}))

(defn- rows-for [state]
  (let [{:keys [obs sigma arc]} state]
    (mapv (fn [{:keys [i t station measured]}]
            (let [{:keys [r phi]} (nth arc i)
                  d   (mapv - r (station-eci (nth stations station) t))
                  rho (Math/sqrt (reduce + (map * d d)))
                  los (mapv #(/ % rho) d)]
              {:H (mapv (fn [j] (reduce + (map-indexed
                                           (fn [ii u] (* u (nth (nth phi ii) j))) los)))
                        (range 6))
               :residual (- measured rho)
               :weight (/ 1.0 (* sigma sigma))}))
          obs)))

(defn- iterate-once
  "One Gauss-Newton step. Returns the state unchanged if the normal matrix is
  not positive definite -- which means the observations leave some direction
  of the orbit undetermined, and no amount of solving will invent it."
  [state]
  (let [rows (rows-for state)
        rms  (est/rms rows)
        err  (Math/sqrt (reduce + (map (fn [a b] (let [d (- a b)] (* d d)))
                                       (subvec (:x state) 0 3) (subvec (:truth state) 0 3))))
        sol  (est/solve rows 6)
        hist (conj (:history state) {:rms rms :error err})]
    (if-not sol
      (assoc state :history hist :singular? true)
      (let [x' (mapv + (:x state) (:correction sol))]
        (assoc state :x x' :arc (arc x' (:times state))
               :iter (inc (:iter state)) :history hist
               :covariance (:covariance sol))))))

;; ------------------------------------------------------------------ drawing

(defn- trail-geometry [n]
  (doto (THREE/BufferGeometry.)
    (.setAttribute "position" (THREE/BufferAttribute. (js/Float32Array. (* n 3)) 3))))

(defn- fill-trail! [^js geo states]
  (let [arr (.-array (.getAttribute geo "position"))]
    (doseq [[i {[x y z] :r}] (map-indexed vector states)]
      (aset arr (* i 3) (* scale x))
      (aset arr (+ (* i 3) 1) (* scale z))
      (aset arr (+ (* i 3) 2) (* scale y)))
    (.setDrawRange geo 0 (count states))
    (set! (.-needsUpdate (.getAttribute geo "position")) true)))

(def ^:private chart-w 230)
(def ^:private chart-h 84)

(defn- draw-chart! [^js ctx history]
  (.clearRect ctx 0 0 chart-w chart-h)
  (set! (.-fillStyle ctx) "rgba(5,7,13,0.75)")
  (.fillRect ctx 0 0 chart-w chart-h)
  (set! (.-strokeStyle ctx) "rgba(255,255,255,0.10)")
  (doseq [d (range 0 5)]
    (let [y (* chart-h (/ d 4.0))]
      (.beginPath ctx) (.moveTo ctx 0 y) (.lineTo ctx chart-w y) (.stroke ctx)))
  ;; both on a shared log scale, decades from 1e-4 to 1e2 km
  (let [lo -4.0 hi 2.0
        plot (fn [key color]
               (when (> (count history) 1)
                 (set! (.-strokeStyle ctx) color)
                 (set! (.-lineWidth ctx) 2)
                 (.beginPath ctx)
                 (doseq [[i h] (map-indexed vector history)]
                   (let [v (js/Math.log10 (max 1e-9 (get h key)))
                         x (* chart-w (/ i (max 1 (dec (count history)))))
                         y (* chart-h (- 1.0 (/ (- v lo) (- hi lo))))]
                     (if (zero? i) (.moveTo ctx x y) (.lineTo ctx x y))))
                 (.stroke ctx)))]
    (plot :rms "#6c8cff")
    (plot :error "#3fd07a"))
  (set! (.-fillStyle ctx) "rgba(148,148,171,0.9)")
  (set! (.-font ctx) "9px ui-monospace, Menlo, monospace")
  (.fillText ctx "log km" 6 11)
  (set! (.-fillStyle ctx) "#6c8cff") (.fillText ctx "residual" 6 (- chart-h 14))
  (set! (.-fillStyle ctx) "#3fd07a") (.fillText ctx "true error" 6 (- chart-h 4)))

(defn init! [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 50 (/ (.-clientWidth container) (.-clientHeight container)) 0.1 5000)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        truth-geo (trail-geometry 4000)
        est-geo   (trail-geometry 4000)
        truth-line (THREE/Line. truth-geo (THREE/LineBasicMaterial. #js {:color 0x3fd07a}))
        est-line   (THREE/Line. est-geo (THREE/LineBasicMaterial.
                                         #js {:color 0x6c8cff :transparent true :opacity 0.9}))
        earth    (THREE/LineSegments.
                  (THREE/WireframeGeometry. (THREE/SphereGeometry. (* scale c/R-earth) 24 16))
                  (THREE/LineBasicMaterial. #js {:color 0x2a4a7a :transparent true :opacity 0.3}))
        site-dots (THREE/Points.
                   (doto (THREE/BufferGeometry.)
                     (.setAttribute "position" (THREE/BufferAttribute. (js/Float32Array. (* 5 3)) 3)))
                   (THREE/PointsMaterial. #js {:color 0xffd479 :size 1.2}))
        readout  (js/document.createElement "div")
        chart    (js/document.createElement "canvas")
        ctx      (.getContext chart "2d")
        running? (atom false)
        tick-fps! (fps/meter! container)
        state    (atom (setup))
        last-solve (atom 0)]
    (set! (.-className readout) "numeric-readout")
    (set! (.-className chart) "numeric-chart")
    (set! (.-width chart) chart-w)
    (set! (.-height chart) chart-h)
    (.appendChild container readout)
    (.appendChild container chart)
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (doseq [^js o [truth-line est-line]] (set! (.-frustumCulled o) false))
    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
    (.setSize renderer (.-clientWidth container) (.-clientHeight container))
    (set! (.-borderRadius (.-style (.-domElement renderer))) "8px")
    (.appendChild container (.-domElement renderer))
    (doseq [o [truth-line est-line earth site-dots]] (.add scene o))
    (.set (.-position camera) 40 26 40)
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 0 0)
      (letfn [(refresh! []
                (fill-trail! truth-geo (:truth-arc @state))
                (fill-trail! est-geo (:arc @state))
                (let [arr (.-array (.getAttribute (.-geometry site-dots) "position"))]
                  (doseq [[i s] (map-indexed vector stations)]
                    (let [[x y z] (station-eci s 0.0)]
                      (aset arr (* i 3) (* scale x))
                      (aset arr (+ (* i 3) 1) (* scale z))
                      (aset arr (+ (* i 3) 2) (* scale y))))
                  (set! (.-needsUpdate (.getAttribute (.-geometry site-dots) "position")) true)))
              (restart! [] (reset! state (setup)) (refresh!))
              (step! [] (swap! state iterate-once) (refresh!))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h))))
              (publish! []
                (let [{:keys [obs iter history covariance singular?]} @state
                      last* (peek history)
                      formal (when covariance
                               (Math/sqrt (+ (nth (nth covariance 0) 0)
                                             (nth (nth covariance 1) 1)
                                             (nth (nth covariance 2) 2))))
                      ;; The dangerous case is not a failed solve but a
                      ;; successful one nobody should trust: residuals at the
                      ;; noise floor, and a covariance quietly saying the
                      ;; orbit is undetermined.
                      weak? (and formal last* (> formal (* 20.0 (:rms last*))))]
                  (draw-chart! ctx history)
                  (set! (.-textContent readout)
                        (str "iteration " iter "   " (count obs) " ranges from "
                             (count (distinct (map :station obs))) " stations"
                             "\\n"
                             (if singular?
                               "normal matrix not positive definite: orbit undetermined"
                               (str "residual "
                                    (if last* (.toFixed (* 1000.0 (:rms last*)) 1) "-") " m"
                                    "   true error "
                                    (if last* (.toFixed (* 1000.0 (:error last*)) 1) "-") " m"
                                    (if formal (str "   formal " (.toFixed (* 1000.0 formal) 1) " m") "")
                                    (when weak?
                                      "\nweakly determined: good fit, untrustworthy orbit")))))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0 (js/performance.now)]
                    ;; The solve is far too heavy for a frame, so it runs at
                    ;; most once a second and only while it is still moving.
                    (when (and (.-autoIterate controls)
                               (< (:iter @state) 6)
                               (> (- t0 @last-solve) 1000.0))
                      (reset! last-solve t0)
                      (step!))
                    (publish!)
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (refresh!)
        (let [gui (GUI. #js {:container container})]
          (-> (.add gui controls "noiseMeters" 0 200 1) (.onChange restart!))
          ;; Down to one hour, so the ill-conditioned failure is reachable
          ;; rather than merely described.
          (-> (.add gui controls "arcHours" 1 10 0.5) (.onChange restart!))
          (-> (.add gui controls "sampleSeconds" 60 300 10) (.onChange restart!))
          (-> (.add gui controls "startErrorKm" 0.1 20 0.1) (.onChange restart!))
          (.add gui controls "autoIterate")
          (.add gui #js {:iterate step!} "iterate")
          (.add gui #js {:restart restart!} "restart"))
        {:start (fn [] (when-not @running? (reset! running? true) (on-resize) (animate)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "orbit-determination")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
