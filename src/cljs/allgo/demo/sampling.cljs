(ns allgo.demo.sampling
  "Canvas demo for `allgo.numerics.quasirandom`: sampling the unit square
  more evenly than chance, after LaValle, *Planning Algorithms*, Cambridge
  University Press, 2006, sections 5.2.3 and 5.2.4.

  Four ways of placing the same number of points, a point at a time:
  uniformly at random, which leaves clumps and holes; Halton's sequence,
  each new point in the largest gap the others leave; Hammersley's set,
  which knows the count in advance and does better still; and Sukharev's
  grid of cell centers, the best of all by dispersion but only at square
  counts. Each panel draws its largest empty square -- the dispersion,
  how far a point of the square can be from the nearest sample under the
  max norm -- and the charts below follow dispersion and star discrepancy
  (the worst gap between a box's share of the points and its area) as
  the points come."
  (:require [allgo.numerics.quasirandom :as qr]
            [allgo.random :as random]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:seed 1 :speed 2 :most 400})

(def ^:private kinds ["random" "Halton" "Hammersley" "Sukharev"])
(def ^:private colors {"random" "#ff8f70" "Halton" "#ffd479" "Hammersley" "#7dffa0" "Sukharev" "#7fd1ff"})

(def ^:private probes (vec (for [i (range 41) j (range 41)] [(/ i 40.0) (/ j 40.0)])))
(def ^:private corners (vec (for [i (range 1 21) j (range 1 21)] [(/ i 20.0) (/ j 20.0)])))

(defn- points
  "The first `n` of each kind; random from `pool`, a long draw."
  [kind n pool]
  (case kind
    "random" (subvec pool 0 n)
    "Halton" (qr/halton n 2)
    "Hammersley" (qr/hammersley n 2)
    "Sukharev" (let [k (long (math/floor (math/sqrt n)))] (if (pos? k) (qr/sukharev k 2) []))))

(defn- largest-gap
  "The probe farthest from every point, and that distance (max norm)."
  [pts]
  (apply max-key second
         (for [q probes]
           [q (reduce min ##Inf (map (fn [p] (max (abs (- (p 0) (q 0))) (abs (- (p 1) (q 1))))) pts))])))

(defn- measures [pool n]
  (into {} (for [kind kinds :let [pts (points kind n pool)] :when (seq pts)]
             [kind {:points pts :gap (largest-gap pts) :discrepancy (qr/star-discrepancy pts corners)}])))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        model (atom nil)]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (reset! model {:pool (qr/uniform 2000 2 (random/rng (.-seed controls))) :n 1 :t 0 :history []}))
            (advance! []
              (swap! model (fn [{:keys [n t pool history] :as m}]
                             (let [t (+ t (.-speed controls))
                                   n' (min (.-most controls) (inc (quot t 6)))]
                               (if (and (> n' n) (<= n' (.-most controls)))
                                 (let [ms (measures pool n')]
                                   (assoc m :t t :n n' :current ms
                                          :history (conj history [n' (into {} (map (fn [[k {:keys [gap discrepancy]}]] [k [(second gap) discrepancy]]) ms))])))
                                 (assoc m :t t))))))
            (text! [s x y color & [size]]
              (set! (.-fillStyle ctx) color) (set! (.-font ctx) (str (or size 12) "px ui-monospace,Menlo,monospace"))
              (.fillText ctx s x y))
            (chart! [x y w h history idx label top]
              (set! (.-fillStyle ctx) "#161b25") (.fillRect ctx x y w h)
              (text! label x (- y 6) "#9fb4d8" 11)
              (let [most (.-most controls)
                    cx (fn [n] (+ x (* w (/ (math/log n) (math/log most)))))
                    cy (fn [v] (- (+ y h) (* h (/ (min v top) top))))]
                (doseq [kind kinds]
                  (set! (.-strokeStyle ctx) (colors kind)) (set! (.-lineWidth ctx) 1.5)
                  (.beginPath ctx)
                  (doseq [[i [n m]] (map-indexed vector (filter (fn [[_ m]] (m kind)) history))
                          :let [v (nth (m kind) idx)]]
                    (if (zero? i) (.moveTo ctx (cx n) (cy v)) (.lineTo ctx (cx n) (cy v))))
                  (.stroke ctx))))
            (draw! []
              (let [{:keys [n current history]} @model
                    cw (.-width canvas) ch (.-height canvas)
                    side (min (/ (- cw 100) 4) (- ch 400))
                    ;; below the controls
                    gap 20 x0 (* 0.5 (- cw (* 4 side) (* 3 gap))) y0 110]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (doseq [[i kind] (map-indexed vector kinds)
                        :let [x (+ x0 (* i (+ side gap)))
                              {:keys [points] [[gx gy] d] :gap discrepancy :discrepancy} (current kind)
                              px (fn [[u v]] [(+ x (* side u)) (+ y0 (* side (- 1 v)))])]]
                  (set! (.-fillStyle ctx) "#161b25") (.fillRect ctx x y0 side side)
                  (text! kind x (- y0 8) (colors kind) 13)
                  (when points
                    (let [[qx qy] (px [gx gy])]
                      ;; the empty square, clipped to the panel
                      (.save ctx)
                      (.beginPath ctx) (.rect ctx x y0 side side) (.clip ctx)
                      (set! (.-fillStyle ctx) "rgba(255,255,255,0.08)")
                      (.fillRect ctx (- qx (* side d)) (- qy (* side d)) (* 2 side d) (* 2 side d))
                      (set! (.-strokeStyle ctx) "rgba(255,255,255,0.7)") (set! (.-lineWidth ctx) 1)
                      (.strokeRect ctx (- qx (* side d)) (- qy (* side d)) (* 2 side d) (* 2 side d))
                      (.restore ctx))
                    (set! (.-fillStyle ctx) (colors kind))
                    (doseq [p points :let [[a b] (px p)]]
                      (.beginPath ctx) (.arc ctx a b 2.2 0 (* 2 math/PI)) (.fill ctx))
                    (text! (str (count points) " points") x (+ y0 side 16) "#c9d2e0" 11)
                    (text! (str "dispersion " (.toFixed d 3) "  discrepancy " (.toFixed discrepancy 3)) x (+ y0 side 30) "#c9d2e0" 11)))
                (let [cy (+ y0 side 62) cw2 (/ (- (* 4 side) (* 3 gap) 20) 2) chh (- ch cy 70)]
                  (when (> chh 40)
                    (chart! x0 cy cw2 chh history 0 "dispersion against the count (log scale)" 0.5)
                    (chart! (+ x0 cw2 20 (* 3 gap)) cy cw2 chh history 1 "star discrepancy against the count (log scale)" 0.5)))
                (set! (.-textContent caption)
                      (str n " points  |  the white square: the largest empty one, half its side the dispersion"
                           "\nred random, gold Halton, green Hammersley, blue Sukharev (at the largest square count so far)"))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (advance!)
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))]
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.add gui controls "speed" 1 20 1)
        (.onChange (.add gui controls "most" 50 1000 10) restart!)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      (swap! model assoc :current (measures (:pool @model) 1))
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "sampling")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
