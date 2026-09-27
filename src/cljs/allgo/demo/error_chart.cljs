(ns allgo.demo.error-chart
  "The strip chart of relative energy error the integrator demos share:
  log10 |dE/E| over the last couple of hundred samples, with a rule every
  three decades. A drift in energy is the integrator's error made visible,
  and the log scale is what lets a 1e-12 method and a 1e-4 one share it.")

(def ^:private w 210)
(def ^:private h 74)
(def ^:private lo -14.0)
(def ^:private hi -1.0)

(defn mount!
  "A chart canvas in `container`'s corner. Returns its 2D context."
  [^js container]
  (let [canvas (js/document.createElement "canvas")]
    (set! (.-className canvas) "numeric-chart")
    (set! (.-width canvas) w)
    (set! (.-height canvas) h)
    (.appendChild container canvas)
    (.getContext canvas "2d")))

(defn push
  "`samples` with relative error `rel` added, the oldest dropped once the
  chart is full."
  [samples rel]
  (let [v (js/Math.log10 (max 1e-16 rel))]
    (if (>= (count samples) w) (conj (subvec samples 1) v) (conj samples v))))

(defn draw! [^js ctx samples]
  (.clearRect ctx 0 0 w h)
  (set! (.-fillStyle ctx) "rgba(5,7,13,0.72)")
  (.fillRect ctx 0 0 w h)
  (set! (.-strokeStyle ctx) "rgba(255,255,255,0.10)")
  (set! (.-lineWidth ctx) 1)
  (doseq [decade (range (js/Math.ceil lo) hi 3)]
    (let [y (* h (- 1.0 (/ (- decade lo) (- hi lo))))]
      (.beginPath ctx) (.moveTo ctx 0 y) (.lineTo ctx w y) (.stroke ctx)))
  (when (> (count samples) 1)
    (set! (.-strokeStyle ctx) "#6c8cff")
    (set! (.-lineWidth ctx) 1.5)
    (.beginPath ctx)
    (doseq [[i v] (map-indexed vector samples)]
      (let [x (* w (/ i (max 1 (dec (count samples)))))
            y (* h (- 1.0 (/ (- v lo) (- hi lo))))]
        (if (zero? i) (.moveTo ctx x y) (.lineTo ctx x y))))
    (.stroke ctx))
  (set! (.-fillStyle ctx) "rgba(148,148,171,0.9)")
  (set! (.-font ctx) "9px ui-monospace, Menlo, monospace")
  (.fillText ctx "log |dE/E|" 6 11))
