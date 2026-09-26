(ns allgo.demo.broad-phase
  "Four ways to avoid comparing everything with everything, racing each
  other. Ten Minute Physics 23 and 24, against the hash from 11.

  Finding which of n objects touch is n(n-1)/2 comparisons done directly,
  which is a thousand objects costing half a million tests a frame. Every
  method here is a different way of throwing most of those out before
  doing them, and they should all arrive at the same answer -- the `pairs`
  readout is the same number whichever is selected, which is the point:
  they are interchangeable.

    brute force   the half-million. Kept as the thing to beat, and as the
                  check that the others are right.
    spatial hash  bin by position. Needs one cell size chosen up front,
                  which is why `sizes: mixed` hurts it: the cell has to
                  suit the largest object, and then the smallest share
                  cells with far more neighbors than they touch.
    sweep+prune   sort by lower edge along one axis and compare only what
                  overlaps on it. Does not care that sizes differ, and
                  barely cares that things moved, because the sorted
                  order is nearly right already.
    Morton BVH    sort by Morton code, split the sorted list down the
                  middle, and query the tree. The sort is what builds the
                  tree; there is no search for a good split anywhere.

  `show: morton order` draws the Z-curve the codes put the objects in --
  the order the BVH is built from, and the reason a run of it makes a
  compact box.

  `sizes: mixed` is the interesting switch. Watch the hash's frame time
  against sweep's as it goes over."
  (:require [allgo.demo.fps :as fps]
            [allgo.spatial.bvh :as bvh]
            [allgo.spatial.hash :as hash]
            [allgo.spatial.morton :as morton]
            [allgo.spatial.sweep :as sweep]
            ["lil-gui" :default GUI]))

(def ^:private ^js controls
  #js {:objects 900
       :sizes   "equal"
       :method  "sweep and prune"
       :speed   1.0
       :show    "collisions"})

(defn- build [[w h]]
  (let [n     (.-objects controls)
        mixed? (= "mixed" (.-sizes controls))
        world-h 10.0
        world-w (* world-h (/ w (max 1 h)))
        pos   (js/Float64Array. (* 3 n))
        vel   (js/Float64Array. (* 3 n))
        rad   (js/Float64Array. n)
        half  (js/Float64Array. (* 3 n))]
    (dotimes [i n]
      (let [r (if mixed?
                ;; A few large, many small -- the distribution a single
                ;; cell size cannot serve.
                (+ 0.04 (* 0.55 (js/Math.pow (js/Math.random) 6)))
                0.11)]
        (aset rad i r)
        (aset pos (* 3 i) (+ r (* (- world-w (* 2 r)) (js/Math.random))))
        (aset pos (+ (* 3 i) 1) (+ r (* (- world-h (* 2 r)) (js/Math.random))))
        (aset pos (+ (* 3 i) 2) 0.0)
        (aset vel (* 3 i) (- (js/Math.random) 0.5))
        (aset vel (+ (* 3 i) 1) (- (js/Math.random) 0.5))
        (dotimes [k 3] (aset half (+ (* 3 i) k) r))
        ;; Flat in z, so every box overlaps on z and the 3D structures
        ;; answer a 2D question without knowing it.
        (aset half (+ (* 3 i) 2) 0.0)))
    {:n n :pos pos :vel vel :rad rad :half half
     :world [world-w world-h]
     :mins (js/Float64Array. (* 3 n))
     :maxs (js/Float64Array. (* 3 n))
     :max-radius (reduce max 0.0 (map #(aget rad %) (range n)))
     :hash (hash/spatial-hash (* 2 (reduce max 0.11 (map #(aget rad %) (range n)))) n)
     :sweep (sweep/sweep n)}))

(defn- integrate! [{:keys [n ^js pos ^js vel ^js rad world]} dt]
  (let [[w h] world]
    (dotimes [i n]
      (let [b (* 3 i)
            r (aget rad i)]
        (aset pos b (+ (aget pos b) (* dt (aget vel b))))
        (aset pos (+ b 1) (+ (aget pos (+ b 1)) (* dt (aget vel (+ b 1)))))
        (when (< (aget pos b) r)
          (aset pos b r) (aset vel b (- (aget vel b))))
        (when (> (aget pos b) (- w r))
          (aset pos b (- w r)) (aset vel b (- (aget vel b))))
        (when (< (aget pos (+ b 1)) r)
          (aset pos (+ b 1) r) (aset vel (+ b 1) (- (aget vel (+ b 1)))))
        (when (> (aget pos (+ b 1)) (- h r))
          (aset pos (+ b 1) (- h r)) (aset vel (+ b 1) (- (aget vel (+ b 1)))))))))

(defn- boxes! [{:keys [n ^js pos ^js half ^js mins ^js maxs]}]
  (dotimes [i n]
    (let [b (* 3 i)]
      (dotimes [k 3]
        (aset mins (+ b k) (- (aget pos (+ b k)) (aget half (+ b k))))
        (aset maxs (+ b k) (+ (aget pos (+ b k)) (aget half (+ b k))))))))

(defn- touching?
  "The narrow phase every method's candidates are put through, so that all
  four end at the same answer."
  [{:keys [^js pos ^js rad]} i j]
  (let [a (* 3 i) b (* 3 j)
        dx (- (aget pos a) (aget pos b))
        dy (- (aget pos (+ a 1)) (aget pos (+ b 1)))
        r (+ (aget rad i) (aget rad j))]
    (< (+ (* dx dx) (* dy dy)) (* r r))))

(defn- candidates
  "The pairs the chosen method thinks are worth testing."
  [{:keys [n pos mins maxs hash sweep max-radius]}]
  (case (.-method controls)
    "brute force" (for [i (range n) j (range (inc i) n)] [i j])
    "spatial hash" (hash/overlapping-pairs hash pos n (* 2 max-radius))
    "sweep and prune" (sweep/overlapping-pairs sweep mins maxs n)
    "morton BVH" (bvh/overlapping-pairs (bvh/build mins maxs n) mins maxs n)
    (sweep/overlapping-pairs sweep mins maxs n)))

(defn- resolve!
  "Separates a touching pair and swaps the velocity along the line between
  them. Equal masses, so they simply exchange."
  [{:keys [^js pos ^js vel ^js rad]} i j]
  (let [a (* 3 i) b (* 3 j)
        dx (- (aget pos b) (aget pos a))
        dy (- (aget pos (+ b 1)) (aget pos (+ a 1)))
        d (js/Math.hypot dx dy)]
    (when (pos? d)
      (let [nx (/ dx d) ny (/ dy d)
            overlap (- (+ (aget rad i) (aget rad j)) d)
            cx (* 0.5 overlap nx) cy (* 0.5 overlap ny)]
        (aset pos a (- (aget pos a) cx))
        (aset pos (+ a 1) (- (aget pos (+ a 1)) cy))
        (aset pos b (+ (aget pos b) cx))
        (aset pos (+ b 1) (+ (aget pos (+ b 1)) cy))
        (let [v1 (+ (* (aget vel a) nx) (* (aget vel (+ a 1)) ny))
              v2 (+ (* (aget vel b) nx) (* (aget vel (+ b 1)) ny))]
          (when (< (- v2 v1) 0)
            (let [dv (- v2 v1)]
              (aset vel a (+ (aget vel a) (* nx dv)))
              (aset vel (+ a 1) (+ (aget vel (+ a 1)) (* ny dv)))
              (aset vel b (- (aget vel b) (* nx dv)))
              (aset vel (+ b 1) (- (aget vel (+ b 1)) (* ny dv))))))))))

(defn- draw! [^js ctx {:keys [n ^js pos ^js rad world]} hits [w h] stats]
  (let [[_ world-h] world
        scale (/ h world-h)]
    (set! (.-fillStyle ctx) "#05070d")
    (.fillRect ctx 0 0 w h)
    (when (= "morton order" (.-show controls))
      ;; The Z curve the BVH is built from.
      (let [order (morton/sorted-points pos n)]
        (set! (.-strokeStyle ctx) "rgba(110,190,255,0.5)")
        (set! (.-lineWidth ctx) 1)
        (.beginPath ctx)
        (dotimes [k n]
          (let [b (* 3 (aget order k))
                x (* scale (aget pos b))
                y (- h (* scale (aget pos (+ b 1))))]
            (if (zero? k) (.moveTo ctx x y) (.lineTo ctx x y))))
        (.stroke ctx)))
    (dotimes [i n]
      (let [b (* 3 i)]
        (set! (.-fillStyle ctx) (if (hits i) "#f2b134" "#3b6ea5"))
        (.beginPath ctx)
        (.arc ctx (* scale (aget pos b)) (- h (* scale (aget pos (+ b 1))))
              (max 1.0 (* scale (aget rad i))) 0 (* 2 js/Math.PI))
        (.fill ctx)))
    (set! (.-fillStyle ctx) "rgba(0,0,0,0.55)")
    (.fillRect ctx 8 8 250 66)
    (set! (.-fillStyle ctx) "#dfe7f5")
    (set! (.-font ctx) "12px ui-monospace, monospace")
    (.fillText ctx (str "method    " (.-method controls)) 18 28)
    (.fillText ctx (str "candidates " (:candidates stats)) 18 45)
    (.fillText ctx (str "pairs      " (:pairs stats)
                        "   broad " (.toFixed (:ms stats) 2) " ms") 18 62)))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx    (.getContext canvas "2d")
        running? (atom false)
        tick-fps! (fps/meter! container)
        state  (atom {:bounds [0 0]})]
    (set! (.-style canvas) "display:block;border-radius:8px")
    (.appendChild container canvas)
    (letfn [(rebuild! []
              (let [{:keys [bounds]} @state]
                (when (pos? (first bounds))
                  (swap! state merge (build bounds)))))
            (resize! []
              (let [w (.-clientWidth container) h (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (set! (.-width canvas) w)
                  (set! (.-height canvas) h)
                  (set! (.-width (.-style canvas)) (str w "px"))
                  (set! (.-height (.-style canvas)) (str h "px"))
                  (swap! state assoc :bounds [w h])
                  (rebuild!))))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (let [t0 (js/performance.now)
                      st @state]
                  (integrate! st (* 0.016 (.-speed controls)))
                  (boxes! st)
                  (let [t1 (js/performance.now)
                        cands (candidates st)
                        ms (- (js/performance.now) t1)
                        hits (volatile! #{})
                        pairs (volatile! 0)]
                    (doseq [[i j] cands]
                      (when (touching? st i j)
                        (vswap! pairs inc)
                        (vswap! hits conj i j)
                        (resolve! st i j)))
                    (draw! ctx st @hits (:bounds st)
                           {:candidates (count cands) :pairs @pairs :ms ms}))
                  (tick-fps! (- (js/performance.now) t0)))))]
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "objects" 100 2500 100) (.onFinishChange rebuild!))
          (-> (.add controls "sizes" #js ["equal" "mixed"]) (.onChange rebuild!))
          (.add controls "method" #js ["brute force" "spatial hash"
                                       "sweep and prune" "morton BVH"])
          (.add controls "speed" 0 3 0.1)
          (.add controls "show" #js ["collisions" "morton order"])
          (.add #js {:reset rebuild!} "reset")))
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "broad-phase")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
