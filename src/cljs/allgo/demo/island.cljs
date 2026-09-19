(ns allgo.demo.island
  "Canvas demo for `allgo.procedural.island`.

  A polygon map, drawn as the polygons it is made of. Flat and
  two-dimensional on purpose: this is a map rather than a landscape, and
  the thing worth looking at is which cell is which, not what it would
  look like from a hillside.

  ## The views are the argument

  Each one is the output of one traversal, and switching between them on
  the same island is the point of the demo:

    `biome`     the finished map. Elevation and moisture through
                Whittaker's table.
    `water`     what the flood fill found -- ocean, lake, land. A lake
                differs from the sea by connectivity alone, and this is
                the view that shows it.
    `elevation` distance from the coast, with the curve applied. It is
                a dome, because that is what distance from the coast is.
    `moisture`  distance from fresh water. Dry interiors and green
                valleys, with nothing modelling either.
    `mesh`      the dual graph itself, with the Delaunay edges between
                cell centres over the Voronoi cells they belong to.
                Worth a look once: every other view is a property
                attached to this.

  ## A note on `points`

  The triangulation underneath is quadratic, so this is where the cost
  is, and `relax` multiplies it because each round is another
  triangulation. Five hundred points and one round is about a second and
  a half in a browser; the top of both sliders is several seconds. They
  are capped well below where it stops being a demo, and the ceiling is
  the triangulation rather than anything here -- see the note on cost in
  `allgo.geometry.dual-mesh`."
  (:require [allgo.geometry.dual-mesh :as dm]
            [allgo.procedural.island :as island]
            ["lil-gui" :default GUI]))

(def ^:private background "#0a0c16")
(def ^:private text-color "#7d8aa3")

(def ^:private ^js controls
  #js {:shape "noise"
       :points 500
       :relax 1
       :seed 1
       :view "biome"
       :rivers true
       :outlines true})

(defn- rgb [[r g b]] (str "rgb(" r "," g "," b ")"))

(defn- lerp3 [[r1 g1 b1] [r2 g2 b2] t]
  (let [t (max 0.0 (min 1.0 t))]
    [(js/Math.round (+ r1 (* t (- r2 r1))))
     (js/Math.round (+ g1 (* t (- g2 g1))))
     (js/Math.round (+ b1 (* t (- b2 b1))))]))

(defn- ramp
  "A colour from a list of stops, `t` in [0, 1]."
  [stops t]
  (let [n (dec (count stops))
        t (max 0.0 (min 0.9999 t))
        i (js/Math.floor (* t n))]
    (lerp3 (nth stops i) (nth stops (inc i)) (- (* t n) i))))

(def ^:private elevation-stops
  [[62 104 78] [122 150 86] [186 176 112] [150 120 84] [140 140 140] [255 255 255]])

(def ^:private moisture-stops
  [[196 176 124] [176 186 120] [120 168 110] [70 140 120] [44 104 140]])

(defn- cell-color [view c]
  (case view
    "biome" (rgb (get island/biome-colors (:biome c) [255 0 255]))
    "water" (cond (:ocean? c) "rgb(58,72,120)"
                  (:lake? c) "rgb(64,118,170)"
                  (:coast? c) "rgb(178,162,132)"
                  :else "rgb(96,132,86)")
    "elevation" (if (:ocean? c)
                  "rgb(44,54,92)"
                  (rgb (ramp elevation-stops (:elevation c))))
    "moisture" (if (:ocean? c)
                 "rgb(44,54,92)"
                 (rgb (ramp moisture-stops (:moisture c))))
    ;; The mesh view wants the cells muted, so the graph on top reads.
    (if (:water? c) "rgb(38,46,76)" "rgb(64,70,92)")))

(defn- draw! [^js ctx m [w h]]
  (set! (.-fillStyle ctx) background)
  (.fillRect ctx 0 0 w h)
  (when m
    (let [{:keys [centers corners edges bounds]} m
          [[x0 y0] [x1 y1]] bounds
          pad 6
          side (- (min w (- h 22)) (* 2 pad))
          sx (fn [x] (+ pad (* side (/ (- x x0) (- x1 x0)))))
          sy (fn [y] (+ pad (* side (/ (- y y0) (- y1 y0)))))
          view (.-view controls)
          mesh? (= view "mesh")
          ring! (fn [c]
                  (let [pts (dm/polygon m c)]
                    (.beginPath ctx)
                    (doseq [[i [x y]] (map-indexed vector pts)]
                      (if (zero? i) (.moveTo ctx (sx x) (sy y)) (.lineTo ctx (sx x) (sy y))))
                    (.closePath ctx)))]
      (doseq [c centers]
        (ring! c)
        (set! (.-fillStyle ctx) (cell-color view c))
        (.fill ctx)
        ;; Filling and stroking the same path closes the hairline seams
        ;; antialiasing leaves between neighbouring polygons.
        (when (or (.-outlines controls) (not mesh?))
          (set! (.-strokeStyle ctx)
                (if (.-outlines controls) "rgba(10,12,22,0.45)" (cell-color view c)))
          (set! (.-lineWidth ctx) 1)
          (.stroke ctx)))
      ;; The coastline, which is every edge with the sea on exactly one
      ;; side. A one-sided edge is the rim of the map and is not a coast.
      (when-not mesh?
        (set! (.-strokeStyle ctx) "rgba(12,16,30,0.85)")
        (set! (.-lineWidth ctx) 1.4)
        (.beginPath ctx)
        (doseq [e edges
                :let [[a b] (:centers e)]
                :when (and a b (not= (:ocean? (centers a)) (:ocean? (centers b))))]
          (let [[v0 v1] (:corners e)
                [ax ay] (:point (corners v0))
                [bx by] (:point (corners v1))]
            (.moveTo ctx (sx ax) (sy ay))
            (.lineTo ctx (sx bx) (sy by))))
        (.stroke ctx))
      (when mesh?
        ;; The Delaunay half of the dual: centre to centre across every
        ;; edge, which is the adjacency every pass actually walks.
        (set! (.-strokeStyle ctx) "rgba(150,180,230,0.35)")
        (set! (.-lineWidth ctx) 0.6)
        (.beginPath ctx)
        (doseq [e edges
                :let [[a b] (:centers e)]
                :when (and a b)]
          (let [[ax ay] (:point (centers a))
                [bx by] (:point (centers b))]
            (.moveTo ctx (sx ax) (sy ay))
            (.lineTo ctx (sx bx) (sy by))))
        (.stroke ctx))
      (when (.-rivers controls)
        (set! (.-strokeStyle ctx) "rgb(74,126,196)")
        (set! (.-lineCap ctx) "round")
        (doseq [e edges
                :let [flow (long (:river e 0))]
                :when (pos? flow)]
          (let [[v0 v1] (:corners e)
                [ax ay] (:point (corners v0))
                [bx by] (:point (corners v1))]
            (set! (.-lineWidth ctx) (min 4.0 (+ 0.9 (* 0.8 (js/Math.sqrt flow)))))
            (.beginPath ctx)
            (.moveTo ctx (sx ax) (sy ay))
            (.lineTo ctx (sx bx) (sy by))
            (.stroke ctx))))
      (let [land (remove :water? centers)
            lakes (count (filter :lake? centers))
            rivers (count (filter #(pos? (long (:river % 0))) edges))]
        (set! (.-fillStyle ctx) text-color)
        (set! (.-font ctx) "11px ui-monospace, SFMono-Regular, Menlo, monospace")
        (.fillText ctx
                   (str (count centers) " cells  ·  "
                        (js/Math.round (* 100.0 (/ (count land) (count centers)))) "% land  ·  "
                        lakes (if (= 1 lakes) " lake" " lakes") "  ·  "
                        rivers " river edges")
                   10 (- h 8))))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        state (atom {:bounds [0 0] :map nil})]
    (set! (.-style canvas) "display:block;border-radius:8px")
    (.appendChild container canvas)
    (letfn [(render! [] (draw! ctx (:map @state) (:bounds @state)))
            (resize! []
              (let [dpr (or js/window.devicePixelRatio 1)
                    w (.-clientWidth container)
                    h (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (set! (.-width canvas) (* w dpr))
                  (set! (.-height canvas) (* h dpr))
                  (set! (.-width (.-style canvas)) (str w "px"))
                  (set! (.-height (.-style canvas)) (str h "px"))
                  (.setTransform ctx dpr 0 0 dpr 0 0)
                  (swap! state assoc :bounds [w h])
                  (render!))))
            (generate! []
              (let [seed (long (.-seed controls))
                    ;; Seeded rather than `rand`, so that nudging a slider
                    ;; and putting it back gives the island you had.
                    stream (atom (+ 1 (* 2654435761 seed)))
                    rng (fn []
                          (let [x (swap! stream #(mod (* 16807 %) 2147483647))]
                            (/ (double x) 2147483647.0)))]
                (swap! state assoc :map
                       (island/generate {:points (long (.-points controls))
                                         :relax (long (.-relax controls))
                                         :shape (keyword (.-shape controls))
                                         :seed seed
                                         :rng rng})))
              (render!))]
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (generate!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "view" #js ["biome" "water" "elevation" "moisture" "mesh"])
              (.onChange render!))
          (-> (.add controls "shape" #js ["noise" "radial"]) (.onFinishChange generate!))
          (-> (.add controls "seed" 1 60 1) (.onFinishChange generate!))
          (-> (.add controls "points" 200 1200 100) (.onFinishChange generate!))
          (-> (.add controls "relax" 0 3 1) (.onFinishChange generate!))
          (-> (.add controls "rivers") (.onChange render!))
          (-> (.add controls "outlines") (.onChange render!))
          (.add #js {:regenerate generate!} "regenerate")))
      {:start (fn [] (resize!) (render!))
       :stop (fn [])})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "island")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
