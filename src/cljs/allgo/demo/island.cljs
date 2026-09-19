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
    `territory` who holds what, from `allgo.procedural.settlement`. Worth
                comparing against `elevation`: the borders sit on the
                ridges, because territory is assigned by travel cost and
                a ridge is expensive to cross from either side. Assigned
                by distance instead this would be a Voronoi diagram of
                the towns and would ignore the land entirely.
    `mesh`      the dual graph itself, with the Delaunay edges between
                cell centres over the Voronoi cells they belong to.
                Worth a look once: every other view is a property
                attached to this.

  Towns, roads, borders and names come from
  `allgo.procedural.settlement`. The names are worth reading rather than
  glancing at: each culture invents its own small phonology, so a run of
  neighbouring towns share a sound and the places where that changes are
  the cultural borders -- which do not line up with the political ones,
  because a culture spans several realms.

  `noisy` is the one to toggle back and forth. It redraws rather than
  regenerates -- the paths are always there -- so the same island appears
  as polygons and then as a coastline, and nothing else about it moves.
  Straight edges are what the model actually computed; the wandering ones
  are the same boundaries with the ink allowed to wander inside the two
  cells that share them.

  ## A note on `points`

  `relax` multiplies the cost, because each round is another
  triangulation of the whole point set. The triangulation itself used to
  be quadratic and was the reason this demo topped out at twelve hundred
  points; now that it indexes its triangles the sliders go further, and
  the cost rises with the point count rather than with its square."
  (:require [allgo.geometry.dual-mesh :as dm]
            [allgo.procedural.island :as island]
            [allgo.procedural.settlement :as settlement]
            ["lil-gui" :default GUI]))

(def ^:private background "#0a0c16")
(def ^:private text-color "#7d8aa3")

(def ^:private ^js controls
  #js {:shape "noise"
       :points 800
       :relax 2
       :seed 1
       :view "biome"
       :towns 14
       :labels true
       :rivers true
       :roads true
       :outlines true
       :noisy true})

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

(defn- territory-color
  "A hue per town, spread by the golden ratio so that neighbouring ids do
  not come out as neighbouring colours."
  [idx]
  (str "hsl(" (js/Math.round (* 360 (mod (* (inc idx) 0.61803398875) 1.0)))
       ",42%,62%)"))

(defn- cell-color [view c town-index]
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
    "territory" (cond
                  (:ocean? c) "rgb(44,54,92)"
                  (:lake? c) "rgb(64,118,170)"
                  (:territory c) (territory-color (get town-index (:territory c) 0))
                  ;; Land no town can walk to. An islet of its own.
                  :else "rgb(96,99,104)")
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
          town-index (into {} (map-indexed (fn [i t] [t i])) (:towns m))
          noisy? (.-noisy controls)
          ;; A render-time switch, not a regeneration. The paths are
          ;; always computed; this only decides whether to follow them,
          ;; which makes the comparison instant and makes the point that
          ;; nothing about the map changed.
          edge-path (fn [e]
                      (let [[v0 v1] (:corners e)]
                        (if (and noisy? (:path e))
                          (:path e)
                          [(:point (corners v0)) (:point (corners v1))])))
          stroke-path (fn [pts]
                        (.beginPath ctx)
                        (doseq [[i [x y]] (map-indexed vector pts)]
                          (if (zero? i) (.moveTo ctx (sx x) (sy y)) (.lineTo ctx (sx x) (sy y))))
                        (.stroke ctx))
          ring! (fn [c]
                  (let [pts (if noisy? (dm/noisy-polygon m c) (dm/polygon m c))]
                    (.beginPath ctx)
                    (doseq [[i [x y]] (map-indexed vector pts)]
                      (if (zero? i) (.moveTo ctx (sx x) (sy y)) (.lineTo ctx (sx x) (sy y))))
                    (.closePath ctx)))]
      (doseq [c centers]
        (ring! c)
        (set! (.-fillStyle ctx) (cell-color view c town-index))
        (.fill ctx)
        ;; Filling and stroking the same path closes the hairline seams
        ;; antialiasing leaves between neighbouring polygons.
        (when (or (.-outlines controls) (not mesh?))
          (set! (.-strokeStyle ctx)
                (if (.-outlines controls) "rgba(10,12,22,0.45)" (cell-color view c town-index)))
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
          (doseq [[i [x y]] (map-indexed vector (edge-path e))]
            (if (zero? i) (.moveTo ctx (sx x) (sy y)) (.lineTo ctx (sx x) (sy y)))))
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
        (set! (.-lineJoin ctx) "round")
        (doseq [e edges
                :let [flow (long (:river e 0))]
                :when (pos? flow)]
          (set! (.-lineWidth ctx) (min 4.0 (+ 0.9 (* 0.8 (js/Math.sqrt flow)))))
          (stroke-path (edge-path e))))
      ;; Territory borders, where two towns' claims meet.
      (when (= view "territory")
        (set! (.-strokeStyle ctx) "rgba(16,16,26,0.9)")
        (set! (.-lineWidth ctx) 1.8)
        (set! (.-lineJoin ctx) "round")
        (doseq [e edges
                :let [[a b] (:centers e)]
                :when (and a b
                           (:territory (centers a))
                           (:territory (centers b))
                           (not= (:territory (centers a)) (:territory (centers b))))]
          (stroke-path (edge-path e))))
      (when (.-roads controls)
        ;; Roads run centre to centre, because a road goes through a
        ;; place rather than along its boundary -- which is exactly the
        ;; other half of the dual from where the rivers are.
        (set! (.-strokeStyle ctx) "rgb(226,206,160)")
        (set! (.-lineCap ctx) "round")
        (doseq [e edges
                :let [traffic (long (:road e 0))
                      [a b] (:centers e)]
                :when (and (pos? traffic) a b)]
          (set! (.-lineWidth ctx) (min 3.0 (+ 0.9 (* 0.5 traffic))))
          (stroke-path [(:point (centers a)) (:point (centers b))])))
      (doseq [t (:towns m)]
        (let [[x y] (:point (centers t))]
          (set! (.-fillStyle ctx) "rgb(24,22,18)")
          (.beginPath ctx) (.arc ctx (sx x) (sy y) 4.5 0 (* 2 js/Math.PI)) (.fill ctx)
          (set! (.-fillStyle ctx) "rgb(250,240,212)")
          (.beginPath ctx) (.arc ctx (sx x) (sy y) 3.0 0 (* 2 js/Math.PI)) (.fill ctx)))
      (when (.-labels controls)
        (set! (.-textAlign ctx) "center")
        (set! (.-lineJoin ctx) "round")
        ;; The bigger rivers, named along their length -- which is what
        ;; having a body rather than a scattering of flowing edges buys.
        ;; Placed a third of the way up from the mouth, where a river is
        ;; wide and the label has somewhere to sit.
        (set! (.-font ctx) "italic 9px ui-serif, Georgia, serif")
        (doseq [r (take 5 (:rivers m))
                :let [path (:path r)]
                :when (and (:name r) (>= (count path) 3))]
          (let [v (nth path (quot (* 2 (count path)) 3))
                [x y] (:point (corners v))]
            (set! (.-strokeStyle ctx) "rgba(8,10,18,0.85)")
            (set! (.-lineWidth ctx) 3)
            (.strokeText ctx (:name r) (sx x) (sy y))
            (set! (.-fillStyle ctx) "rgb(150,196,238)")
            (.fillText ctx (:name r) (sx x) (sy y)))))
      (when (and (.-labels controls) (seq (:towns m)))
        (set! (.-font ctx) "10px ui-sans-serif, system-ui, sans-serif")
        (set! (.-textAlign ctx) "center")
        (set! (.-lineJoin ctx) "round")
        (doseq [t (:towns m)]
          (let [c (centers t)
                [x y] (:point c)]
            (when-let [nm (:name c)]
              ;; Drawn twice: a dark stroke under a light fill, so a name
              ;; stays readable over pale sand and dark forest alike
              ;; without a box behind it.
              (set! (.-strokeStyle ctx) "rgba(8,10,18,0.9)")
              (set! (.-lineWidth ctx) 3)
              (.strokeText ctx nm (sx x) (- (sy y) 7))
              (set! (.-fillStyle ctx) "rgb(246,240,226)")
              (.fillText ctx nm (sx x) (- (sy y) 7)))))
        (set! (.-textAlign ctx) "left"))
      (let [land (remove :water? centers)
            lakes (count (filter :lake? centers))
            rivers (count (filter #(pos? (long (:river % 0))) edges))]
        (set! (.-fillStyle ctx) text-color)
        (set! (.-font ctx) "11px ui-monospace, SFMono-Regular, Menlo, monospace")
        (.fillText ctx
                   (str (count centers) " cells  ·  "
                        (js/Math.round (* 100.0 (/ (count land) (count centers)))) "% land  ·  "
                        lakes (if (= 1 lakes) " lake" " lakes") "  ·  "
                        rivers " river edges  ·  "
                        (count (:towns m)) " towns")
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
                       (cond-> (island/generate {:points (long (.-points controls))
                                                 :relax (long (.-relax controls))
                                                 :shape (keyword (.-shape controls))
                                                 :seed seed
                                                 :rng rng})
                         (pos? (long (.-towns controls)))
                         (settlement/populate {:towns (long (.-towns controls))}))))
              (render!))]
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (generate!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "view"
                    #js ["biome" "water" "elevation" "moisture" "territory" "mesh"])
              (.onChange render!))
          (-> (.add controls "shape" #js ["noise" "radial"]) (.onFinishChange generate!))
          (-> (.add controls "seed" 1 60 1) (.onFinishChange generate!))
          (-> (.add controls "points" 200 3000 100) (.onFinishChange generate!))
          (-> (.add controls "relax" 0 3 1) (.onFinishChange generate!))
          (-> (.add controls "towns" 0 30 1) (.onFinishChange generate!))
          (-> (.add controls "noisy") (.onChange render!))
          (-> (.add controls "labels") (.onChange render!))
          (-> (.add controls "rivers") (.onChange render!))
          (-> (.add controls "roads") (.onChange render!))
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
