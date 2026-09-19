(ns allgo.demo.dla
  "Canvas demo for `allgo.procedural.dla`.

  Particles wander in from the edge and stick where they first touch what
  is already there. Watching it is the point: nothing in the rule
  mentions branching, and a branching thing grows anyway. It happens
  because the cluster shadows itself -- a walker is far likelier to meet
  a tip sticking out than to find its way into a gap between two -- so
  tips outrun gaps and every tip becomes a branch with tips of its own.

  The whole cluster is grown once and then revealed in the order the
  particles arrived, rather than being grown a frame at a time. It comes
  to the same picture and it means `speed` scrubs an existing sequence
  instead of reseeding the aggregate every time it is dragged.

  ## The two views

  `skeleton` is the tree, drawn from each cell to the one it stuck to and
  thickened by how much hangs off it. That is what makes this terrain
  rather than a texture: the cluster is not a set of cells but a tree,
  and a tree can be measured. A cell carrying half the cluster is a trunk
  and should be high; a cell carrying nothing is a twig.

  `terrain` is that measurement shaded as ground -- the same heights
  `allgo.procedural.dla/heightmap` hands to the rest of the repository,
  blurred and refined so a one-cell skeleton has flanks. Ridges branch
  the way real ones do, with a direction out from the trunk, which is the
  thing `allgo.procedural.fractal` cannot produce at any setting."
  (:require [allgo.procedural.dla :as dla]
            ["lil-gui" :default GUI]))

(def ^:private background "#07080f")
(def ^:private text-color "#7d8aa3")

(def ^:private ^js controls
  #js {:dim 129
       :density 6
       :seed 1
       :view "skeleton"
       :speed 60
       :grow true})

(defn- hillshade-color [n dzdx dzdy]
  (let [l (max 0.15 (min 1.0 (/ (+ 1.0 (* -0.6 dzdx) (* -0.6 dzdy)) 2.2)))
        v (js/Math.round (* 255 (min 1.0 (* (+ 0.28 (* 0.72 n)) l))))]
    (str "rgb(" v "," (js/Math.round (* 0.98 v)) "," (js/Math.round (* 0.93 v)) ")")))

(defn- draw-terrain! [^js ctx {:keys [^js heights dim]} pad side]
  (let [dim (long dim)
        n (* dim dim)
        [lo hi] (loop [i 0 lo ##Inf hi ##-Inf]
                  (if (= i n) [lo hi]
                      (let [v (aget heights i)] (recur (inc i) (min lo v) (max hi v)))))
        span (max 1e-9 (- hi lo))
        at (fn [i j] (aget heights (+ (* (min (dec dim) (max 0 i)) dim)
                                      (min (dec dim) (max 0 j)))))
        px (/ (double side) dim)
        step (js/Math.ceil px)]
    (dotimes [i dim]
      (dotimes [j dim]
        (let [v (/ (- (at i j) lo) span)
              dzdx (* (- (at i (inc j)) (at i (dec j))) (/ 1.0 span) 30.0)
              dzdy (* (- (at (inc i) j) (at (dec i) j)) (/ 1.0 span) 30.0)]
          (set! (.-fillStyle ctx) (hillshade-color v dzdx dzdy))
          (.fillRect ctx (+ pad (* j px)) (+ pad (* i px)) step step))))))

(defn- draw-skeleton! [^js ctx cluster ^js sizes shown pad side]
  (let [dim (long (:dim cluster))
        ^js parent (:parent cluster)
        seed (:seed cluster)
        total (aget sizes seed)
        top (js/Math.log (+ 1.0 total))
        sc (/ (double side) dim)
        px (fn [c] (+ pad (* sc (rem c dim))))
        py (fn [c] (+ pad (* sc (quot c dim))))]
    (set! (.-lineCap ctx) "round")
    (doseq [cell shown
            :when (not= cell seed)]
      (let [p (aget parent cell)
            t (/ (js/Math.log (+ 1.0 (aget sizes cell))) top)
            v (js/Math.round (+ 70 (* 185 t)))]
        (set! (.-strokeStyle ctx)
              (str "rgb(" v "," (js/Math.round (* 0.9 v)) "," (js/Math.round (* 0.72 v)) ")"))
        (set! (.-lineWidth ctx) (+ 0.8 (* 3.2 t)))
        (.beginPath ctx)
        (.moveTo ctx (px cell) (py cell))
        (.lineTo ctx (px p) (py p))
        (.stroke ctx)))))

(defn- draw! [^js ctx state [w h]]
  (set! (.-fillStyle ctx) background)
  (.fillRect ctx 0 0 w h)
  (let [{:keys [cluster sizes terrain revealed]} state]
    (when cluster
      (let [pad 6
            side (- (min w (- h 22)) (* 2 pad))
            order (:order cluster)
            shown (take revealed order)
            done? (>= revealed (count order))]
        (if (and (= (.-view controls) "terrain") done? terrain)
          (draw-terrain! ctx terrain pad side)
          (draw-skeleton! ctx cluster sizes shown pad side))
        (set! (.-fillStyle ctx) text-color)
        (set! (.-font ctx) "11px ui-monospace, SFMono-Regular, Menlo, monospace")
        (.fillText ctx
                   (str (min revealed (count order)) " / " (count order) " particles"
                        (when (and (= (.-view controls) "terrain") (not done?))
                          "  ·  terrain once it finishes"))
                   10 (- h 8))))))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        state (atom {:bounds [0 0] :cluster nil :revealed 0})
        running? (atom false)]
    (set! (.-style canvas) "display:block;border-radius:8px")
    (.appendChild container canvas)
    (letfn [(render! [] (draw! ctx @state (:bounds @state)))
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
                    ;; Park-Miller, so a seed means the same cluster here
                    ;; as on the JVM.
                    stream (atom (+ 1 (* 2654435761 seed)))
                    rng (fn [] (/ (double (swap! stream #(mod (* 16807 %) 2147483647)))
                                  2147483647.0))
                    cluster (dla/aggregate {:dim (long (.-dim controls))
                                            :density (/ (double (.-density controls)) 100.0)
                                            :rng rng})]
                (swap! state assoc
                       :cluster cluster
                       :sizes (dla/subtree-sizes cluster)
                       ;; Refined once: enough to give the skeleton sides
                       ;; without quadrupling the pixels to shade.
                       :terrain (dla/heightmap cluster {:refine 1 :blur 2})
                       :revealed (if (.-grow controls) 1 ##Inf)))
              (render!))
            (animate []
              (when @running?
                (js/requestAnimationFrame animate)
                (let [{:keys [cluster revealed]} @state]
                  (when (and cluster (< revealed (count (:order cluster))))
                    (swap! state update :revealed + (long (.-speed controls)))
                    (render!)))))]
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (generate!)
      (let [gui (GUI. #js {:container container})]
        (doto gui
          (-> (.add controls "view" #js ["skeleton" "terrain"]) (.onChange render!))
          (-> (.add controls "dim" #js [65 97 129 193 257]) (.onFinishChange generate!))
          (-> (.add controls "density" 1 15 1) (.onFinishChange generate!))
          (-> (.add controls "seed" 1 60 1) (.onFinishChange generate!))
          (-> (.add controls "speed" 5 400 5))
          (-> (.add controls "grow") (.onChange generate!))
          (.add #js {:regrow generate!} "regrow")))
      {:start (fn [] (when-not @running?
                       (reset! running? true)
                       (resize!)
                       (animate)))
       :stop (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "dla")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
