(ns procedurals.boids-voronoi
  "A Voronoi diagram of a flock, rebuilt every frame.

  The mesh is recomputed from scratch rather than maintained incrementally.
  Incremental schemes -- Lawson flip repair, or a kinetic data structure
  scheduling flips from known trajectories -- pay off when a small part of
  the structure changes per step. Here every site moves every frame, so local
  repair would touch the whole mesh anyway; and boids wrap toroidally, which
  breaks the small-displacement assumption flip repair relies on. A rebuild
  is immune to that by construction and costs ~1-3ms at these flock sizes."
  (:require [procedurals.boids :as boids]
            [procedurals.delaunay :as delaunay]
            ["lil-gui" :default GUI]))

(def ^:private background "#05070d")

(def ^:private controls
  #js {:boids         40
       :separation    1.6
       :alignment     1.0
       :cohesion      0.9
       :perception    70
       :personalSpace 26
       :speed         1.8
       :cells         true
       :mesh          false})

(defn- params []
  {:separation-weight (.-separation controls)
   :alignment-weight  (.-alignment controls)
   :cohesion-weight   (.-cohesion controls)
   :perception-radius (.-perception controls)
   :separation-radius (.-personalSpace controls)
   :max-speed         (.-speed controls)})

(defn- heading-hue [[vx vy]]
  (mod (+ 180 (/ (* 180 (js/Math.atan2 vy vx)) js/Math.PI)) 360))

(defn- resize-to-flock [flock n bounds max-speed]
  (let [have (count flock)]
    (cond
      (= have n) flock
      (< n have) (subvec flock 0 n)
      :else      (into flock (repeatedly (- n have) #(boids/random-boid bounds max-speed))))))

(defn- trace-polygon! [^js ctx [[x0 y0] & rest]]
  (.beginPath ctx)
  (.moveTo ctx x0 y0)
  (doseq [[x y] rest] (.lineTo ctx x y))
  (.closePath ctx))

(defn- draw-cells! [^js ctx flock cells]
  (doseq [{:keys [pos vel]} flock
          :let [cell (cells pos)]
          :when (>= (count cell) 3)]
    (let [hue (heading-hue vel)]
      (trace-polygon! ctx cell)
      (set! (.-fillStyle ctx) (str "hsl(" hue " 62% 19%)"))
      (.fill ctx)
      (set! (.-strokeStyle ctx) (str "hsl(" hue " 70% 52%)"))
      (set! (.-lineWidth ctx) 1)
      (.stroke ctx))))

(defn- draw-mesh! [^js ctx triangles]
  (set! (.-strokeStyle ctx) "rgba(150, 170, 230, 0.28)")
  (set! (.-lineWidth ctx) 1)
  (doseq [{:keys [points]} triangles]
    (trace-polygon! ctx points)
    (.stroke ctx)))

(defn- draw-boids! [^js ctx flock]
  (doseq [{[x y] :pos vel :vel} flock]
    (set! (.-fillStyle ctx) (str "hsl(" (heading-hue vel) " 95% 72%)"))
    (.save ctx)
    (.translate ctx x y)
    (.rotate ctx (js/Math.atan2 (second vel) (first vel)))
    (.beginPath ctx)
    (.moveTo ctx 6 0)
    (.lineTo ctx -4 3)
    (.lineTo ctx -2 0)
    (.lineTo ctx -4 -3)
    (.closePath ctx)
    (.fill ctx)
    (.restore ctx)))

(defn- draw! [^js ctx flock [w h]]
  (set! (.-fillStyle ctx) background)
  (.fillRect ctx 0 0 w h)
  (let [sites (mapv :pos flock)]
    (when (.-cells controls)
      (draw-cells! ctx flock (delaunay/bounded-cells sites [[0 0] [w h]])))
    (when (.-mesh controls)
      ;; Of the boids alone -- edges to the cell sentinels would fly off-screen.
      (draw-mesh! ctx (delaunay/triangulate sites)))
    (draw-boids! ctx flock)))

(defn init! [^js container]
  (let [canvas   (js/document.createElement "canvas")
        ctx      (.getContext canvas "2d")
        running? (atom false)
        state    (atom {:bounds [0 0] :flock []})]
    (set! (.-style canvas) "display:block;border-radius:8px")
    (.appendChild container canvas)
    (letfn [(resize! []
              (let [dpr (or js/window.devicePixelRatio 1)
                    w   (.-clientWidth container)
                    h   (.-clientHeight container)]
                (when (and (pos? w) (pos? h))
                  (set! (.-width canvas) (* w dpr))
                  (set! (.-height canvas) (* h dpr))
                  (set! (.-width (.-style canvas)) (str w "px"))
                  (set! (.-height (.-style canvas)) (str h "px"))
                  (.setTransform ctx dpr 0 0 dpr 0 0)
                  (swap! state assoc :bounds [w h]))))
            (reset-flock! []
              (swap! state assoc :flock (boids/flock (.-boids controls) (:bounds @state) (params))))
            (tick []
              (when @running?
                (js/requestAnimationFrame tick)
                (let [{:keys [bounds]} @state
                      p     (params)
                      flock (-> (:flock @state)
                                (resize-to-flock (.-boids controls) bounds (:max-speed p))
                                (boids/step bounds p))]
                  (swap! state assoc :flock flock)
                  (draw! ctx flock bounds))))]
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (reset-flock!)
      (let [gui (GUI. #js {:container container})]
        (.add gui controls "boids" 10 120 5)
        (.add gui controls "separation" 0 3 0.1)
        (.add gui controls "alignment" 0 3 0.1)
        (.add gui controls "cohesion" 0 3 0.1)
        (.add gui controls "perception" 20 140 5)
        (.add gui controls "personalSpace" 5 60 1)
        (.add gui controls "speed" 0.3 4 0.1)
        (.add gui controls "cells")
        (.add gui controls "mesh")
        (.add gui #js {:reset reset-flock!} "reset"))
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (tick)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "boids-voronoi")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
