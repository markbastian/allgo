(ns procedurals.terrain-shape
  (:require [procedurals.mesh :as mesh]
            [procedurals.terrain :as terrain]
            [procedurals.tin :as tin]
            [quil.applet :as applet]
            [quil.core :as q])
  (:import [processing.core PApplet PConstants PShape]))

(defn terrain-color [height-t slope-t noise-t]
  (let [[r g b] (mesh/terrain-rgb height-t slope-t noise-t)]
    (q/color r g b)))

(defn build-mesh-shape
  "Bakes the heightmap into a single retained PShape (GPU-buffered) so `draw`
  can render the whole terrain with one `q/shape` call instead of re-issuing
  tens of thousands of vertex/fill calls every frame."
  [grid dim lo hi cell-scale height-scale noise-scale]
  (let [offset     (* 0.5 (dec dim) cell-scale)
        slopes     (for [i (range dim) j (range dim)] (mesh/slope-at grid dim i j cell-scale height-scale))
        max-slope  (max 1e-6 (apply max slopes))
        ^PShape sh (.createShape ^PApplet (applet/current-applet))]
    (.beginShape sh PConstants/QUADS)
    (.noStroke sh)
    (doseq [i (range (dec dim))
            j (range (dec dim))]
      (doseq [[ii jj] [[i j] [i (inc j)] [(inc i) (inc j)] [(inc i) j]]]
        (let [h            (get-in grid [ii jj])
              height-t     (if (= lo hi) 0.5 (/ (double (- h lo)) (- hi lo)))
              slope-t      (/ (mesh/slope-at grid dim ii jj cell-scale height-scale) max-slope)
              wx           (- (* jj cell-scale) offset)
              wz           (- (* ii cell-scale) offset)
              noise-t      (mesh/noise-at wx wz noise-scale)
              [nx ny nz]   (mesh/normal-at grid dim ii jj cell-scale height-scale -1.0)]
          (.fill sh (unchecked-int (terrain-color height-t slope-t noise-t)))
          (.normal sh (float nx) (float ny) (float nz))
          (.vertex sh (float wx) (float (- (* h height-scale))) (float wz)))))
    (.endShape sh)
    sh))

(defn- diamond-square-grid [{:keys [iterations width]}]
  (let [grid-data (terrain/generate {:width width :iterations iterations
                                     :corners [0.0 (rand) (rand) (rand)]})]
    {:grid (terrain/cells->grid grid-data) :dim (:dim grid-data)}))

(defn- tin-grid [{:keys [tin-points dim smooth-passes] :or {tin-points 500 dim 129 smooth-passes 3}}]
  {:grid (-> (tin/generate {:n tin-points :size 1.0}) (tin/sample-grid dim) (tin/smooth-grid smooth-passes))
   :dim  dim})

(defn build-terrain
  "Builds a heightmap via `:generator` (`:diamond-square`, the default
  regular-grid midpoint-displacement fractal, or `:tin`, a Delaunay-
  triangulated irregular network per `procedurals.tin`) and bakes it into
  a PShape."
  [{:keys [generator height-scale cell-scale noise-scale] :or {generator :diamond-square noise-scale 180.0} :as config}]
  (let [{:keys [grid dim]} (case generator
                             :diamond-square (diamond-square-grid config)
                             :tin (tin-grid config))
        hs (flatten grid)
        lo (apply min hs)
        hi (apply max hs)]
    {:dim dim
     :terrain-shape (build-mesh-shape grid dim lo hi cell-scale height-scale noise-scale)}))
