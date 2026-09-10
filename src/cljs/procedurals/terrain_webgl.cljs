(ns procedurals.terrain-webgl
  (:require [procedurals.mesh :as mesh]
            [procedurals.terrain :as terrain]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def config
  {:iterations   8
   :width        1.0
   :height-scale 200
   :cell-scale   8
   :noise-scale  180.0})

(defn- set-vec3! [^js arr idx x y z]
  (let [base (* idx 3)]
    (aset arr base x)
    (aset arr (+ base 1) y)
    (aset arr (+ base 2) z)))

(defn build-terrain-data [{:keys [iterations width]}]
  (let [grid-data (terrain/generate {:width width :iterations iterations
                                     :corners [0.0 (rand) (rand) (rand)]})
        hs        (terrain/heights grid-data)]
    {:grid (terrain/cells->grid grid-data)
     :dim  (:dim grid-data)
     :lo   (apply min hs)
     :hi   (apply max hs)}))

(defn build-geometry
  "Bakes the heightmap into an indexed THREE.BufferGeometry: one vertex per
  grid cell (shared across its adjacent triangles) so normals interpolate
  smoothly, matching `procedurals.terrain-shape/build-mesh-shape` but using
  the standard y-up world (+height, not the Quil renderer's -height)."
  [{:keys [grid dim lo hi]} {:keys [cell-scale height-scale noise-scale]}]
  (let [n          (* dim dim)
        offset     (* 0.5 (dec dim) cell-scale)
        positions  (js/Float32Array. (* n 3))
        normals    (js/Float32Array. (* n 3))
        colors     (js/Float32Array. (* n 3))
        slopes     (for [i (range dim) j (range dim)] (mesh/slope-at grid dim i j cell-scale height-scale))
        max-slope  (max 1e-6 (apply max slopes))]
    (doseq [i (range dim) j (range dim)]
      (let [idx        (+ (* i dim) j)
            h          (get-in grid [i j])
            height-t   (if (= lo hi) 0.5 (/ (double (- h lo)) (- hi lo)))
            slope-t    (/ (mesh/slope-at grid dim i j cell-scale height-scale) max-slope)
            wx         (- (* j cell-scale) offset)
            wz         (- (* i cell-scale) offset)
            noise-t    (mesh/noise-at wx wz noise-scale)
            [nx ny nz] (mesh/normal-at grid dim i j cell-scale height-scale 1.0)
            [r g b]    (mesh/terrain-rgb height-t slope-t noise-t)]
        (set-vec3! positions idx wx (* h height-scale) wz)
        (set-vec3! normals idx nx ny nz)
        (set-vec3! colors idx (/ r 255) (/ g 255) (/ b 255))))
    (let [quads (* (dec dim) (dec dim))
          index (js/Uint32Array. (* quads 6))]
      (doseq [i (range (dec dim)) j (range (dec dim))]
        (let [a    (+ (* i dim) j)
              b    (+ (* i dim) (inc j))
              c    (+ (* (inc i) dim) j)
              d    (+ (* (inc i) dim) (inc j))
              base (* 6 (+ (* i (dec dim)) j))]
          (aset index base a) (aset index (+ base 1) c) (aset index (+ base 2) b)
          (aset index (+ base 3) b) (aset index (+ base 4) c) (aset index (+ base 5) d)))
      (doto (THREE/BufferGeometry.)
        (.setAttribute "position" (THREE/BufferAttribute. positions 3))
        (.setAttribute "normal" (THREE/BufferAttribute. normals 3))
        (.setAttribute "color" (THREE/BufferAttribute. colors 3))
        (.setIndex (THREE/BufferAttribute. index 1))))))

(defn build-mesh [config]
  (let [terrain-data (build-terrain-data config)
        geometry     (build-geometry terrain-data config)
        material     (THREE/MeshStandardMaterial. #js {:vertexColors true
                                                       :side         THREE/DoubleSide
                                                       :roughness    0.9
                                                       :metalness    0.0})]
    {:span (* (:cell-scale config) (dec (:dim terrain-data)))
     :mesh (THREE/Mesh. geometry material)}))

(defn- dispose-mesh! [^js terrain-mesh]
  (.dispose (.-geometry terrain-mesh))
  (.dispose (.-material terrain-mesh)))

(defn init! [^js container config]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 55 (/ (.-clientWidth container) (.-clientHeight container)) 1 100000)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        {:keys [span mesh]} (build-mesh config)
        state    (atom {:terrain-mesh mesh})]
    (set! (.-background scene) (THREE/Color. 0x0f0f19))
    (.setSize renderer (.-clientWidth container) (.-clientHeight container))
    (.appendChild container (.-domElement renderer))
    (.add scene mesh)
    (.add scene (THREE/AmbientLight. 0xffffff 0.5))
    (let [sun (THREE/DirectionalLight. 0xfff4e0 1.2)]
      (.set (.-position sun) (* span 0.4) (* span 1.0) (* span 0.3))
      (.add scene sun))
    (.set (.-position camera) (* span 0.9) (* span 0.7) (* span 0.9))
    (let [controls (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping controls) true)
      (.set (.-target controls) 0 0 0)
      (letfn [(regenerate! []
                (dispose-mesh! (:terrain-mesh @state))
                (.remove scene (:terrain-mesh @state))
                (let [{:keys [mesh]} (build-mesh config)]
                  (.add scene mesh)
                  (swap! state assoc :terrain-mesh mesh)))
              (on-key-down [^js e]
                (when (= (.-key e) "r") (regenerate!)))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (set! (.-aspect camera) (/ w h))
                  (.updateProjectionMatrix camera)
                  (.setSize renderer w h)))
              (animate []
                (js/requestAnimationFrame animate)
                (.update controls)
                (.render renderer scene camera))]
        (.addEventListener js/window "keydown" on-key-down)
        (.addEventListener js/window "resize" on-resize)
        (animate)))))

(defn main []
  (when-let [container (js/document.getElementById "terrain-app")]
    (init! container config)))

(main)
