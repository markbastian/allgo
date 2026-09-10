(ns procedurals.boids-voronoi-3d
  "The 3D Voronoi diagram of a flock, rebuilt every frame.

  Cells are drawn as wireframe edges rather than solid polyhedra: filled
  cells occlude each other completely, so a solid rendering shows only the
  outer shell of the box. Faces are available translucently as an option,
  with depth writes off so the interior still reads."
  (:require [procedurals.boids :as boids]
            [procedurals.fps :as fps]
            [procedurals.voronoi3d :as v3]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private world [100 100 100])
(def ^:private box [[0 0 0] [100 100 100]])
(def ^:private max-verts 60000)

(def ^:private controls
  #js {:boids         28
       :separation    1.6
       :alignment     1.0
       :cohesion      0.9
       :perception    30
       :personalSpace 12
       :speed         0.5
       :cells         true
       :faces         false})

(defn- params []
  {:separation-weight (.-separation controls)
   :alignment-weight  (.-alignment controls)
   :cohesion-weight   (.-cohesion controls)
   :perception-radius (.-perception controls)
   :separation-radius (.-personalSpace controls)
   :max-speed         (.-speed controls)
   :max-force         0.02})

(defn- resize-to-flock [flock n max-speed]
  (let [have (count flock)]
    (cond
      (= have n) flock
      (< n have) (subvec flock 0 n)
      :else      (into flock (repeatedly (- n have) #(boids/random-boid world max-speed))))))

(defn- heading-rgb
  "Hue from compass bearing, brightness from climb -- the tint used by the
  other flocking demos, so a flock in alignment reads as one colour."
  [^js color [vx vy vz]]
  (let [m (js/Math.sqrt (+ (* vx vx) (* vy vy) (* vz vz)))
        n (if (zero? m) 1 m)]
    (.setHSL color
             (+ 0.5 (/ (js/Math.atan2 (/ vz n) (/ vx n)) (* 2 js/Math.PI)))
             0.8
             (+ 0.58 (* 0.16 (/ vy n))))))

(defn- write-edges!
  "Pack every cell's face loops into `positions`/`colors` as line segments,
  returning the vertex count written."
  [^js positions ^js colors diagram flock ^js color]
  (let [[ox oy oz] (mapv #(/ % 2.0) world)]
    (loop [[boid & more] flock, i 0]
      (if (or (nil? boid) (> i (- max-verts 64)))
        i
        (let [{:keys [faces]} (get diagram (:pos boid))]
          (heading-rgb color (:vel boid))
          (recur more
                 (reduce
                  (fn [i face]
                    (let [n (count face)]
                      (reduce
                       (fn [i k]
                         (let [[ax ay az] (nth face k)
                               [bx by bz] (nth face (rem (inc k) n))]
                           (doseq [[j [x y z]] [[i [ax ay az]] [(inc i) [bx by bz]]]]
                             (aset positions (* j 3) (- x ox))
                             (aset positions (+ (* j 3) 1) (- y oy))
                             (aset positions (+ (* j 3) 2) (- z oz))
                             (aset colors (* j 3) (.-r color))
                             (aset colors (+ (* j 3) 1) (.-g color))
                             (aset colors (+ (* j 3) 2) (.-b color)))
                           (+ i 2)))
                       i (range n))))
                  i faces)))))))

(defn- write-faces!
  "Fan-triangulate every face into `positions`/`colors`."
  [^js positions ^js colors diagram flock ^js color]
  (let [[ox oy oz] (mapv #(/ % 2.0) world)]
    (loop [[boid & more] flock, i 0]
      (if (or (nil? boid) (> i (- max-verts 64)))
        i
        (let [{:keys [faces]} (get diagram (:pos boid))]
          (heading-rgb color (:vel boid))
          (recur more
                 (reduce
                  (fn [i face]
                    (reduce
                     (fn [i k]
                       (doseq [[j [x y z]] [[i (nth face 0)]
                                            [(inc i) (nth face k)]
                                            [(+ i 2) (nth face (inc k))]]]
                         (aset positions (* j 3) (- x ox))
                         (aset positions (+ (* j 3) 1) (- y oy))
                         (aset positions (+ (* j 3) 2) (- z oz))
                         (aset colors (* j 3) (.-r color))
                         (aset colors (+ (* j 3) 1) (.-g color))
                         (aset colors (+ (* j 3) 2) (.-b color)))
                       (+ i 3))
                     i (range 1 (dec (count face)))))
                  i faces)))))))

(defn- write-sites!
  "One vertex per boid, so the sites themselves stay visible inside the mesh."
  [^js positions ^js colors flock ^js color]
  (let [[ox oy oz] (mapv #(/ % 2.0) world)]
    (reduce (fn [i {[x y z] :pos vel :vel}]
              (heading-rgb color vel)
              (aset positions (* i 3) (- x ox))
              (aset positions (+ (* i 3) 1) (- y oy))
              (aset positions (+ (* i 3) 2) (- z oz))
              (aset colors (* i 3) (.-r color))
              (aset colors (+ (* i 3) 1) (.-g color))
              (aset colors (+ (* i 3) 2) (.-b color))
              (inc i))
            0 flock)))

(defn- buffers []
  (let [geo (THREE/BufferGeometry.)]
    (doto geo
      (.setAttribute "position" (THREE/BufferAttribute. (js/Float32Array. (* max-verts 3)) 3))
      (.setAttribute "color" (THREE/BufferAttribute. (js/Float32Array. (* max-verts 3)) 3)))))

(defn- flush-geometry!
  "Publish the first `n` vertices. No bounding sphere: it would be computed
  over the whole buffer including stale slots, so these objects opt out of
  frustum culling instead."
  [^js geo n]
  (.setDrawRange geo 0 n)
  (set! (.-needsUpdate (.getAttribute geo "position")) true)
  (set! (.-needsUpdate (.getAttribute geo "color")) true))

(defn init! [^js container]
  (let [scene     (THREE/Scene.)
        camera    (THREE/PerspectiveCamera. 55 (/ (.-clientWidth container) (.-clientHeight container)) 0.5 5000)
        renderer  (THREE/WebGLRenderer. #js {:antialias true})
        edge-geo  (buffers)
        face-geo  (buffers)
        lines     (THREE/LineSegments. edge-geo (THREE/LineBasicMaterial.
                                                 #js {:vertexColors true :transparent true :opacity 0.55}))
        shells    (THREE/Mesh. face-geo (THREE/MeshBasicMaterial.
                                         #js {:vertexColors true :transparent true :opacity 0.09
                                              :depthWrite false :side THREE/DoubleSide}))
        site-geo  (buffers)
        dots      (THREE/Points. site-geo (THREE/PointsMaterial.
                                           #js {:vertexColors true :size 2.4 :sizeAttenuation true}))
        color     (THREE/Color.)
        running?  (atom false)
        tick-fps! (fps/meter! container)
        state     (atom {:flock (boids/flock (.-boids controls) world (params))})]
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
    (.setSize renderer (.-clientWidth container) (.-clientHeight container))
    (set! (.-borderRadius (.-style (.-domElement renderer))) "8px")
    (.appendChild container (.-domElement renderer))
    (doseq [^js obj [lines shells dots]]
      (set! (.-frustumCulled obj) false)
      (.add scene obj))
    (.add scene (THREE/LineSegments.
                 (THREE/EdgesGeometry. (let [[w h d] world] (THREE/BoxGeometry. w h d)))
                 (THREE/LineBasicMaterial. #js {:color 0x2c3a5c})))
    (.set (.-position camera) 118 78 118)
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 0 0)
      (letfn [(reset-flock! []
                (swap! state assoc :flock (boids/flock (.-boids controls) world (params))))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0    (js/performance.now)
                        p     (params)
                        flock (-> (:flock @state)
                                  (resize-to-flock (.-boids controls) (:max-speed p))
                                  (boids/step world p))
                        diagram (v3/diagram (mapv :pos flock) box)]
                    (swap! state assoc :flock flock)
                    (set! (.-visible lines) (.-cells controls))
                    (set! (.-visible shells) (.-faces controls))
                    (when (.-cells controls)
                      (flush-geometry!
                       edge-geo (write-edges! (.-array (.getAttribute edge-geo "position"))
                                              (.-array (.getAttribute edge-geo "color"))
                                              diagram flock color)))
                    (when (.-faces controls)
                      (flush-geometry!
                       face-geo (write-faces! (.-array (.getAttribute face-geo "position"))
                                              (.-array (.getAttribute face-geo "color"))
                                              diagram flock color)))
                    (flush-geometry!
                     site-geo (write-sites! (.-array (.getAttribute site-geo "position"))
                                            (.-array (.getAttribute site-geo "color"))
                                            flock color))
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (let [gui (GUI. #js {:container container})]
          (.add gui controls "boids" 8 60 2)
          (.add gui controls "separation" 0 3 0.1)
          (.add gui controls "alignment" 0 3 0.1)
          (.add gui controls "cohesion" 0 3 0.1)
          (.add gui controls "perception" 8 60 1)
          (.add gui controls "personalSpace" 2 30 1)
          (.add gui controls "speed" 0.2 2 0.1)
          (.add gui controls "cells")
          (.add gui controls "faces")
          (.add gui #js {:reset reset-flock!} "reset"))
        {:start (fn [] (when-not @running? (reset! running? true) (on-resize) (animate)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "boids-voronoi-3d")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
