(ns procedurals.boids-3d
  "Three.js demo for `procedurals.boids`, driven by the same rules as the 2D
  canvas viewer -- the flocking namespace is arity-generic, so going to three
  dimensions is purely a matter of passing 3-element bounds.

  Boids are instanced cones oriented along their heading and tinted by it
  (hue from compass bearing, brightness from climb), so a flock settling into
  alignment reads as the swarm converging on a single colour."
  (:require [procedurals.boids :as boids]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def ^:private world [100 100 100])
(def ^:private capacity 400)

(def ^:private controls
  #js {:boids         60
       :separation    1.6
       :alignment     1.0
       :cohesion      0.9
       :perception    22
       :personalSpace 9
       :speed         0.8})

(defn- params []
  {:separation-weight (.-separation controls)
   :alignment-weight  (.-alignment controls)
   :cohesion-weight   (.-cohesion controls)
   :perception-radius (.-perception controls)
   :separation-radius (.-personalSpace controls)
   :max-speed         (.-speed controls)
   :max-force         0.02})

(defn- boid-geometry
  "A cone nose-up rotated onto +Z, the axis `Object3D.lookAt` aligns."
  []
  (doto (THREE/ConeGeometry. 0.9 2.8 6)
    (.rotateX (/ js/Math.PI 2))))

(defn- bounds-box [[w h d]]
  (THREE/LineSegments.
   (THREE/EdgesGeometry. (THREE/BoxGeometry. w h d))
   (THREE/LineBasicMaterial. #js {:color 0x2c3a5c})))

(defn- resize-to-flock [flock n max-speed]
  (let [have (count flock)]
    (cond
      (= have n) flock
      (< n have) (subvec flock 0 n)
      :else      (into flock (repeatedly (- n have) #(boids/random-boid world max-speed))))))

(defn- write-instances!
  "Push flock state into the InstancedMesh: one matrix and colour per boid.
  Boid space is [0,size); the mesh is centred on the origin, hence the offset."
  [^js mesh ^js scratch ^js color flock [w h d]]
  (let [offset [(/ w 2) (/ h 2) (/ d 2)]]
    (doseq [[i {:keys [pos vel]}] (map-indexed vector flock)]
      (let [[x y z]    (boids/v- pos offset)
            [nx ny nz] (boids/normalize vel)]
        (.set (.-position scratch) x y z)
        (.lookAt scratch (+ x nx) (+ y ny) (+ z nz))
        (.updateMatrix scratch)
        (.setMatrixAt mesh i (.-matrix scratch))
        (.setHSL color
                 (+ 0.5 (/ (js/Math.atan2 nz nx) (* 2 js/Math.PI)))
                 0.8
                 (+ 0.55 (* 0.18 ny)))
        (.setColorAt mesh i color)))
    (set! (.-count mesh) (count flock))
    (set! (.-needsUpdate (.-instanceMatrix mesh)) true)
    (when-let [ic (.-instanceColor mesh)] (set! (.-needsUpdate ic) true))))

(defn init! [^js container]
  (let [scene    (THREE/Scene.)
        camera   (THREE/PerspectiveCamera. 55 (/ (.-clientWidth container) (.-clientHeight container)) 0.5 5000)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        mesh     (THREE/InstancedMesh. (boid-geometry)
                                       (THREE/MeshLambertMaterial.)
                                       capacity)
        scratch  (THREE/Object3D.)
        color    (THREE/Color.)
        state    (atom {:flock (boids/flock (.-boids controls) world (params))})]
    (set! (.-background scene) (THREE/Color. 0x05070d))
    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
    (.setSize renderer (.-clientWidth container) (.-clientHeight container))
    (set! (.-borderRadius (.-style (.-domElement renderer))) "8px")
    (.appendChild container (.-domElement renderer))
    (.add scene mesh)
    (.add scene (bounds-box world))
    (.add scene (THREE/AmbientLight. 0xffffff 0.65))
    (let [sun (THREE/DirectionalLight. 0xfff4e0 0.9)]
      (.set (.-position sun) 80 120 60)
      (.add scene sun))
    (.set (.-position camera) 115 75 115)
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (.set (.-target orbit) 0 0 0)
      (letfn [(reset-flock! []
                (swap! state assoc :flock (boids/flock (.-boids controls) world (params))))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (set! (.-aspect camera) (/ w h))
                  (.updateProjectionMatrix camera)
                  (.setSize renderer w h)))
              (animate []
                (js/requestAnimationFrame animate)
                (let [p     (params)
                      flock (-> (:flock @state)
                                (resize-to-flock (.-boids controls) (:max-speed p))
                                (boids/step world p))]
                  (swap! state assoc :flock flock)
                  (write-instances! mesh scratch color flock world))
                (.update orbit)
                (.render renderer scene camera))]
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (let [gui (GUI. #js {:container container})]
          (.add gui controls "boids" 10 capacity 10)
          (.add gui controls "separation" 0 3 0.1)
          (.add gui controls "alignment" 0 3 0.1)
          (.add gui controls "cohesion" 0 3 0.1)
          (.add gui controls "perception" 8 50 1)
          (.add gui controls "personalSpace" 2 30 1)
          (.add gui controls "speed" 0.2 2.5 0.1)
          (.add gui #js {:reset reset-flock!} "reset"))
        (animate)))))

(defn main []
  (when-let [container (js/document.getElementById "boids-3d")]
    (init! container)))

(main)
