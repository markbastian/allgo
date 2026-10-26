(ns allgo.demo.city
  "Three.js demo for `allgo.procedural.city`: the whole of Parish and
  Mueller's pipeline, in three dimensions.

  Land, roads, blocks, lots, a building on every lot, and facades on the
  buildings -- each stage the input of the next, from nothing but a
  seed. The buildings are grown by their L-system from bounding boxes,
  and `detail` is how many rounds of it to run: at 0 every building is
  its box, the paper's coarsest level of detail, and each round refines
  every box into what is inside it -- the base and shaft of a tower, its
  setbacks, its antenna; an office's sections; a house and its roof.
  What is still to be rewritten is drawn plain gray.

  The facades are the layered grids of `allgo.procedural.city.facade`,
  evaluated into a texture per style -- a stretch of upper floors that
  repeats both ways, and the ground floor, which only repeats along the
  street, because it is a different layer. Drag to orbit, scroll to
  zoom, right-drag to pan."
  (:require [allgo.demo.fps :as fps]
            [allgo.geometry.polygon :as poly]
            [allgo.procedural.city.buildings :as buildings]
            [allgo.procedural.city.facade :as facade]
            [allgo.procedural.city.land :as land]
            [allgo.procedural.city.lots :as lots]
            [allgo.procedural.city.roads :as roads]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]
            [clojure.math :as math]))

(def ^:private ^js controls
  #js {:seed     2
       :pattern  "mixed"
       :detail   8
       :facades  true})

;; ---------------------------------------------------------------------------
;; Facade textures

(def ^:private tiles
  "For each style: how wide a tile is (a whole number of its bays, so it
  repeats without a seam), its floor height, and how many floors the
  upper tile holds."
  {:tower  {:w 9.6 :floor 3.8 :floors 4}
   :base   {:w 18.0 :floor 3.8 :floors 2}
   :office {:w 12.8 :floor 3.6 :floors 4}
   :shops  {:w 18.0 :floor 3.6 :floors 2}
   :house  {:w 6.8 :floor 3.0 :floors 2}
   :plant  {:w 8.0 :floor 3.0 :floors 1}})

(def ^:private px-per-meter 14)

(defn- bake
  "A canvas of `painter` over `[0, w] x [t0, t1]`, s along and t up."
  [painter w t0 t1]
  (let [c (js/document.createElement "canvas")
        pw (max 2 (math/round (* w px-per-meter)))
        ph (max 2 (math/round (* (- t1 t0) px-per-meter)))
        _ (set! (.-width c) pw)
        _ (set! (.-height c) ph)
        ctx (.getContext c "2d")
        img (.createImageData ctx pw ph)
        data (.-data img)]
    (dotimes [py ph]
      (dotimes [px pw]
        (let [s (* w (/ (+ px 0.5) pw))
              t (- t1 (* (- t1 t0) (/ (+ py 0.5) ph)))
              [r g b] (painter s t)
              o (* 4 (+ px (* pw py)))]
          (aset data o (* 255 r)) (aset data (+ o 1) (* 255 g))
          (aset data (+ o 2) (* 255 b)) (aset data (+ o 3) 255))))
    (.putImageData ctx img 0 0)
    c))

(defn- texture [canvas]
  (doto (THREE/CanvasTexture. canvas)
    (set! -wrapS THREE/RepeatWrapping)
    (set! -wrapT THREE/RepeatWrapping)
    (set! -anisotropy 4)))

(defn- facade-materials
  "Two materials a style: `[style :ground]` and `[style :upper]`."
  [seed]
  (into {}
        (mapcat (fn [[style {:keys [w floor floors]}]]
                  ;; A tall wall, so the upper layers are all there; the
                  ;; ground floor and the floors above it are cut from it.
                  (let [p (facade/painter style w 200.0 floor seed)]
                    [[[style :ground] (THREE/MeshLambertMaterial.
                                       #js {:map (texture (bake p w 0.0 floor)) :side THREE/DoubleSide})]
                     [[style :upper] (THREE/MeshLambertMaterial.
                                      #js {:map (texture (bake p w floor (* floor (inc floors)))) :side THREE/DoubleSide})]]))
                tiles)))

;; ---------------------------------------------------------------------------
;; Geometry

(defn- buffers [] {:pos #js [] :nrm #js [] :uv #js []})

(defn- push-tri!
  "One triangle of world points `[x y z]` (z up), with a normal and uvs."
  [{:keys [^js pos ^js nrm ^js uv]} pts n uvs]
  (let [[nx ny nz] n]
    (doseq [[[x y z] [u v]] (map vector pts uvs)]
      ;; Ours is x east, y north, z up; three's is y up, z south.
      (.push pos x z (- y))
      (.push nrm nx nz (- ny))
      (.push uv u v))))

(defn- quad! [b [p0 p1 p2 p3] n [u0 v0 u1 v1]]
  (push-tri! b [p0 p1 p2] n [[u0 v0] [u1 v0] [u1 v1]])
  (push-tri! b [p0 p2 p3] n [[u0 v0] [u1 v1] [u0 v1]]))

(defn- geometry [{:keys [^js pos ^js nrm ^js uv]}]
  (doto (THREE/BufferGeometry.)
    (.setAttribute "position" (THREE/Float32BufferAttribute. pos 3))
    (.setAttribute "normal" (THREE/Float32BufferAttribute. nrm 3))
    (.setAttribute "uv" (THREE/Float32BufferAttribute. uv 2))))

(defn- walls!
  "The walls of a prism standing on `ground`, split into the ground floor
  and the floors above it when it stands on the street, each side's uvs
  in tiles of its style."
  [bufs style {:keys [footprint z h]} ground floor-tile]
  (let [{:keys [w floor]} (tiles style (tiles :office))
        upper-h (* floor (:floors (tiles style (tiles :office))))
        z0 (+ ground z) z1 (+ ground z h)
        split (if (zero? z) (min z1 (+ z0 floor)) z0)]
    (doseq [[[ax ay] [bx by]] (poly/edges footprint)]
      (let [len (poly/distance [ax ay] [bx by])
            n (poly/unit [(- by ay) (- ax bx)])
            n [(first n) (second n) 0.0]
            u1 (/ len w)]
        (when (> split z0)
          (quad! (bufs [style :ground]) [[ax ay z0] [bx by z0] [bx by split] [ax ay split]] n
                 [0.0 0.0 u1 (/ (- split z0) floor)]))
        (when (> z1 split)
          (quad! (bufs (if floor-tile [style :upper] :plain))
                 [[ax ay split] [bx by split] [bx by z1] [ax ay z1]] n
                 ;; Upper floors counted from the first floor up.
                 [0.0 (/ (- split z0 (if (zero? z) floor 0.0)) upper-h) u1 (/ (- z1 z0 (if (zero? z) floor 0.0)) upper-h)]))))))

(defn- cap!
  "A flat convex lid at height `z`."
  [b pts z]
  (let [c (first pts)]
    (doseq [[p q] (partition 2 1 (rest pts))]
      (push-tri! b [(conj c z) (conj p z) (conj q z)] [0.0 0.0 1.0] [[0 0] [0 0] [0 0]]))))

(defn- roof!
  "A gable or a hipped roof over the rectangle `[a b c d]`, `a`-`b` the
  long side."
  [b {:keys [kind rect z h]} ground]
  (let [[p0 p1 p2 p3] rect
        z (+ ground z) top (+ z h)
        m0 (poly/lerp p0 p3 0.5) m1 (poly/lerp p1 p2 0.5)
        short (poly/distance p0 p3)
        ;; A hip pulls the ridge in from each end by half the width.
        [r0 r1] (if (= kind :hip)
                  (let [u (poly/unit (poly/sub m1 m0)) k (min (* 0.5 short) (* 0.45 (poly/distance m0 m1)))]
                    [(poly/add m0 (poly/scale u k)) (poly/sub m1 (poly/scale u k))])
                  [m0 m1])
        at (fn [p zz] (conj (vec p) zz))
        face! (fn [pts]
                (let [[a bb c] pts
                      e1 (mapv - bb a) e2 (mapv - c a)
                      [x y zz] [(- (* (e1 1) (e2 2)) (* (e1 2) (e2 1)))
                                (- (* (e1 2) (e2 0)) (* (e1 0) (e2 2)))
                                (- (* (e1 0) (e2 1)) (* (e1 1) (e2 0)))]
                      l (max 1e-9 (math/sqrt (+ (* x x) (* y y) (* zz zz))))
                      n (if (neg? zz) [(- (/ x l)) (- (/ y l)) (- (/ zz l))] [(/ x l) (/ y l) (/ zz l)])]
                  (push-tri! b pts n [[0 0] [0 0] [0 0]])))]
    (face! [(at p0 z) (at p1 z) (at r1 top)]) (face! [(at p0 z) (at r1 top) (at r0 top)])
    (face! [(at p2 z) (at p3 z) (at r0 top)]) (face! [(at p2 z) (at r0 top) (at r1 top)])
    (face! [(at p1 z) (at p2 z) (at r1 top)]) (face! [(at p3 z) (at p0 z) (at r0 top)])))

(defn- antenna! [b {:keys [at z h]} ground]
  (let [[x y] at r 0.5
        sq [[(- x r) (- y r)] [(+ x r) (- y r)] [(+ x r) (+ y r)] [(- x r) (+ y r)]]]
    (doseq [[[ax ay] [bx by]] (poly/edges sq)]
      (let [n (poly/unit [(- by ay) (- ax bx)])]
        (quad! b [[ax ay (+ ground z)] [bx by (+ ground z)] [bx by (+ ground z h)] [ax ay (+ ground z h)]]
               [(first n) (second n) 0.0] [0 0 0 0])))))

(defn- build-city
  "Everything that does not depend on `detail`: the land, the roads, the
  lots and each lot's stages of growth."
  [seed pattern]
  (let [l (land/land {:seed seed :pattern pattern})
        r (roads/grow l {:seed seed})
        blocks (lots/blocks (roads/graph r))
        ls (lots/lots blocks (:zone l) {:seed seed})]
    {:land l :roads r :blocks blocks
     :lots (mapv (fn [lot]
                   {:lot lot
                    :ground (reduce min (map (fn [[x y]] ((:elevation l) x y)) (:points lot)))
                    :stages (buildings/rounds lot (:max-height l) {:seed seed})})
                 ls)}))

;; ---------------------------------------------------------------------------
;; The scene

(defn- terrain-mesh [l]
  (let [n 160 size (:size l)
        geo (THREE/PlaneGeometry. (* 2 size) (* 2 size) n n)
        _ (.rotateX geo (- (/ math/PI 2)))
        pos (.. geo -attributes -position)
        colors (js/Float32Array. (* 3 (.-count pos)))]
    (dotimes [i (.-count pos)]
      (let [x (.getX pos i) y (- (.getZ pos i))
            water? ((:water? l) x y)
            e ((:elevation l) x y)
            t (min 1.0 (/ e 70.0))]
        (.setY pos i (if water? -4.0 e))
        (aset colors (* 3 i) (+ 0.33 (* 0.2 t)))
        (aset colors (+ 1 (* 3 i)) (+ 0.42 (* 0.15 t)))
        (aset colors (+ 2 (* 3 i)) (+ 0.25 (* 0.12 t)))))
    (.setAttribute geo "color" (THREE/BufferAttribute. colors 3))
    (.computeVertexNormals geo)
    (THREE/Mesh. geo (THREE/MeshLambertMaterial. #js {:vertexColors true}))))

(defn- road-mesh [l r]
  (let [b (buffers)
        elev (fn [[x y]] (max 0.0 ((:elevation l) x y)))]
    (doseq [{:keys [a b-pt kind bridge?] :as s} (map #(assoc % :b-pt (:b %)) (roads/segments r))]
      (let [half (if (= kind :highway) 8.0 4.5)
            d (poly/unit (poly/sub b-pt a))
            side (poly/scale [(- (second d)) (first d)] half)
            len (poly/distance a b-pt)
            pieces (max 1 (math/ceil (/ len 25.0)))
            ;; A bridge runs straight between its ends, above the water.
            za (+ 0.6 (elev a)) zb (+ 0.6 (elev b-pt))]
        (dotimes [k pieces]
          (let [p (poly/lerp a b-pt (/ k pieces)) q (poly/lerp a b-pt (/ (inc k) pieces))
                zp (if bridge? (+ za (* (- zb za) (/ k pieces)) 3.0) (+ 0.6 (elev p)))
                zq (if bridge? (+ za (* (- zb za) (/ (inc k) pieces)) 3.0) (+ 0.6 (elev q)))
                [p1 p2] [(poly/sub p side) (poly/add p side)]
                [q1 q2] [(poly/sub q side) (poly/add q side)]]
            (quad! b [(conj p1 zp) (conj q1 zq) (conj q2 zq) (conj p2 zp)] [0.0 0.0 1.0]
                   [0 0 0 0])))
        s))
    (THREE/Mesh. (geometry b) (THREE/MeshLambertMaterial. #js {:color 0x4a4d52 :side THREE/DoubleSide}))))

(defn- building-meshes
  "The buildings of every lot after `detail` rounds, one mesh per
  material."
  [{:keys [lots]} detail facades? materials]
  (let [bufs (into {:plain (buffers) :massing (buffers) :roof (buffers) :tile-roof (buffers)}
                   (map (fn [k] [k (buffers)]) (keys materials)))]
    (doseq [{:keys [stages ground]} lots
            m (buildings/at-detail stages detail)]
      (cond
        (:sym m) (do (walls! (fn [_] (:massing bufs)) :office m ground false)
                     (cap! (:massing bufs) (:footprint m) (+ ground (:z m) (:h m))))
        (= :prism (:op m)) (let [style (:style m)]
                             (if facades?
                               (walls! bufs style m ground true)
                               (walls! (fn [_] (:plain bufs)) style m ground false))
                             (cap! (:roof bufs) (:footprint m) (+ ground (:z m) (:h m))))
        (= :roof (:op m)) (roof! (:tile-roof bufs) m ground)
        (= :antenna (:op m)) (antenna! (:roof bufs) m ground)))
    (into []
          (keep (fn [[k b]]
                  (when (pos? (.-length ^js (:pos b)))
                    (THREE/Mesh. (geometry b)
                                 (or (get materials k)
                                     (THREE/MeshLambertMaterial.
                                      #js {:color (case k :massing 0xb9bcc2 :roof 0x6f7176 :tile-roof 0x8e4a36 0xc8c3b8)
                                           :side THREE/DoubleSide}))))))
          bufs)))

(defn init! [^js container]
  (let [scene (THREE/Scene.)
        camera (THREE/PerspectiveCamera. 45 1 2 12000)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        caption (js/document.createElement "div")
        running? (atom false)
        dirty? (atom true)
        tick-fps! (fps/meter! container)
        city (atom nil)
        materials (atom nil)
        static (atom [])
        dynamic (atom [])
        sky 0xaec6dc]
    (set! (.-background scene) (THREE/Color. sky))
    (set! (.-fog scene) (THREE/Fog. sky 2500 7000))
    (.setPixelRatio renderer (or js/window.devicePixelRatio 1))
    (.appendChild container (.-domElement renderer))
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;padding:3px 8px;border-radius:4px;background:rgba(10,14,22,0.7);"
               "font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;pointer-events:none"))
    (.appendChild container caption)
    (.add scene (THREE/HemisphereLight. 0xe8f0ff 0x4a4030 0.7))
    (let [sun (THREE/DirectionalLight. 0xfff1d8 0.75)]
      (.set (.-position sun) -600 900 400)
      (.add scene sun))
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (set! (.-maxPolarAngle orbit) 1.5)
      (letfn [(clear! [meshes]
                (doseq [^js m @meshes]
                  (.remove scene m)
                  (.dispose (.-geometry m)))
                (reset! meshes []))
              (show-buildings! []
                (clear! dynamic)
                (let [ms (building-meshes @city (.-detail controls) (.-facades controls) @materials)]
                  (doseq [m ms] (.add scene m))
                  (reset! dynamic ms))
                (reset! dirty? true))
              (regenerate! []
                (set! (.-textContent caption) "growing the city...")
                (js/setTimeout
                 (fn []
                   (let [t0 (js/performance.now)
                         c (build-city (.-seed controls) (keyword (.-pattern controls)))
                         l (:land c)
                         t1 (js/performance.now)
                         water (doto (THREE/Mesh. (THREE/PlaneGeometry. (* 2 (:size l)) (* 2 (:size l)))
                                                  (THREE/MeshLambertMaterial. #js {:color 0x3d6b8f}))
                                 (.rotateX (- (/ math/PI 2))))]
                     (.set (.-position water) 0 -0.5 0)
                     (reset! city c)
                     (clear! static)
                     (let [ms [(terrain-mesh l) water (road-mesh l (:roads c))]]
                       (doseq [m ms] (.add scene m))
                       (reset! static ms))
                     (when-not @materials (reset! materials (facade-materials 7)))
                     (show-buildings!)
                     (let [[cx cy] (:center l)]
                       (.set (.-target orbit) cx 20 (- cy))
                       (.set (.-position camera) (+ cx 900) 650 (- 700 cy))
                       (.update orbit))
                     (set! (.-textContent caption)
                           (str (count (roads/segments (:roads c))) " road segments, " (count (:blocks c)) " blocks, "
                                (count (:lots c)) " lots and buildings in " (.toFixed (- t1 t0) 0) " ms"))))
                 30))
              (on-resize []
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h)
                    (reset! dirty? true))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [moved? (.update orbit)]
                    (when (or moved? @dirty?)
                      (reset! dirty? false)
                      (let [t0 (js/performance.now)]
                        (.render renderer scene camera)
                        (tick-fps! (- (js/performance.now) t0)))))))]
        (.addEventListener orbit "change" #(reset! dirty? true))
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (let [gui (GUI. #js {:container container})]
          (-> (.add gui controls "seed" 1 99 1) (.onFinishChange regenerate!))
          (-> (.add gui controls "pattern" #js ["mixed" "grid" "radial" "organic" "terrain"]) (.onChange regenerate!))
          (-> (.add gui controls "detail" 0 8 1) (.onChange show-buildings!))
          (-> (.add gui controls "facades") (.onChange show-buildings!)))
        (on-resize)
        (regenerate!)
        {:start (fn [] (when-not @running? (reset! running? true) (on-resize) (animate)))
         :stop  (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "city")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
