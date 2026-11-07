(ns allgo.demo.vision
  "Canvas demo for `allgo.vision`: the computer vision of Russell and
  Norvig, *Artificial Intelligence: A Modern Approach*, 4th edition,
  chapter 25, each method a step at a time on images drawn so that the
  right answer is known.

  Edges: Canny's stages -- smoothing, gradient, thinning, and the edges
  growing out from the strong ones. Segmentation: the regions merging an
  edge at a time, cheapest first; or k-means on the colors. Optical flow:
  a textured disk moving over a drifting backdrop, and the flow Lucas and
  Kanade find, as arrows. Stereo: the left and right views of boxes at
  different depths and the disparity found a row at a time, brighter
  nearer. Shading: an ellipsoid's normals found from its brightness and
  outline, an iteration a step, shown relit from another side beside the
  true ones. Texture: the texels of a tiled ground, each ellipse and the
  normal it implies. Contour: overlapping shapes, their T-junctions, and
  the depth order those give. Structure from motion: points on a house
  tracked through views of it turning, then factored into its shape,
  shown turning against the truth. Detection: a convolutional network
  learning disks from squares, a few examples a frame; at the end of each
  epoch it is slid over a new scene, and after a few, the windows it
  wrongly fires on are mined and added to its negatives. Selective
  search: regions merging a step, every one a proposal. Anchors and
  regression: objectness learned an epoch a frame, then the box
  regression an example at a time, and the proposals moved onto their
  objects. Boundaries: Canny's edges on regions told apart by texture,
  then the learned Pb's. Normalized cuts: superpixels, then the image
  cut where the second eigenvector says, a cut a step. Time to contact:
  flying at a textured wall, the flow, its focus of expansion and the
  time to contact it gives. Horizon: pedestrian detections placing the
  horizon and the horizon ruling out the impossible ones. Light and
  color: a shiny sphere and colored patches under orange light, the
  light's color estimated and divided out beside the scene under white
  light."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.graph :as graph]
            [allgo.learning.neural :as neural]
            [allgo.learning.unsupervised :as unsupervised]
            [allgo.random :as random]
            [allgo.vision.boundaries :as bd]
            [allgo.vision.camera :as camera]
            [allgo.vision.detection :as det]
            [allgo.vision.edges :as edges]
            [allgo.vision.flow :as flow]
            [allgo.vision.image :as img]
            [allgo.vision.proposals :as proposals]
            [allgo.vision.segmentation :as seg]
            [allgo.vision.shape :as shape]
            [allgo.vision.stereo :as stereo]
            [allgo.vision.structure :as structure]
            ["lil-gui" :default GUI]
            [clojure.math :as math]))

(def ^:private ^js controls #js {:method "edges" :seed 1 :speed 3})

;; ---------------------------------------------------------------------------
;; Scenes

(defn- scene
  "A color scene of a few shapes over a faint pattern, with a little noise."
  [rng]
  (let [palette [[0.85 0.35 0.25] [0.3 0.75 0.4] [0.95 0.8 0.3] [0.35 0.5 0.9] [0.8 0.4 0.8]]
        back (img/pattern rng 4 20.0)
        shapes (for [i (range 5)]
                 (let [c (palette i) x (+ 8 (* 80 (rng))) y (+ 8 (* 48 (rng)))]
                   (if (< (rng) 0.5)
                     {:kind :disk :center [x y] :radius (+ 6 (* 8 (rng))) :value c}
                     (let [a (+ 6 (* 10 (rng))) b (+ 6 (* 8 (rng)))]
                       {:kind :rect :box [(- x a) (- y b) (+ x a) (+ y b)] :value c}))))]
    (img/add-noise (img/render 96 64 (fn [x y] (let [v (+ 0.1 (* 0.12 (back x y)))] [v v (+ v 0.05)])) shapes) 0.02 rng)))

(def ^:private det-size 12)

(defn- det-scene [rng]
  (let [slots (vec (random/shuffle rng (for [x [8 20 32 44] y [8 22]] [x y])))
        disks (mapv (fn [[x y]] [(+ x (- (* 2 (rng)) 1.0)) (+ y (- (* 2 (rng)) 1.0))]) (take 3 slots))
        squares (take 3 (drop 3 slots))
        shapes (concat (map (fn [c] {:kind :disk :center c :radius (+ 3.5 (rng)) :value 0.9}) disks)
                       (map (fn [[x y]] (let [r (+ 3 (rng))] {:kind :rect :box [(- x r) (- y r) (+ x r) (+ y r)] :value 0.9})) squares))]
    {:image (img/add-noise (img/render 56 32 0.1 shapes) 0.05 rng)
     :boxes (mapv (fn [[x y]] (let [x0 (math/round (- x 6)) y0 (math/round (- y 6))] [x0 y0 (+ x0 det-size) (+ y0 det-size)])) disks)}))

;; ---------------------------------------------------------------------------
;; Scenes for shape, structure from motion and proposals

(def ^:private light (v/normalize [0.3 -0.2 1.0]))
(def ^:private relight (v/normalize [-0.8 -0.3 0.5]))

(defn- ellipsoid
  "A matte ellipsoid of semi-axes `a` and `b` pixels across and `c` deep,
  albedo 1, lit from `light`: its image, mask and true normals."
  [w h a b c]
  (let [cx (/ w 2.0) cy (/ h 2.0)
        normal (fn [x y] (let [dx (/ (- (+ x 0.5) cx) a) dy (/ (- (+ y 0.5) cy) b) q (+ (* dx dx) (* dy dy))]
                           (when (< q 1.0) (v/normalize [(/ dx a) (/ dy b) (/ (math/sqrt (- 1.0 q)) c)]))))]
    {:image (img/image w h (fn [x y] (if-let [n (normal x y)] (max 0.0 (v/dot n light)) 0.0)))
     :mask (img/image w h (fn [x y] (some? (normal x y))))
     :truth (img/image w h normal)}))

(def ^:private texture-camera
  {:f 120.0 :center [80.0 60.0] :camera (camera/look-at [0.0 3.0 0.0] [0.0 0.0 6.0] 120.0 [80.0 60.0])})

(defn- ground
  "The camera's view of the ground plane, round texels on a lattice."
  [cam w h spacing radius]
  (img/image w h (fn [x y]
                   (let [[o d] (camera/ray cam [(+ x 0.5) (+ y 0.5)])]
                     (if (>= (d 1) -1e-9)
                       0.15
                       (let [[px _ pz] (v/add-scaled o d (/ (- (o 1)) (d 1)))
                             gx (* spacing (math/round (/ px spacing))) gz (* spacing (math/round (/ pz spacing)))]
                         (if (<= (math/hypot (- px gx) (- pz gz)) radius) 0.9 0.25)))))))

(defn- overlapping
  "Five shapes drawn back to front, each partly over the one before: a
  label image, the labels 1 to 5 in drawing order, 0 the background."
  [rng]
  (let [shapes (loop [out [] x 14.0]
                 (if (= 5 (count out))
                   out
                   (let [y (+ 22 (* 22 (rng))) s (+ 9 (* 6 (rng)))
                         shape (if (even? (count out))
                                 {:kind :disk :center [x y] :radius s}
                                 {:kind :rect :box [(- x s) (- y (* 0.8 s)) (+ x s) (+ y (* 0.8 s))]})]
                     (recur (conj out (assoc shape :value (inc (count out)))) (+ x (* 1.3 s))))))]
    (img/render 96 64 0 shapes)))

(def ^:private house
  "A house of points: a box's corners, the ridge of its roof, a door."
  (vec (concat (for [x [-1.0 1.0] y [-0.6 0.6] z [-0.8 0.8]] [x y z])
               [[-1.0 -1.3 0.0] [1.0 -1.3 0.0]]
               [[-0.25 0.6 0.8] [0.25 0.6 0.8] [0.25 0.05 0.8] [-0.25 0.05 0.8]])))

(def ^:private house-edges
  (concat (for [i (range 8) j (range i 8)
                :when (= 1 (count (filter false? (map = (house i) (house j)))))]
            [i j])
          [[8 9] [8 0] [8 1] [9 4] [9 5] [10 13] [13 12] [12 11]]))

(defn- house-views
  "Orthographic views of the house turning about a near-vertical axis,
  a little noise in the images: `[{:image :axes}]`."
  [rng]
  (let [axis (v/normalize [0.15 1.0 0.25]) tilt (q/from-axis-angle [1.0 0.0 0.0] 0.35)]
    (vec (for [f (range 30)
               :let [r (q/mul (q/from-axis-angle axis (* f 0.07)) tilt)
                     i (q/rotate r [1.0 0.0 0.0]) j (q/rotate r [0.0 1.0 0.0])]]
           {:image (mapv (fn [p] [(+ (* 14.0 (v/dot i p)) 40.0 (* 0.2 f) (random/gaussian rng 0.0 0.15))
                                  (+ (* 14.0 (v/dot j p)) 30.0 (random/gaussian rng 0.0 0.15))])
                         house)
            :axes [i j]}))))

(defn- two-toned
  "Three objects apart from one another, each lighter above and darker
  below, on a shaded background with noise: `{:image :truth}`."
  [rng]
  (let [w 80 h 60
        candidate (fn [] (let [s (+ 10 (* 14 (rng))) s2 (* s (+ 0.6 (* 0.8 (rng))))
                               x (+ 2 (* (- w s 4) (rng))) y (+ 2 (* (- h s2 4) (rng)))]
                           {:box [x y (+ x s) (+ y s2)] :disk? (< (rng) 0.5)}))
        grown (fn [[x0 y0 x1 y1]] [(- x0 3) (- y0 3) (+ x1 3) (+ y1 3)])
        objects (loop [out []]
                  (if (= 3 (count out))
                    out
                    (let [c (candidate)]
                      (recur (if (every? #(zero? (det/iou (grown (:box %)) (:box c))) out) (conj out c) out)))))
        truth (mapv (fn [{[x0 y0 x1 y1] :box disk? :disk?}]
                      (if disk?
                        (let [r (/ (min (- x1 x0) (- y1 y0)) 2) cx (/ (+ x0 x1) 2) cy (/ (+ y0 y1) 2)]
                          [(- cx r) (- cy r) (+ cx r) (+ cy r)])
                        [x0 y0 x1 y1]))
                    objects)
        shapes (mapcat (fn [{disk? :disk?} [x0 y0 x1 y1]]
                         (let [v (+ 0.6 (* 0.35 (rng))) v2 (- v 0.25) ym (+ y0 (* (- y1 y0) (+ 0.35 (* 0.3 (rng)))))]
                           (if disk?
                             [{:kind :disk :center [(/ (+ x0 x1) 2) (/ (+ y0 y1) 2)] :radius (/ (- x1 x0) 2)
                               :texture (fn [_ y] (if (< y (/ (+ y0 y1) 2)) v v2))}]
                             [{:kind :rect :box [x0 y0 x1 ym] :value v} {:kind :rect :box [x0 ym x1 y1] :value v2}])))
                       objects truth)]
    {:image (img/add-noise (img/render w h (fn [x _] (+ 0.1 (* 0.15 (/ x w)))) shapes) 0.03 rng)
     :truth truth}))

(defn- found [truth boxes t] (count (filter (fn [g] (some #(>= (det/iou g %) t) boxes)) truth)))

(defn- textured
  "Horizontal stripes behind a square of vertical ones, a square of dots
  and a plain disk: `[image labels]`."
  [rng]
  (let [stripes (fn [vertical?] (fn [x y] (if (< (mod (if vertical? x y) 4) 2) 0.25 0.75)))
        dots (fn [x y] (if (and (< (mod x 5) 2) (< (mod y 5) 2)) 0.9 0.4))
        shapes (vec (for [t [(stripes true) dots 0.8]]
                      (let [x (+ 4 (long (* (rng) 24))) y (+ 4 (long (* (rng) 24))) s (+ 12 (long (* (rng) 8)))]
                        (if (fn? t)
                          {:kind :rect :box [x y (+ x s) (+ y s)] :texture t}
                          {:kind :disk :center [(+ x 6) (+ y 6)] :radius 7 :value t}))))]
    [(img/add-noise (img/render 48 48 (stripes false) shapes) 0.03 rng)
     (img/render 48 48 0 (map-indexed (fn [i s] (-> s (dissoc :texture) (assoc :value (inc i)))) shapes))]))

(def ^:private pb-thetas (mapv #(* % (/ math/PI 8)) (range 8)))

(defn- gray-shapes
  "Three gray shapes on a darker ground, with noise."
  [rng]
  (img/add-noise (img/render 48 36 0.15 [{:kind :rect :box [(+ 3 (random/below rng 6)) (+ 3 (random/below rng 6)) (+ 18 (random/below rng 6)) (+ 24 (random/below rng 8))] :value 0.75}
                                         {:kind :disk :center [(+ 32 (random/below rng 6)) (+ 14 (random/below rng 8))] :radius (+ 6 (random/below rng 3)) :value 0.45}
                                         {:kind :rect :box [(+ 24 (random/below rng 4)) (+ 26 (random/below rng 3)) (+ 44 (random/below rng 3)) (+ 33 (random/below rng 3))] :value 0.6}])
                 0.05 rng))

(defn- wall
  "A textured wall seen from `z`, `z0` the distance its pattern is drawn
  at: nearer, it fills more of the view."
  [pat z0 z]
  (let [s (/ z0 z)]
    (img/image 64 64 (fn [x y] (pat (+ 32 (/ (- (+ x 0.5) 32) s)) (+ 32 (/ (- (+ y 0.5) 32) s)))))))

(defn- street
  "Detections of pedestrians 1.75 m tall from a camera 1.6 m up with focal
  length 300 and its horizon at row `v0`, and two no pedestrian could
  make."
  [rng v0]
  (concat (for [z (repeatedly 10 #(+ 4.0 (* 26 (rng))))]
            {:foot (random/gaussian rng (+ v0 (/ (* 300 1.6) z)) 0.5) :height (random/gaussian rng (/ (* 300 1.75) z) 0.5) :x (* 150 (rng))})
          [{:foot (+ v0 12) :height 70.0 :x (* 150 (rng))} {:foot (+ v0 100) :height 15.0 :x (* 150 (rng))}]))

(defn- lit-scene
  "Colored patches facing the viewer and a shiny sphere before them, lit
  by `color` from `dir`: patches of the given albedos in a grid."
  [albedos dir color]
  (img/image 72 48 (fn [x y]
                     (let [dx (/ (- x 36) 16.0) dy (/ (- y 24) 16.0) r2 (+ (* dx dx) (* dy dy))]
                       (if (< r2 1.0)
                         (camera/shade {:albedo [0.85 0.85 0.85] :specular 0.4 :shininess 25} [dx dy (math/sqrt (- 1.0 r2))] dir [0.0 0.0 1.0] color 0.1)
                         (camera/shade {:albedo (albedos (+ (* 6 (quot y 8)) (quot x 12)))} [0.0 0.0 1.0] dir [0.0 0.0 1.0] color 0.1))))))

(defn- jitter
  "A box near `g`, overlapping it by at least a half."
  [rng g]
  (let [wobble #(- (* 2 (rng)) 1)]
    (first (filter #(>= (det/iou % g) 0.5)
                   (repeatedly (fn [] (let [[x0 y0 x1 y1] g w (- x1 x0) h (- y1 y0)
                                            cx (+ (/ (+ x0 x1) 2) (* w 0.25 (wobble))) cy (+ (/ (+ y0 y1) 2) (* h 0.25 (wobble)))
                                            sw (* w (math/exp (* 0.3 (wobble)))) sh (* h (math/exp (* 0.3 (wobble))))]
                                        [(- cx (/ sw 2)) (- cy (/ sh 2)) (+ cx (/ sw 2)) (+ cy (/ sh 2))])))))))

;; ---------------------------------------------------------------------------
;; Each method as a seq of frames to show

(defn- region-colors
  "A color for each region root, by hashing it."
  [root]
  (let [h (hash root)]
    [(+ 0.25 (* 0.7 (/ (mod h 97) 97.0))) (+ 0.25 (* 0.7 (/ (mod (quot h 97) 89) 89.0))) (+ 0.25 (* 0.7 (/ (mod (quot h 8633) 83) 83.0)))]))

(defn- frames [method rng]
  (case method
    "edges"
    (let [color (scene rng)
          gray (img/gray color)]
      (map (fn [st] {:left color :stage st}) (edges/canny-seq gray :sigma 1.2 :low 0.03 :high 0.08)))

    "segmentation"
    (let [color (scene rng)
          smoothed (img/smooth color 0.8)
          states (vec (seg/graph-based-seq smoothed 1.5))
          ;; every 200th edge, then the specks tidied into their neighbors
          shown (concat (take-nth 200 states) [(peek states) (seg/merge-small (peek states) smoothed 20)])]
      (map (fn [st] {:left color :forest st :smoothed smoothed}) shown))

    "k-means"
    (let [color (scene rng)
          points (vec (:px color))]
      (map (fn [st] {:left color :kmeans st}) (take 15 (rest (unsupervised/k-means-seq points 6 rng)))))

    "optical flow"
    (let [back (img/pattern rng 6 10.0) disk (img/pattern rng 6 5.0)
          frame (fn [t] (img/image 72 48 (fn [x y]
                                           (let [dx (- x (+ 20 (* 1.2 t))) dy (- y (+ 18 (* 0.6 t)))]
                                             (if (< (+ (* dx dx) (* dy dy)) 100)
                                               (disk (- x (* 1.2 t)) (- y (* 0.6 t)))
                                               (* 0.6 (back (+ x (* 0.5 t)) y)))))))]
      (map (fn [t] (let [a (frame t) b (frame (inc t))]
                     {:left a :flow (flow/flow-field a b 6 6 4 3 :lucas-kanade)}))
           (range 0 28)))

    "stereo"
    (let [layers [{:texture (img/pattern rng 8 6.0) :disparity 1}
                  {:texture (img/pattern rng 8 4.0) :disparity 5 :box [8 8 30 40]}
                  {:texture (img/pattern rng 8 4.0) :disparity 9 :box [40 14 62 34]}]
          [l r] (stereo/stereo-pair layers 72 48)]
      (map (fn [d] {:left l :right r :disparity d}) (stereo/disparity-seq l r 3 12)))

    "shading"
    (let [{:keys [image mask truth]} (ellipsoid 72 56 30.0 20.0 24.0)]
      (map (fn [st] {:left image :mask mask :truth truth :shading st})
           (take-nth 8 (take 801 (shape/shading-seq image mask 1.0 light)))))

    "texture"
    (let [{:keys [f center camera]} texture-camera
          g (ground camera 160 120 1.0 0.3)
          ts (vec (sort-by (comp second :pixel) (shape/texels g f center 0.6 40)))
          up (camera/to-camera (assoc camera :position [0.0 0.0 0.0]) [0.0 1.0 0.0])]
      (map (fn [n] {:left g :texels (subvec ts 0 n) :all ts :up up}) (range (inc (count ts)))))

    "contour"
    (let [labels (overlapping rng)
          votes (vec (shape/occlusions labels 4))]
      (map (fn [n] {:labels labels :votes (subvec votes 0 n) :pairs (when (= n (count votes)) (shape/in-front labels 4))})
           (range (inc (count votes)))))

    "structure from motion"
    (let [vs (house-views rng)
          frames (mapv :image vs)
          [i1 j1] (:axes (first vs))
          k1 (v/cross i1 j1)
          truth (let [ps (mapv (fn [p] [(v/dot i1 p) (v/dot j1 p) (v/dot k1 p)]) house)
                      m (v/scale (reduce v/add ps) (/ 1.0 (count ps)))]
                  (mapv #(v/scale (v/sub % m) 14.0) ps))
          states (structure/factorization-seq frames)
          {:keys [points residual]} (structure/factorize frames)
          ;; the mirror image chosen to match the truth: the views cannot tell
          mirror (if (neg? (reduce + (map #(* (%1 2) (%2 2)) points truth))) -1.0 1.0)
          recovered (mapv (fn [[x y z]] [x y (* mirror z)]) points)]
      (concat (for [f (range (count frames))] {:frames frames :f f})
              (for [st states] {:frames frames :f (dec (count frames)) :stage st})
              (for [t (range 126)] {:frames frames :f (dec (count frames)) :recovered recovered :truth truth
                                    :mirror mirror :residual residual :angle (* t 0.05)})))

    "selective search"
    (let [{:keys [image truth]} (two-toned rng)
          labels (seg/graph-based image 3 8)
          states (proposals/selective-search-seq image labels)
          ;; which region each starting one is now part of
          parents (rest (reductions (fn [m st] (if-let [[a b n] (:merged st)] (assoc m a n b n) m)) {} states))]
      (map (fn [st parent] {:left image :truth truth :labels labels :ss st :parent parent}) states parents))

    "boundaries"
    (let [w (nth (bd/train-seq (bd/examples (vec (repeatedly 3 #(textured rng))) 50 3 5 pb-thetas 8 1.0 rng) 1.0) 300)
          [image labels] (textured rng)
          truth (bd/true-boundary labels)]
      ;; Pb only once the edges have been shown: it takes a moment
      (cons {:left image :truth truth :edges (edges/canny image)}
            (lazy-seq [{:left image :truth truth :edges (edges/canny image) :pb (bd/thin (bd/pb w image 5 pb-thetas 8 1.0))}])))

    "normalized cuts"
    (let [image (gray-shapes rng)
          graph (seg/superpixel-graph image 0.05 8 0.1)]
      (cons {:left image :superpixels graph}
            (map (fn [st] {:left image :ncut st}) (seg/normalized-cuts-seq image {}))))

    "time to contact"
    (let [pat (img/pattern rng 12 9.0)]
      (map (fn [t] (let [z (- 20.0 (* 0.5 t)) a (wall pat 20.0 z) b (wall pat 20.0 (- z 0.5))
                         field (flow/flow-field a b 6 10 4 3 :lucas-kanade)
                         foe (flow/focus-of-expansion field)]
                     {:left a :flow field :foe foe :truth (/ (- z 0.5) 0.5)
                      :ttc (when foe (let [ts (for [[p v] field :when (> (math/hypot (- (p 0) (foe 0)) (- (p 1) (foe 1))) 8)]
                                                (flow/time-to-contact foe p v))]
                                       (/ (reduce + ts) (max 1 (count ts)))))}))
           (range 0 28)))

    "horizon"
    (let [v0 (+ 30 (* 30 (rng)))
          ds (street rng v0)]
      (concat [{:detections ds :v0 v0}]
              (map (fn [st] {:detections ds :v0 v0 :fit st}) (shape/horizon-seq ds 0.2))))

    "light and color"
    (let [albedos (vec (repeatedly 36 (fn [] (vec (repeatedly 3 #(+ 0.1 (* 0.8 (rng))))))))
          tint [1.0 0.75 0.45] white (/ (reduce + tint) 3)]
      (map (fn [t] (let [a (* t 0.12) dir (v/normalize [(math/cos a) (* 0.6 (math/sin a)) 1.2])
                         tinted (lit-scene albedos dir tint)
                         {:keys [illuminant image]} (camera/gray-world tinted)]
                     {:tinted tinted :corrected image :white (lit-scene albedos dir [white white white])
                      :illuminant illuminant :tint (mapv #(/ % white) tint)}))
           (range 52)))

    "anchors and regression" nil

    "detection" nil))

;; ---------------------------------------------------------------------------
;; Detection, trained a few examples a frame

(defn- detection-init [rng]
  (let [train (vec (repeatedly 6 #(det-scene rng)))
        ex (map #(det/training-windows (:image %) (:boxes %) det-size 10 0.2 rng) train)
        pos (vec (mapcat :positives ex)) neg (vec (mapcat :negatives ex))
        [X Y] (det/examples pos neg det-size)]
    {:train train :pos pos :neg neg :X X :Y Y :order (random/shuffle rng (range (count X))) :i 0
     :net (det/classifier det-size 4 8 rng) :epoch 0 :round 0 :test (det-scene rng) :dets nil :rng rng}))

(defn- detection-advance [{:keys [X Y order i net epoch round train pos neg test rng] :as st} n]
  (let [k (min n (- (count order) i))
        net (reduce (fn [net j] (neural/sgd-step net (X j) (Y j) 0.05 :cross-entropy)) net (subvec order i (+ i k)))
        i (+ i k)]
    (if (< i (count order))
      (assoc st :net net :i i)
      ;; an epoch done: detect on the test scene; every few, mine
      (let [epoch (inc epoch)
            dets (det/detect net (:image test) det-size 2 0.5 0.3)
            mine? (and (zero? (mod epoch 6)) (< round 2))]
        (if mine?
          (let [hard (vec (mapcat #(det/hard-negatives net (:image %) (:boxes %) det-size 2 0.3 0.2) train))
                neg (into neg hard)
                [X Y] (det/examples pos neg det-size)]
            (assoc st :net net :epoch epoch :dets dets :neg neg :X X :Y Y :i 0 :round (inc round) :mined (count hard)
                   :order (random/shuffle rng (range (count X)))))
          (assoc st :net net :epoch epoch :dets dets :i 0 :order (random/shuffle rng order)))))))

;; ---------------------------------------------------------------------------
;; Anchors and box regression, learned a little a frame

(def ^:private anchor-k 10)
(def ^:private anchor-context 1.5)

(defn- anchors-of [image] (proposals/anchors image 4 [12 18 26] [0.6 1.0 1.6]))

(defn- anchors-init [rng]
  {:phase :examples :train (vec (repeatedly 12 #(two-toned rng))) :test (two-toned rng) :rng rng
   :done 0 :X [] :ys [] :epoch 0})

(defn- propose [weights {:keys [image]}]
  (map :box (proposals/propose weights image (anchors-of image) anchor-k anchor-context 0.5 0.5)))

(defn- anchors-advance
  "A training scene's anchors sampled; or `n` epochs of objectness, the
  test scene's proposals redone every ten; or `4 n` examples added to
  the box regression, then the proposals refined."
  [{:keys [phase train test rng done X ys epoch] :as st} n]
  (case phase
    :examples (let [{:keys [image truth]} (train done)
                    [x y] (proposals/objectness-examples image truth (anchors-of image) anchor-k anchor-context 0.6 0.3 64 rng)
                    st (assoc st :X (into X x) :ys (into ys y) :done (inc done))]
                (if (= (inc done) (count train))
                  (assoc st :phase :objectness :weights (proposals/objectness-seq [(:X st) (:ys st)] 2.0))
                  st))
    :objectness (let [ws (nthrest (:weights st) n)
                      e (+ epoch n)
                      st (assoc st :weights ws :epoch e)
                      st (if (not= (quot epoch 10) (quot e 10)) (assoc st :proposals (propose (first ws) test)) st)]
                  (if (>= e 50)
                    (let [examples (vec (for [{:keys [image truth]} train g truth _ (range 10)] [image (jitter rng g) g]))]
                      (assoc st :phase :regression :final (first ws) :proposals (propose (first ws) test)
                             :equations (rest (proposals/box-regression-seq examples anchor-k anchor-context))
                             :total (count examples)))
                    st))
    :regression (let [s (drop (dec (* 4 n)) (:equations st))]
                  (if (next s)
                    (assoc st :equations (rest s) :eqs (first s))
                    (let [regressor (proposals/regressor (or (first s) (:eqs st)) 1.0)]
                      (assoc st :phase :done :eqs (or (first s) (:eqs st)) :equations nil
                             :refined (mapv #(proposals/refine regressor (:image test) %) (:proposals st))))))
    :done st))

;; ---------------------------------------------------------------------------
;; Drawing

(defn- to-rgb [v]
  (cond (vector? v) v
        (nil? v) [0.0 0.0 0.0]
        :else [v v v]))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        scratch (js/document.createElement "canvas")
        running? (atom false)
        model (atom nil)]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! []
              (let [rng (random/rng (.-seed controls)) method (.-method controls)]
                (reset! model (case method
                                "detection" {:det (detection-init rng) :clock 0.0}
                                "anchors and regression" {:anc (anchors-init rng) :clock 0.0}
                                (let [fs (frames method rng)] {:rest fs :frame (first fs) :step 0 :clock 0.0})))))
            (blit! [{:keys [w h px]} f x y side-w side-h]
              ;; an image drawn to fit side-w by side-h at (x, y), pixels square
              (set! (.-width scratch) w) (set! (.-height scratch) h)
              (let [c2 (.getContext scratch "2d") data (.createImageData c2 w h) d (.-data data)]
                (dotimes [i (* w h)]
                  (let [[r g b] (to-rgb (f (px i))) o (* 4 i)]
                    (aset d o (* 255 (max 0 (min 1 r)))) (aset d (+ o 1) (* 255 (max 0 (min 1 g))))
                    (aset d (+ o 2) (* 255 (max 0 (min 1 b)))) (aset d (+ o 3) 255)))
                (.putImageData c2 data 0 0))
              (let [s (min (/ side-w w) (/ side-h h))]
                (set! (.-imageSmoothingEnabled ctx) false)
                (.drawImage ctx scratch x y (* s w) (* s h))
                s))
            (panels [n]
              ;; n side-by-side panels below the controls: [[x y w h] ...]
              (let [cw (.-width canvas) ch (.-height canvas)
                    gap 16 pw (/ (- cw (* gap (inc n))) n) ph (- ch 200)]
                (mapv (fn [i] [(+ gap (* i (+ pw gap))) 120 pw ph]) (range n))))
            (label! [text x y]
              (set! (.-fillStyle ctx) "#9fb4d8") (set! (.-font ctx) "11px ui-monospace,Menlo,monospace")
              (.fillText ctx text x (- y 8)))
            (advance! []
              (let [method (.-method controls) n (.-speed controls)]
                (case method
                  "detection" (swap! model update :det detection-advance (* 8 n))
                  "anchors and regression" (swap! model update :anc anchors-advance n)
                  (let [{:keys [clock]} @model
                        rate (case method "segmentation" 7.0 "stereo" 12.0 "optical flow" 4.0 "shading" 6.0
                                   "structure from motion" 6.0 "selective search" 0.6 "texture" 4.0 "time to contact" 3.0
                                   "light and color" 6.0 "boundaries" 0.5 "normalized cuts" 0.8 "horizon" 0.8 2.0)
                        clock (+ clock (* rate n 0.016))]
                    (if (>= clock 1.0)
                      (swap! model (fn [{:keys [rest step] :as m}]
                                     (if (next rest)
                                       (assoc m :rest (next rest) :frame (fnext rest) :step (inc step) :clock 0.0)
                                       ;; flow loops; the rest hold their last frame
                                       (if (= method "optical flow")
                                         (let [fs (frames method (random/rng (.-seed controls)))] (assoc m :rest fs :frame (first fs) :step 0 :clock 0.0))
                                         (assoc m :clock 0.0)))))
                      (swap! model assoc :clock clock))))))
            (draw! []
              (let [method (.-method controls)
                    cw (.-width canvas) ch (.-height canvas)]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (set! (.-textContent caption)
                      (str method "  |  "
                           (case method
                             "edges" (draw-edges!)
                             "segmentation" (draw-segmentation!)
                             "k-means" (draw-kmeans!)
                             "optical flow" (draw-flow!)
                             "stereo" (draw-stereo!)
                             "detection" (draw-detection!)
                             "shading" (draw-shading!)
                             "texture" (draw-texture!)
                             "contour" (draw-contour!)
                             "structure from motion" (draw-sfm!)
                             "selective search" (draw-selective!)
                             "anchors and regression" (draw-anchors!)
                             "boundaries" (draw-boundaries!)
                             "normalized cuts" (draw-ncut!)
                             "time to contact" (draw-ttc!)
                             "horizon" (draw-horizon!)
                             "light and color" (draw-light!))))))
            (draw-boundaries! []
              (let [{{:keys [left truth edges pb]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2] [x3 y3 w3 h3]] (panels 3)
                    {:keys [w h]} left
                    near? (fn [[x y]] (some #(truth [(+ x (% 0)) (+ y (% 1))]) (for [dx (range -2 3) dy (range -2 3)] [dx dy])))
                    marks (fn [on?] (img/image w h (fn [x y] (cond (not (on? x y)) 0.06 (near? [x y]) [0.5 1.0 0.63] :else [1.0 0.56 0.44]))))
                    share (fn [ps] (let [ps (seq ps)] (if ps (str (.toFixed (* 100 (/ (count (filter near? ps)) (count ps))) 0) "%") "-")))
                    pb-pixels (when pb (for [y (range h) x (range w) :when (> (img/at pb x y) 0.5)] [x y]))]
                (blit! left identity x1 y1 w1 h1) (label! "regions told apart by texture, not brightness" x1 y1)
                (blit! (marks (fn [x y] (edges [x y]))) identity x2 y2 w2 h2) (label! "Canny's edges" x2 y2)
                (if pb
                  (do (blit! (marks (fn [x y] (> (img/at pb x y) 0.5))) identity x3 y3 w3 h3) (label! "Pb, learned, thinned" x3 y3))
                  (label! "Pb: comparing the halves of a disk round every pixel ..." x3 y3))
                (str "on a boundary (green, within two pixels of one): Canny " (share edges) (when pb (str ", Pb " (share pb-pixels)))
                     "\nPb compares the brightness and texture histograms of a disk's two halves at eight orientations; logistic regression, learned from scenes whose boundaries are known, turns the differences into a probability")))
            (draw-ncut! []
              (let [{{:keys [left superpixels ncut]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2]] (panels 2)
                    {:keys [w h]} left
                    paint (fn [groups] (let [color (into {} (for [[k ps] (map-indexed vector groups) p ps] [p (region-colors [:ncut k])]))]
                                         (img/image w h (fn [x y] (get color [x y] [0.1 0.1 0.1])))))]
                (blit! left identity x1 y1 w1 h1) (label! "the image" x1 y1)
                (if superpixels
                  (do (blit! (paint (map :pixels (:nodes superpixels))) identity x2 y2 w2 h2)
                      (label! (str (count (:nodes superpixels)) " superpixels: the graph's nodes") x2 y2))
                  (do (blit! (paint (concat (:regions ncut) (map #(mapcat :pixels (:nodes %)) (:pending ncut)))) identity x2 y2 w2 h2)
                      (label! "regions so far" x2 y2)))
                (str (if superpixels
                       "an over-segmentation first, as the chapter suggests: a few dozen nodes in place of 1,728 pixels"
                       (str (count (:regions ncut)) " regions finished, " (count (:pending ncut)) " still to try"))
                     "\nShi and Malik: each part cut where the second eigenvector of the normalized graph gives the least normalized cut, while that is under 0.05")))
            (draw-ttc! []
              (let [{{:keys [left flow foe ttc truth]} :frame} @model
                    [[x1 y1 w1 h1]] (panels 1)
                    s (blit! left identity x1 y1 w1 h1)]
                (set! (.-strokeStyle ctx) "#ffe08a") (set! (.-lineWidth ctx) 2)
                (doseq [[[x y] [u v]] flow :let [px (+ x1 (* s (+ x 0.5))) py (+ y1 (* s (+ y 0.5)))]]
                  (.beginPath ctx) (.moveTo ctx px py) (.lineTo ctx (+ px (* 4 s u)) (+ py (* 4 s v))) (.stroke ctx))
                (when foe
                  (let [[fx fy] foe px (+ x1 (* s (+ fx 0.5))) py (+ y1 (* s (+ fy 0.5)))]
                    (set! (.-strokeStyle ctx) "#7dffa0") (set! (.-lineWidth ctx) 3)
                    (.beginPath ctx) (.moveTo ctx (- px 10) py) (.lineTo ctx (+ px 10) py) (.moveTo ctx px (- py 10)) (.lineTo ctx px (+ py 10)) (.stroke ctx)))
                (str "time to contact from the flow: " (if ttc (.toFixed ttc 1) "-") " frames  |  true: " (.toFixed truth 1)
                     "\nflying at a textured wall: the flow (four times its length) radiates from the focus of expansion, green; each point's distance from it over its speed is the time to contact")))
            (draw-horizon! []
              (let [{{:keys [detections v0 fit]} :frame} @model
                    [[x1 y1 w1 h1]] (panels 1)
                    s (min (/ w1 160) (/ h1 200))
                    kept (set (:kept fit))
                    at (fn [x v] [(+ x1 (* s (+ 5 x))) (+ y1 (* s v))])]
                (set! (.-fillStyle ctx) "#161b25") (.fillRect ctx x1 y1 (* s 160) (* s 200))
                (let [[ax ay] (at -5 v0) [bx _] (at 155 v0)]
                  (set! (.-strokeStyle ctx) "rgba(255,255,255,0.3)") (set! (.-lineWidth ctx) 1)
                  (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx ay) (.stroke ctx))
                (when-let [hz (:horizon fit)]
                  (let [[ax ay] (at -5 hz) [bx _] (at 155 hz)]
                    (set! (.-strokeStyle ctx) "#ffd479") (set! (.-lineWidth ctx) 2)
                    (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx ay) (.stroke ctx)))
                (doseq [{:keys [foot height x] :as d} detections
                        :let [[px py] (at x foot) bw (* s 0.3 height)]]
                  (set! (.-strokeStyle ctx) (cond (nil? fit) "#c9d2e0" (kept d) "#7dffa0" :else "#ff8f70"))
                  (set! (.-lineWidth ctx) 2)
                  (.strokeRect ctx (- px (/ bw 2)) (- py (* s height)) bw (* s height)))
                (str (if fit
                       (str "horizon at row " (.toFixed (:horizon fit) 1) " (true " (.toFixed v0 1) ")  |  " (count (:dropped fit)) " detections ruled out")
                       (str (count detections) " pedestrian detections"))
                     "\nthe same height on the same ground: the nearer, the lower the feet and the taller -- a line, zero at the horizon. Yellow: the horizon found; faint: the true one; red: detections no pedestrian could be")))
            (draw-light! []
              (let [{{:keys [tinted corrected white illuminant tint]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2] [x3 y3 w3 h3]] (panels 3)]
                (blit! tinted identity x1 y1 w1 h1) (label! "under orange light" x1 y1)
                (blit! corrected identity x2 y2 w2 h2) (label! "gray world: the light's color divided out" x2 y2)
                (blit! white identity x3 y3 w3 h3) (label! "under white light" x3 y3)
                (str "light estimated " (apply str (interpose " " (map #(.toFixed % 2) illuminant))) "  |  true " (apply str (interpose " " (map #(.toFixed % 2) tint)))
                     "\ndiffuse surfaces with a specularity, Lambert's law a channel at a time; color constancy assumes the scene averages to gray")))
            (draw-edges! []
              (let [{{:keys [left stage]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2]] (panels 2)]
                (blit! left identity x1 y1 w1 h1) (label! "the scene" x1 y1)
                (case (:stage stage)
                  :smoothed (do (blit! (:image stage) identity x2 y2 w2 h2) (label! "smoothed by a Gaussian" x2 y2))
                  :gradient (let [m (apply max 1e-9 (:px (:image stage)))]
                              (blit! (:image stage) #(/ % m) x2 y2 w2 h2) (label! "gradient magnitude" x2 y2))
                  :thinned (let [m (apply max 1e-9 (:px (:image stage)))]
                             (blit! (:image stage) #(/ % m) x2 y2 w2 h2) (label! "thinned: maxima across the edge" x2 y2))
                  :hysteresis (let [e (:edges stage) {:keys [w h]} left]
                                (blit! (img/image w h (fn [x y] (if (e [x y]) 1.0 0.08))) identity x2 y2 w2 h2)
                                (label! "edges, grown from the strong ones" x2 y2)))
                (str (name (:stage stage)) (when-let [e (:edges stage)] (str ", " (count e) " edge pixels"))
                     "\nleft: the scene  right: Canny's current stage -- of brightness only, so colors as bright as each other meet without an edge")))
            (draw-segmentation! []
              (let [{{:keys [left forest smoothed]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2]] (panels 2)
                    {:keys [w h]} smoothed
                    ;; colored by each region's root, which stays put as
                    ;; regions merge, so the colors do not flicker
                    roots (img/image w h (fn [x y] (graph/canonical (:forest forest) [x y])))
                    colors (into {} (map (fn [r] [r (region-colors r)]) (distinct (:px roots))))]
                (blit! left identity x1 y1 w1 h1) (label! "the scene" x1 y1)
                (blit! roots colors x2 y2 w2 h2) (label! "regions, merged cheapest edge first" x2 y2)
                (str (count colors) " regions"
                     "\nright: each region a color; pixels merge across the weakest boundaries first")))
            (draw-kmeans! []
              (let [{{:keys [left kmeans]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2]] (panels 2)
                    {:keys [centers assignment cost]} kmeans]
                (blit! left identity x1 y1 w1 h1) (label! "the scene" x1 y1)
                (blit! (assoc left :px assignment) centers x2 y2 w2 h2) (label! "each pixel painted its cluster's center" x2 y2)
                (str "cost " (.toFixed cost 2)
                     "\nright: 6 clusters of color, no regard for where pixels are")))
            (draw-flow! []
              (let [{{:keys [left flow]} :frame} @model
                    [[x1 y1 w1 h1]] (panels 1)
                    s (blit! left identity x1 y1 w1 h1)]
                (set! (.-strokeStyle ctx) "#ffe08a") (set! (.-fillStyle ctx) "#ffe08a") (set! (.-lineWidth ctx) 2.5)
                (doseq [[[x y] [u v]] flow
                        :let [px (+ x1 (* s (+ x 0.5))) py (+ y1 (* s (+ y 0.5)))]]
                  (.beginPath ctx) (.moveTo ctx px py) (.lineTo ctx (+ px (* 4 s u)) (+ py (* 4 s v))) (.stroke ctx)
                  (.fillRect ctx (- px 1) (- py 1) 2 2))
                (str (count flow) " windows"
                     "\narrows: the flow by Lucas-Kanade, four times its length -- the disk moving, the backdrop drifting")))
            (draw-stereo! []
              (let [{{:keys [left right disparity]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2] [x3 y3 w3 h3]] (panels 3)]
                (blit! left identity x1 y1 w1 h1) (label! "left eye" x1 y1)
                (blit! right identity x2 y2 w2 h2) (label! "right eye" x2 y2)
                (blit! disparity (fn [d] (let [t (/ d 12.0)] [(* 1.2 t) (* 0.8 t) (- 1.0 t)])) x3 y3 w3 h3)
                (label! "disparity: nearer is warmer" x3 y3)
                "a row at a time: each pixel's best match along its row in the other eye"))
            (draw-detection! []
              (let [{{:keys [test dets epoch round i X mined]} :det} @model
                    [[x1 y1 w1 h1]] (panels 1)
                    s (blit! (:image test) identity x1 y1 w1 h1)
                    rect! (fn [[a b c d] color width]
                            (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
                            (.strokeRect ctx (+ x1 (* s a)) (+ y1 (* s b)) (* s (- c a)) (* s (- d b))))]
                (doseq [b (:boxes test)] (rect! b "rgba(255,255,255,0.35)" 1))
                (doseq [{:keys [box]} dets]
                  (rect! box (if (some #(> (det/iou box %) 0.5) (:boxes test)) "#7dffa0" "#ff8f70") 2.5))
                (str "epoch " epoch ", example " i " of " (count X)
                     (when (pos? round) (str "  |  " round " round" (when (> round 1) "s") " of hard negatives mined"
                                             (when mined (str " (" mined " last time)"))))
                     "\nfaint: the true disks  green: detections that hit one  red: false alarms among the squares")))
            (draw-shading! []
              (let [{{:keys [left mask truth shading]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2] [x3 y3 w3 h3]] (panels 3)
                    ns (shape/normals shading mask)
                    relit (fn [n] (if n (max 0.0 (v/dot n relight)) 0.0))
                    errs (keep (fn [[a b]] (when (and a b) (math/to-degrees (math/acos (min 1.0 (v/dot a b)))))) (map vector (:px ns) (:px truth)))]
                (blit! left identity x1 y1 w1 h1) (label! "the image: lit from the upper right" x1 y1)
                (blit! ns relit x2 y2 w2 h2) (label! "normals found, lit from the left" x2 y2)
                (blit! truth relit x3 y3 w3 h3) (label! "true normals, lit from the left" x3 y3)
                (str "iteration " (:iteration shading) "  |  mean error " (.toFixed (/ (reduce + errs) (max 1 (count errs))) 1) " degrees"
                     "\nIkeuchi and Horn: each normal the mean of its neighbors', moved to fit the brightness; the contour's normals, square to the view, fixed")))
            (draw-texture! []
              (let [{{:keys [left texels all up]} :frame} @model
                    {:keys [f]} texture-camera
                    [[x1 y1 w1 h1]] (panels 1)
                    s (blit! left identity x1 y1 w1 h1)
                    slant-error (fn [{:keys [ray slant]}] (abs (- (math/to-degrees slant) (math/to-degrees (math/acos (abs (v/dot up ray)))))))]
                (set! (.-lineWidth ctx) 2)
                (doseq [{[u v] :pixel :keys [major minor slant tilt]} texels
                        :let [px (+ x1 (* s u)) py (+ y1 (* s v))
                              ;; the needle along the normal's image, drawn upward
                              dir (if (pos? (math/sin tilt)) (+ tilt math/PI) tilt)
                              len (* s 9 (math/sin slant))]]
                  (set! (.-strokeStyle ctx) "#ffd479")
                  (.beginPath ctx) (.ellipse ctx px py (* 0.5 s f major) (* 0.5 s f minor) (+ tilt (/ math/PI 2)) 0 (* 2 math/PI)) (.stroke ctx)
                  (set! (.-strokeStyle ctx) "#7dffa0")
                  (.beginPath ctx) (.moveTo ctx px py) (.lineTo ctx (+ px (* len (math/cos dir))) (+ py (* len (math/sin dir)))) (.stroke ctx))
                (str (count texels) " of " (count all) " texels"
                     (when (seq texels) (str "  |  slant within " (.toFixed (apply max (map slant-error texels)) 1) " degrees of the truth"))
                     (when (= (count texels) (count all))
                       (str "  |  from their sizes alone, the ground's normal within "
                            (.toFixed (math/to-degrees (math/acos (min 1.0 (abs (v/dot up (shape/plane-from-texels all)))))) 1) " degrees"))
                     "\nyellow: each texel's ellipse, foreshortened by the cosine of the slant  green: the normal it implies, longer the more slanted")))
            (draw-contour! []
              (let [{{:keys [labels votes pairs]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2]] (panels 2)
                    gray (fn [l] (if (zero? l) 0.06 (+ 0.2 (* 0.13 l))))
                    s (blit! labels gray x1 y1 w1 h1)]
                (label! "five shapes, drawn back to front" x1 y1)
                (doseq [{[jx jy] :at :keys [front behind]} votes]
                  (set! (.-fillStyle ctx) (if (every? #(> front %) behind) "#7dffa0" "#ff8f70"))
                  (.beginPath ctx) (.arc ctx (+ x1 (* s jx)) (+ y1 (* s jy)) 4 0 (* 2 math/PI)) (.fill ctx))
                (when pairs
                  ;; each region's layer: how many regions, one in front of the next, are in front of it
                  (let [layer (fn layer [r] (reduce max 0 (for [[a b] pairs :when (= b r)] (inc (long (layer a))))))
                        deepest (reduce max 1 (map layer (range 6)))]
                    (blit! labels (fn [l] (if (zero? l) [0.03 0.03 0.12]
                                              (let [t (- 1.0 (/ (layer l) deepest))] [(+ 0.25 (* 0.7 t)) (+ 0.2 (* 0.55 t)) (+ 0.35 (* 0.45 t))])))
                           x2 y2 w2 h2)
                    (label! "depth from the junctions: the nearer, the brighter" x2 y2)))
                (str (count votes) " T-junctions decide"
                     (when pairs
                       (let [shapes (remove (fn [[_ b]] (zero? b)) pairs)]
                         (str "  |  in front: " (apply str (interpose ", " (map (fn [[a b]] (str (if (zero? a) "the background" a) " of " b)) (sort shapes))))
                              "  |  " (count (filter (fn [[a b]] (> a b)) shapes)) " of " (count shapes) " as drawn")))
                     "\ndots: green where the shape whose outline runs on through the junction was drawn later, as it should be; junctions at corners say nothing")))
            (draw-sfm! []
              (let [{{:keys [frames f stage recovered truth mirror residual angle]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2]] (panels 2)
                    s (min (/ w1 96) (/ h1 64))
                    edges! (fn [pts at color width]
                             (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
                             (doseq [[a b] house-edges :let [[ax ay] (at (pts a)) [bx by] (at (pts b))]]
                               (.beginPath ctx) (.moveTo ctx ax ay) (.lineTo ctx bx by) (.stroke ctx)))]
                (set! (.-fillStyle ctx) "#161b25") (.fillRect ctx x1 y1 (* s 96) (* s 64))
                (label! (str "view " (inc f) " of " (count frames) ": the points tracked") x1 y1)
                (set! (.-fillStyle ctx) "rgba(255,212,121,0.35)")
                (doseq [g (range (inc f)) [u v] (frames g)] (.fillRect ctx (+ x1 (* s u) -1) (+ y1 (* s v) -1) 2 2))
                (edges! (frames f) (fn [[u v]] [(+ x1 (* s u)) (+ y1 (* s v))]) "#ffd479" 1.5)
                (let [cx (+ x2 (* 0.5 w2)) cy (+ y2 (* 0.45 h2)) k (* 0.9 s)
                      turn (fn [[x y z]] (let [c (math/cos angle) sn (math/sin angle)] [(+ cx (* k (+ (* c x) (* sn z)))) (+ cy (* k y))]))]
                  (cond recovered (do (edges! truth turn "rgba(255,255,255,0.25)" 8)
                                      (edges! recovered turn "#7dffa0" 1.5)
                                      (label! "the shape found, turning: green; the truth faint" x2 y2))
                        (= :factored (:stage stage)) (let [S (:S stage) ps (mapv (fn [i] (mapv #(% i) S)) (range (count (first S))))
                                                           m (apply max (map #(apply max (map abs %)) ps))]
                                                       (edges! (mapv #(v/scale % (/ 30.0 m)) ps) turn "#ff8f70" 1.5)
                                                       (label! "factored, not yet metric: right up to a skew" x2 y2))
                        :else (label! "the views so far" x2 y2)))
                (str (cond recovered (str "shape and motion found: what rank three leaves over, " (.toFixed (/ residual (math/sqrt (* 2 (count frames) (count (frames 0))))) 2) " pixels a coordinate"
                                          (when (neg? mirror) "  |  shown mirrored in depth to match: the views cannot tell"))
                           stage (str "stage: " (name (:stage stage)))
                           :else (str (count (frames 0)) " points tracked through " (inc f) " views"))
                     "\nTomasi and Kanade: the tracks, less their means, factored into the camera's turning and the points, rank three")))
            (draw-selective! []
              (let [{{:keys [left truth labels ss parent]} :frame} @model
                    [[x1 y1 w1 h1] [x2 y2 w2 h2]] (panels 2)
                    s (blit! left identity x1 y1 w1 h1)
                    root (fn [l] (loop [l l] (if-let [p (parent l)] (recur p) l)))
                    rect! (fn [[a b c d] color width]
                            (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
                            (.strokeRect ctx (+ x1 (* s a)) (+ y1 (* s b)) (* s (- c a)) (* s (- d b))))]
                (label! "proposals: every region met, at every scale" x1 y1)
                (doseq [b truth] (rect! b "rgba(255,255,255,0.45)" 3))
                (doseq [b (:proposals ss)] (rect! b "rgba(255,212,121,0.6)" 1))
                (when-let [[_ _ n] (:merged ss)] (rect! (get-in ss [:regions n :box]) "#7dffa0" 2.5))
                ;; a small number hashes to itself, so its neighbors would get its color: hash a vector
                (blit! (img/pixels labels root) #(region-colors [:region %]) x2 y2 w2 h2) (label! "the regions now" x2 y2)
                (str (count (:regions ss)) " regions, " (count (distinct (:proposals ss))) " proposals  |  objects proposed to within 0.7 IoU: "
                     (found truth (:proposals ss) 0.7) " of " (count truth)
                     "\nthe two most alike neighbors merged a step -- by color, size, and how well they fill their joint box")))
            (draw-anchors! []
              (let [{{:keys [phase test done train epoch proposals refined eqs total]} :anc} @model
                    [[x1 y1 w1 h1]] (panels 1)
                    {:keys [image truth]} test
                    s (blit! image identity x1 y1 w1 h1)
                    rect! (fn [[a b c d] color width]
                            (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) width)
                            (.strokeRect ctx (+ x1 (* s a)) (+ y1 (* s b)) (* s (- c a)) (* s (- d b))))]
                (doseq [b truth] (rect! b "rgba(255,255,255,0.45)" 3))
                ;; the best eight, by objectness: all of them would hide the scene
                (doseq [b (take 8 proposals)] (rect! b (if refined "rgba(255,143,112,0.6)" "#ffd479") 1.5))
                (doseq [b (take 8 refined)] (rect! b "#7dffa0" 2))
                (str (case phase
                       :examples (str "sampling 64 anchors from training scene " (inc done) " of " (count train))
                       :objectness (str "objectness, epoch " epoch " of 50  |  " (count proposals) " proposals from " (count (anchors-of image)) " anchors"
                                        "  |  objects found to 0.5 IoU: " (found truth proposals 0.5) ", to 0.7: " (found truth proposals 0.7))
                       :regression (str "box regression: " (:n eqs 0) " of " total " examples  |  objects found to 0.7 IoU: " (found truth proposals 0.7))
                       :done (str "refined: objects found to 0.7 IoU " (found truth proposals 0.7) " before, " (found truth refined 0.7) " after"))
                     "\nfaint: the true boxes  yellow, then red: the best eight anchors scored as objects, after non-maximum suppression  green: those moved by box regression")))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (advance!)
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))]
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "method" #js ["edges" "segmentation" "k-means" "optical flow" "stereo" "shading" "texture" "contour"
                                                    "structure from motion" "detection" "selective search" "anchors and regression"
                                                    "boundaries" "normalized cuts" "time to contact" "horizon" "light and color"]) restart!)
        (.onChange (.add gui controls "seed" 1 50 1) restart!)
        (.add gui controls "speed" 1 10 1)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "vision")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
