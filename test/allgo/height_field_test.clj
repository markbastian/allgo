(ns allgo.height-field-test
  (:require [allgo.physics.height-field :as hf]
            [clojure.test :refer [deftest is testing]]))

(def ^:private dt (/ 1.0 30.0))

(defn- pond
  ([] (pond {}))
  ([opts] (hf/surface (merge {:size-x 2.5 :size-z 3.0 :spacing 0.04 :depth 0.8} opts))))

(defn- heights [{:keys [nx nz] :as s}]
  (for [i (range nx) j (range nz)] (hf/height s i j)))

(defn- float-ball
  "A ball dropped in and left to settle."
  [density {:keys [radius steps world] :or {radius 0.2 steps 1500 world {}}}]
  (let [s (pond)
        w (atom {:surface s
                 :bodies [(hf/sphere {:pos [0.0 0.9 0.0] :radius radius :density density})]})]
    (dotimes [_ steps] (swap! w hf/step! world))
    @w))

(defn- cap-volume
  "Volume of a sphere of radius `r` submerged to depth `h`."
  [r h]
  (let [h (max 0.0 (min (* 2.0 r) h))]
    (* Math/PI (/ (* h h (- (* 3.0 r) h)) 3.0))))

(defn- submerged-fraction [surface body]
  (let [[x cy z] (:pos body)
        r (:radius body)
        water (hf/sample surface x z)]
    (/ (cap-volume r (- water (- cy r))) (hf/volume body))))

(deftest surface-test
  (let [s (pond {:size-x 1.0 :size-z 2.0 :spacing 0.1 :depth 0.5})]
    (testing "the grid spans the requested size"
      (is (= 11 (:nx s)))
      (is (= 21 (:nz s)))
      (is (= 231 (:n s))))

    (testing "it starts flat at the requested depth"
      (is (every? #(== 0.5 %) (heights s))))

    (testing "the grid is centred on the origin"
      (is (== 0.0 (hf/column-x s (:cx s))))
      (is (== 0.0 (hf/column-z s (:cz s))))
      (is (= [(:cx s) (:cz s)] (hf/nearest-column s 0.0 0.0))))

    (testing "sampling a flat surface gives the depth everywhere, columns or not"
      (is (== 0.5 (hf/sample s 0.0 0.0)))
      (is (== 0.5 (hf/sample s 0.03 -0.07))))

    (testing "sampling interpolates between columns instead of stepping"
      (hf/set-height! s (:cx s) (:cz s) 1.0)
      (let [half (hf/sample s (* 0.5 (:spacing s)) 0.0)]
        (is (< 0.5 half 1.0) "a point between two columns reads between them")))))

(deftest still-water-test
  (testing "flat water stays flat and keeps its volume"
    (let [s (pond)
          v0 (hf/total-volume s)]
      (dotimes [_ 300] (hf/step-surface! s dt))
      (is (every? #(< (abs (- % 0.8)) 1e-5) (heights s)))
      (is (< (abs (- (hf/total-volume s) v0)) 1e-4)))))

(deftest waves-test
  (testing "a splash spreads out and dies down"
    (let [s (pond)
          peak (fn [] (reduce max (heights s)))]
      (hf/splash! s 0.0 0.0 0.2 0.3)
      (let [p0 (peak)]
        (dotimes [_ 30] (hf/step-surface! s dt))
        (let [p1 (peak)]
          (dotimes [_ 300] (hf/step-surface! s dt))
          (is (< p1 p0) "the peak falls as the splash spreads")
          (is (< (peak) p1) "and keeps falling")
          (is (< (abs (- (peak) 0.8)) 0.02) "until the water is nearly flat again")))))

  (testing "a splash moves water without creating or destroying it"
    (let [s (pond)
          v0 (hf/total-volume s)]
      (hf/splash! s 0.0 0.0 0.2 0.3)
      (let [v1 (hf/total-volume s)]
        (is (> v1 v0) "the splash adds the water it drops in")
        (dotimes [_ 300] (hf/step-surface! s dt))
        (is (< (abs (- (hf/total-volume s) v1)) 1e-3)
            "and stepping the surface neither adds nor removes any"))))

  (testing "waves reflect off the edges rather than draining out"
    (let [s (pond)]
      (hf/splash! s 0.0 0.0 0.2 0.3)
      (let [v1 (hf/total-volume s)]
        (dotimes [_ 600] (hf/step-surface! s dt))
        (is (< (abs (- (hf/total-volume s) v1)) 1e-3))))))

(deftest cfl-test
  (testing "the wave speed a grid can carry falls as the step grows"
    (let [s (pond)]
      (is (> (hf/cfl-wave-speed s 0.01) (hf/cfl-wave-speed s 0.1)))))

  (testing "asking for a wave faster than the grid can carry does not blow it up"
    ;; Information cannot cross more than a cell a step; unclamped, the
    ;; stencil answers by oscillating harder every step until it overflows.
    (let [s (pond)]
      (hf/splash! s 0.0 0.0 0.2 0.3)
      (dotimes [_ 400] (hf/step-surface! s dt {:wave-speed 500.0}))
      (is (every? #(Double/isFinite (double %)) (heights s)))
      (is (< (reduce max (heights s)) 2.0)
          "the surface stays in the same order of magnitude it started in")))

  (testing "the clamp does not lower the surface's stored wave speed"
    ;; Clamping in place would mean one big step quietly slowed every step
    ;; after it.
    (let [s (pond)]
      (hf/step-surface! s 1.0 {:wave-speed 5.0})
      (is (= 5.0 (:wave-speed (merge hf/default-surface-world {:wave-speed 5.0})))))))

(deftest buoyancy-test
  (testing "density decides whether a body floats or sinks"
    (let [{:keys [bodies]} (float-ball 0.5 {})
          [_ y _] (:pos (first bodies))]
      (is (> y 0.5) "lighter than water floats near the surface"))
    (let [{:keys [bodies]} (float-ball 5.0 {})
          [_ y _] (:pos (first bodies))]
      (is (< (abs (- y 0.2)) 1e-6) "heavier than water sinks to the floor")))

  (testing "a floating body displaces its own mass of water"
    ;; Archimedes: at rest the submerged fraction is the relative density.
    (doseq [d [0.4 0.6 0.8 0.95]]
      (let [{:keys [surface bodies]} (float-ball d {})
            b (first bodies)
            frac (submerged-fraction surface b)]
        (is (< (abs (double (nth (:vel b) 1))) 1e-3)
            (str "density " d " has come to rest"))
        (is (< (abs (- frac d)) 0.05)
            (str "density " d " floats " (format "%.3f" frac) " submerged")))))

  (testing "where a body floats does not depend on how draggy the water is"
    ;; Damping the buoyant impulse but not the gravity applied after it
    ;; biases the resting depth by 1/damping, which would make Archimedes
    ;; a function of the drag coefficient.
    (let [depth (fn [drag]
                  (let [{:keys [surface bodies]} (float-ball 0.6 {:world {:drag drag}})]
                    (submerged-fraction surface (first bodies))))]
      (is (< (abs (- (depth 8.0) (depth 40.0))) 0.02)
          "the same fraction is submerged whether drag is 8 or 40")))

  (testing "a body pushes up the water it displaces"
    (let [s (pond)
          v0 (hf/total-volume s)
          w (atom {:surface s :bodies [(hf/sphere {:pos [0.0 0.85 0.0] :radius 0.2 :density 0.6})]})]
      (dotimes [_ 60] (swap! w hf/step!))
      (is (> (hf/total-volume s) v0)
          "the water level rises once something is floating in it"))))

(deftest tank-test
  (testing "a body stays inside the tank"
    (let [s (pond)
          tank {:size [2.5 1.0 3.0] :border 0.03}
          b (-> (hf/sphere {:pos [0.0 0.9 0.0] :radius 0.2 :density 3.0})
                (assoc :vel [8.0 0.0 -6.0]))
          w (atom {:surface s :bodies [b]})]
      (dotimes [_ 300] (swap! w hf/step! {:tank tank}))
      (let [[x y z] (:pos (first (:bodies @w)))
            r 0.2
            wx (- (* 0.5 2.5) r (* 0.5 0.03))
            wz (- (* 0.5 3.0) r (* 0.5 0.03))]
        (is (<= (- wx) x wx))
        (is (<= (- wz) z wz))
        (is (>= y (- r 1e-6)))))))

(deftest collision-test
  (testing "overlapping bodies are pushed apart"
    (let [a (hf/sphere {:pos [0.0 0.5 0.0] :radius 0.2 :density 1.0})
          b (hf/sphere {:pos [0.1 0.5 0.0] :radius 0.2 :density 1.0})
          [a' b'] (hf/resolve-collisions [a b])
          [ax _ _] (:pos a') [bx _ _] (:pos b')]
      (is (< (abs (- bx ax)) 0.41))
      (is (>= (abs (- bx ax)) 0.399) "separated to exactly touching")))

  (testing "bodies far apart are left alone"
    (let [a (hf/sphere {:pos [-1.0 0.5 0.0] :radius 0.2 :density 1.0})
          b (hf/sphere {:pos [1.0 0.5 0.0] :radius 0.2 :density 1.0})]
      (is (= [a b] (hf/resolve-collisions [a b])))))

  (testing "a head-on collision exchanges momentum"
    (let [a (-> (hf/sphere {:pos [-0.15 0.5 0.0] :radius 0.2 :density 1.0})
                (assoc :vel [1.0 0.0 0.0]))
          b (-> (hf/sphere {:pos [0.15 0.5 0.0] :radius 0.2 :density 1.0})
                (assoc :vel [-1.0 0.0 0.0]))
          [a' b'] (hf/resolve-collisions [a b])]
      (is (<= (first (:vel a')) 0.0) "the one moving right is turned back")
      (is (>= (first (:vel b')) 0.0) "and so is the one moving left")
      (is (< (abs (+ (* (:mass a') (first (:vel a')))
                     (* (:mass b') (first (:vel b'))))) 1e-9)
          "momentum is conserved"))))

(deftest shape-test
  (testing "a sphere reports the chord of the circle it cuts at each column"
    (let [b (hf/sphere {:pos [0.0 0.5 0.0] :radius 0.2 :density 1.0})]
      (is (< (abs (- (hf/half-height-at b 0.0 0.0) 0.2)) 1e-9)
          "at the centre, the full radius")
      (is (zero? (hf/half-height-at b 0.3 0.0))
          "outside the footprint, nothing")
      (is (< (abs (- (hf/half-height-at b 0.2 0.0) 0.0)) 1e-9)
          "and it tapers to nothing at the rim")))

  (testing "a sphere's footprint bounds it in the plane"
    (let [b (hf/sphere {:pos [1.0 0.5 -2.0] :radius 0.25 :density 1.0})]
      (is (= [0.75 -2.25 1.25 -1.75] (hf/footprint b)))))

  (testing "volume is what it displaces fully submerged"
    (let [b (hf/sphere {:pos [0.0 0.0 0.0] :radius 0.3 :density 1.0})]
      (is (< (abs (- (hf/volume b) (* 4.0 (/ Math/PI 3.0) 0.027))) 1e-9)))))
