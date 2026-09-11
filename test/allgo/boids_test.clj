(ns allgo.boids-test
  (:require [allgo.geometry.gjk :as gjk]
            [allgo.simulation.boids :as b]
            [clojure.test :refer [deftest is testing]]))

(defn- seeded-flock
  "Build a flock from a fixed seed. `b/flock` draws on `rand`, which would
  make the alignment threshold below flaky from run to run."
  [seed n bounds max-speed]
  (let [g (java.util.Random. seed)]
    (vec (repeatedly n
                     #(let [dir (mapv (fn [_] (- (.nextDouble g) 0.5)) bounds)]
                        {:pos (mapv (fn [hi] (* hi (.nextDouble g))) bounds)
                         :vel (b/with-magnitude dir (* max-speed (+ 0.5 (* 0.5 (.nextDouble g)))))})))))

(defn- order
  "Polarisation: 0 when headings are scattered, 1 when the flock is aligned."
  [flock]
  (b/mag (b/v* (reduce b/v+ (map (comp b/normalize :vel) flock)) (/ 1.0 (count flock)))))

(deftest flocking-emerges
  (testing "an initially scattered flock aligns"
    (let [bounds [400 400]
          f0     (seeded-flock 20260910 60 bounds (:max-speed b/defaults))
          run    (iterate #(b/step % bounds) f0)
          start  (order f0)
          end    (order (nth run 600))]
      (is (< start 0.4) (str "starts disordered, was " start))
      (is (> end 0.6) (str "ends aligned, was " end))
      (is (> end (* 2 start)) "alignment increases substantially"))))

(deftest flock-stays-well-formed
  (let [bounds [400 400]
        final  (nth (iterate #(b/step % bounds) (seeded-flock 1 40 bounds (:max-speed b/defaults))) 200)]
    (is (= 40 (count final)))
    (testing "wrapping keeps every boid inside the world"
      (doseq [{[x y] :pos} final]
        (is (and (<= 0 x 400) (<= 0 y 400)))))
    (testing "speed never exceeds the limit"
      (doseq [{:keys [vel]} final]
        (is (<= (b/mag vel) (+ (:max-speed b/defaults) 1e-9)))))
    (testing "no boid goes non-finite"
      (doseq [{:keys [pos vel]} final]
        (is (every? #(Double/isFinite (double %)) (concat pos vel)))))))

(deftest rules-are-dimension-agnostic
  (testing "the same rules drive a 3D flock"
    (let [bounds [100 100 100]
          params {:perception-radius 22 :separation-radius 9 :max-speed 0.8 :max-force 0.02}
          f0     (seeded-flock 2 40 bounds (:max-speed params))
          final  (nth (iterate #(b/step % bounds params) f0) 400)]
      (is (= 3 (count (:pos (first f0)))))
      (is (> (order final) (order f0)) "alignment increases in 3D too")
      (doseq [{:keys [pos]} final]
        (is (every? #(<= 0 % 100) pos))))))

(deftest step-does-not-depend-on-update-order
  (testing "every boid sees the same previous state, so reversing the flock only permutes the result"
    ;; Compared within a tolerance, not bitwise. The invariant is real, but
    ;; reversing the flock also reverses the order neighbour contributions
    ;; are summed in, and floating-point addition is not associative -- the
    ;; two answers can differ in the last unit in the last place.
    (let [bounds [400 400]
          f      (seeded-flock 4242 30 bounds (:max-speed b/defaults))
          forward (sort-by :pos (b/step f bounds))
          reverse* (sort-by :pos (b/step (vec (reverse f)) bounds))]
      (is (= (count forward) (count reverse*)))
      (doseq [[x y] (map vector forward reverse*)]
        (is (< (b/mag (b/v- (:pos x) (:pos y))) 1e-9) "same positions")
        (is (< (b/mag (b/v- (:vel x) (:vel y))) 1e-9) "same velocities")))))

(defn- lift [pos]
  (if (= 3 (count pos)) (vec pos) [(nth pos 0) (nth pos 1) 0.0]))

(defn- penetrating?
  "Allow a hair of slack: contact resolution puts a boid exactly on the
  surface, and asking whether a sphere of its full radius touches would then
  be true by construction."
  [{:keys [pos]} obstacle radius]
  (gjk/intersects? (gjk/sphere (lift pos) (* 0.98 radius)) (:support obstacle)))

(deftest boids-do-not-enter-obstacles
  (testing "in the plane"
    (let [bounds  [400 400]
          radius  3.0
          obs     [(b/sphere-obstacle [200 200] 55)
                   (b/box-obstacle [40 300] [140 360])
                   (b/sphere-obstacle [320 90] 40)]
          params  {:edges :bounce :obstacles obs :boid-radius radius}
          run     (iterate #(b/step % bounds params)
                           (seeded-flock 101 60 bounds (:max-speed b/defaults)))]
      (doseq [t [50 200 600]]
        (doseq [boid (nth run t), o obs]
          (is (not (penetrating? boid o radius)) (str "tick " t))))))
  (testing "in three dimensions"
    (let [world  [100 100 100]
          radius 1.5
          obs    [(b/sphere-obstacle [50 50 50] 18)
                  (b/box-obstacle [10 10 60] [40 40 90])]
          params {:edges :bounce :obstacles obs :boid-radius radius
                  :perception-radius 30 :separation-radius 12
                  :max-speed 0.7 :max-force 0.02 :avoid-radius 14 :avoid-weight 3.0}
          run    (iterate #(b/step % world params)
                          (seeded-flock 202 40 world 0.7))]
      (doseq [t [50 200 600]]
        (doseq [boid (nth run t), o obs]
          (is (not (penetrating? boid o radius)) (str "tick " t)))))))

(deftest bouncing-keeps-the-flock-inside-without-wrapping
  (let [bounds [400 400]
        params {:edges :bounce}
        run    (iterate #(b/step % bounds params)
                        (seeded-flock 7 40 bounds (:max-speed b/defaults)))]
    (doseq [t [10 100 400]]
      (doseq [{[x y] :pos} (nth run t)]
        (is (and (<= 0 x 400) (<= 0 y 400)) (str "tick " t))))
    (testing "a bounced boid reverses rather than teleporting across"
      (let [[p v] (b/bounce [-3.0 200.0] [-2.0 1.0] [400 400])]
        (is (= [3.0 200.0] p) "mirrored back inside, not wrapped to the far edge")
        (is (= [2.0 1.0] v) "and heading back in")))))

(deftest wrapping-remains-the-default
  (let [bounds [400 400]
        far    (b/step [{:pos [399.0 200.0] :vel [5.0 0.0]}] bounds)]
    (is (< (first (:pos (first far))) 100)
        "with no :edges given a boid leaving the right edge reappears on the left")))

(deftest obstacle-constructors
  (testing "a planar box becomes a prism with depth, since GJK needs volume"
    (let [{:keys [support radius]} (b/box-obstacle [0 0] [10 10])]
      (is (pos? radius))
      (is (not (gjk/intersects? (gjk/sphere [20.0 5.0 0.0] 1.0) support)) "clear of it")
      (is (gjk/intersects? (gjk/sphere [5.0 5.0 0.0] 1.0) support) "inside it")))
  (testing "a sphere obstacle meets the plane as its equator"
    (let [{:keys [support]} (b/sphere-obstacle [0 0] 10)]
      (is (gjk/intersects? (gjk/sphere [9.0 0.0 0.0] 0.5) support))
      (is (not (gjk/intersects? (gjk/sphere [11.0 0.0 0.0] 0.5) support)))))
  (testing "the bounding sphere used for broad phase really encloses the body"
    (doseq [{:keys [support centre radius]} [(b/sphere-obstacle [1 2 3] 5)
                                             (b/box-obstacle [0 0 0] [4 6 8])]]
      (is (not (gjk/intersects? (gjk/sphere centre (* 0.999 radius))
                                (gjk/translate support [(* 3 radius) 0.0 0.0])))))))

(deftest flat-box-bounding-sphere-measures-the-rectangle
  ;; A flat box is extruded into a prism so GJK's simplex never lands on a
  ;; degenerate face, but that depth is invented and the prism straddles the
  ;; flock's plane. Letting it into the radius inflates the sphere by more
  ;; than 2x and the broad-phase cull stops culling.
  (let [{:keys [centre radius]} (b/box-obstacle [0 0] [6 6])]
    (is (= 0.0 (nth centre 2)) "the prism is centred on the flock's plane")
    (is (< (Math/abs (- radius (* 0.5 (Math/hypot 6 6)))) 1e-9)
        "radius is the half-diagonal of the rectangle, not of the prism"))

  (testing "it still encloses the rectangle the flock actually meets"
    (let [{:keys [centre radius]} (b/box-obstacle [2 3] [10 9])
          [cx cy] centre]
      (doseq [corner [[2 3] [10 3] [2 9] [10 9]]]
        (is (<= (Math/hypot (- (first corner) cx) (- (second corner) cy))
                (+ radius 1e-9)))))))

(deftest obstacle-index-agrees-with-the-linear-scan
  (let [obstacles (vec (for [x (range 0 200 20) y (range 0 200 20)]
                         (b/box-obstacle [x y] [(+ x 8) (+ y 8)])))
        base      (merge b/defaults {:obstacles obstacles :avoid-radius 12
                                     :boid-radius 1.5 :max-speed 2.0})
        indexed   (assoc base :obstacle-index
                         (b/index-obstacles obstacles (+ 12 1.5)))]
    (testing "avoidance is unchanged by the presence of an index"
      (doseq [boid [{:pos [37.0 41.0] :vel [1.0 0.0]}
                    {:pos [4.0 4.0] :vel [0.0 1.0]}
                    {:pos [1000.0 1000.0] :vel [1.0 1.0]}
                    {:pos [101.0 99.0] :vel [-1.0 0.5]}]]
        (is (= (b/avoidance boid base) (b/avoidance boid indexed))
            (str "at " (:pos boid)))))

    (testing "a boid far from everything is steered by nothing either way"
      (let [far {:pos [5000.0 5000.0] :vel [1.0 0.0]}]
        (is (= [0.0 0.0] (b/avoidance far indexed)))))

    (testing "an obstacle wider than a cell is filed in each cell it covers"
      (let [wide (b/box-obstacle [0 0] [100 100])
            idx  (b/index-obstacles [wide] 5.0)]
        (is (< 1 (count (:cells idx))) "it spans several cells")
        (is (every? #(= 1 (count %)) (vals (:cells idx)))
            "but only once within any one of them")
        ;; Lookup still has to drop duplicates across the 3x3 it reads, and
        ;; that it does is what the avoidance comparison above establishes:
        ;; a doubly-counted obstacle would steer twice as hard.
        (let [boid {:pos [50.0 50.0] :vel [1.0 0.0]}
              base (merge b/defaults {:obstacles [wide] :avoid-radius 12 :boid-radius 1.5})]
          (is (= (b/avoidance boid base)
                 (b/avoidance boid (assoc base :obstacle-index idx)))))))))

(deftest stepping-with-obstacles-is-index-independent
  ;; `step` builds an index internally; the flock it produces must match
  ;; what the unindexed path produces, or the optimisation changed the sim.
  (let [obstacles (vec (for [x (range 0 300 30)] (b/box-obstacle [x 100] [(+ x 12) 112])))
        params    {:obstacles obstacles :avoid-radius 15 :boid-radius 2.0
                   :edges :bounce :perception-radius 40}
        flock     (vec (for [i (range 25)]
                         {:pos [(* 11.0 i) (+ 60.0 (* 3.0 (mod i 7)))]
                          :vel [1.0 0.4]}))
        indexed   (nth (iterate #(b/step % [300.0 300.0] params) flock) 12)
        unindexed (nth (iterate (fn [f]
                                  (let [idx (b/index-flock f 40)
                                        p   (merge b/defaults params)]
                                    (mapv #(b/step-boid % idx [300.0 300.0] p) f)))
                                flock)
                       12)]
    (is (= indexed unindexed))))
