(ns allgo.chapra-pde-test
  "Partial differential equations (Chapra and Canale, *Numerical Methods
  for Engineers*, 7th ed., chapters 29-31): the heated plate by Liebmann's
  method with fixed, insulated and irregular edges and its heat flux; the
  heated rod by the explicit, implicit and Crank-Nicolson methods; the
  plate in time by ADI; and finite elements in one and two dimensions."
  (:require [allgo.numerics.fem :as fem]
            [allgo.numerics.linear-systems :as ls]
            [allgo.numerics.pde :as pde]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (every? true? (map #(<= (abs (- %1 %2)) tol) (flatten [a]) (flatten [b]))))

(def ^:private plate {:nx 4 :ny 4 :boundary {:top 100.0 :bottom 0.0 :left 75.0 :right 50.0}})
(defn- interior [T] (for [j [1 2 3] i [1 2 3]] (get-in T [j i])))

(deftest heated-plate
  (testing "Example 29.1: Liebmann's first sweep, over-relaxed by 1.5"
    (let [{:keys [T]} (pde/liebmann (assoc plate :max-iter 1))]
      (is (close? (interior T) [28.125 10.54688 22.70508 38.67188 18.45703 34.18579 80.12696 74.46900 96.99554] 1e-5))))
  (testing "converged, the solution of the nine difference equations, as the book's Fig. 29.5 gives it"
    (let [{:keys [T]} (pde/liebmann plate)
          ;; the nine equations 4 T_ij - neighbors = boundary terms, solved directly
          A (vec (for [j [1 2 3] i [1 2 3]]
                   (vec (for [jj [1 2 3] ii [1 2 3]]
                          (cond (and (= i ii) (= j jj)) 4.0 (= 1 (+ (abs (- i ii)) (abs (- j jj)))) -1.0 :else 0.0)))))
          b (vec (for [j [1 2 3] i [1 2 3]] (+ (if (= i 1) 75.0 0.0) (if (= i 3) 50.0 0.0) (if (= j 3) 100.0 0.0))))
          direct (ls/gauss A b)]
      (is (close? (interior T) direct 1e-8))
      ;; the book's figure is its iterate at a 1% stopping criterion, not converged
      (is (close? (interior T) [42.86 33.26 33.93 63.21 56.25 52.46 78.59 76.06 69.71] 0.1))
      (testing "Example 29.2: the heat flux at (1, 1), 1.856 cal/(cm^2 s) at -56.58 degrees"
        (let [{:keys [qx qy qn theta]} (pde/flux T 1 1 0.49 10.0 10.0)]
          (is (close? [qx qy qn] [1.022 -1.549 1.856] 2e-3))
          (is (close? (math/to-degrees theta) -56.58 0.1))))))
  (testing "Example 29.3: the lower edge insulated"
    (let [{:keys [T]} (pde/liebmann (assoc plate :derivative {:bottom 0.0}))]
      (is (close? (for [j [0 1 2 3] i [1 2 3]] (get-in T [j i]))
                  [71.91 67.01 59.54 72.81 68.31 60.57 76.01 72.84 64.42 83.41 82.63 74.26] 6e-3))))
  (testing "Example 29.4: an irregular lower-left boundary at 75, cutting node (1, 1)'s arms to 0.732"
    (let [{:keys [T]} (pde/liebmann (assoc plate :boundary {:top 100.0 :bottom 75.0 :left 75.0 :right 50.0}
                                           :arms (fn [i j dir] (when (and (= [i j] [1 1]) (#{:left :down} dir)) [0.732 75.0]))))]
      ;; the book's figure has 77.23 at (1, 2), its text 74.23; the equations
      ;; give the figure's. At (1, 1) the book's 74.98 is from its rounded
      ;; coefficients, 0.8453 and 1.1547
      (is (close? (interior T) [74.98 72.76 66.07 77.23 75.00 66.52 83.93 83.48 75.00] 0.025)))))

(def ^:private rod {:k 0.835 :n 5 :dx 2.0 :dt 0.1 :initial 0.0 :left 100.0 :right 50.0})
(defn- inner [T] (subvec T 1 5))

(deftest heated-rod
  (testing "Example 30.1: the explicit method's first two steps"
    (let [traj (pde/heat-1d (assoc rod :t-end 0.2))]
      (is (close? (inner (second (traj 1))) [2.0875 0.0 0.0 1.0438] 1e-4))
      (is (close? (inner (second (traj 2))) [4.0878 0.043577 0.021788 2.0439] 1e-4))))
  (testing "Example 30.2: the simple implicit method at 0.1 s"
    (is (close? (inner (second ((pde/heat-1d (assoc rod :t-end 0.1 :scheme :implicit)) 1))) [2.0047 0.0406 0.0209 1.0023] 1e-4)))
  (testing "Example 30.3: Crank-Nicolson at 0.1 s"
    (is (close? (inner (second ((pde/heat-1d (assoc rod :t-end 0.1 :scheme :crank-nicolson)) 1))) [2.0450 0.0210 0.0107 1.0225] 1e-4)))
  (testing "all three reach the steady straight line between the ends"
    (doseq [s [:explicit :implicit :crank-nicolson]]
      (is (close? (second (peek (pde/heat-1d (assoc rod :t-end 400.0 :scheme s))))
                  [100.0 90.0 80.0 70.0 60.0 50.0] 1e-6) (name s))))
  (testing "an insulated right end: the whole rod comes to the left end's temperature"
    (doseq [s [:explicit :implicit :crank-nicolson]]
      (is (close? (second (peek (pde/heat-1d (assoc rod :t-end 800.0 :scheme s :derivative {:right 0.0}))))
                  (repeat 6 100.0) 1e-4) (name s)))))

(deftest plate-in-time
  (let [adi {:k 0.835 :nx 4 :ny 4 :dx 10.0 :dt 10.0 :boundary {:top 100.0 :bottom 0.0 :left 75.0 :right 50.0}}]
    (testing "Example 30.5: one ADI step of 10 s, the neighbors taken as the book's example takes them"
      (let [[_ T] ((pde/heat-2d (assoc adi :t-end 10.0 :updated-neighbors? true)) 1)]
        (is (close? (interior T) [5.5855 0.4782 3.7388 6.1683 0.8238 4.2359 13.1120 8.3207 11.3606] 1e-4))))
    (testing "ADI and the explicit scheme both settle on Liebmann's steady plate"
      (let [steady (interior (:T (pde/liebmann plate)))]
        (is (close? (interior (second (peek (pde/heat-2d (assoc adi :t-end 3000.0))))) steady 1e-6))
        (is (close? (interior (second (peek (pde/heat-2d (assoc adi :t-end 3000.0 :dt 20.0 :scheme :explicit))))) steady 1e-6))))))

(deftest finite-elements
  (testing "Examples 31.1-31.2: the heated rod, T'' = -10 between 40 and 200, on four elements --
            linear elements exact at the nodes, the end fluxes 66 and -34"
    (let [xs [0.0 2.5 5.0 7.5 10.0]
          {:keys [u left-flux right-flux]} (fem/line-elements xs {:f (constantly 10.0) :left 40.0 :right 200.0})]
      (is (close? u (map #(+ (* -5.0 % %) (* 66.0 %) 40.0) xs) 1e-10))
      (is (close? [left-flux right-flux] [66.0 -34.0] 1e-10))
      (testing "and given the left flux instead of the temperature, the same rod"
        (is (close? (:u (fem/line-elements xs {:f (constantly 10.0) :left-flux 66.0 :right 200.0})) u 1e-9)))))
  (let [mesh (fn [n]
               (let [idx (fn [i j] (+ i (* j (inc n))))
                     nodes (vec (for [j (range (inc n)) i (range (inc n))] [(/ i n) (/ j n)]))
                     elements (vec (mapcat (fn [[i j]] [[(idx i j) (idx (inc i) j) (idx (inc i) (inc j))]
                                                        [(idx i j) (idx (inc i) (inc j)) (idx i (inc j))]])
                                           (for [j (range n) i (range n)] [i j])))
                     edge (for [j (range (inc n)) i (range (inc n)) :when (or (#{0 n} i) (#{0 n} j))] (idx i j))]
                 [nodes elements edge]))]
    (testing "triangles: a linear field is reproduced exactly (the patch test)"
      (let [[nodes elements edge] (mesh 4)
            g (fn [[x y]] (+ 1.0 (* 2.0 x) (* 3.0 y)))
            u (fem/triangles nodes elements {:fixed (into {} (map (fn [k] [k (g (nodes k))]) edge))})]
        (is (close? u (map g nodes) 1e-12))))
    (testing "and the Poisson problem -lap u = 2 pi^2 sin(pi x) sin(pi y) converges at second order"
      (let [err (fn [n]
                  (let [[nodes elements edge] (mesh n)
                        exact (fn [[x y]] (* (math/sin (* math/PI x)) (math/sin (* math/PI y))))
                        u (fem/triangles nodes elements {:f (fn [x y] (* 2.0 math/PI math/PI (exact [x y])))
                                                         :fixed (into {} (map (fn [k] [k 0.0]) edge))})]
                    (apply max (map #(abs (- %1 (exact %2))) u nodes))))
            e8 (err 8) e16 (err 16)]
        (is (< e16 0.01))
        (is (< 3.0 (/ e8 e16) 5.0))))))
