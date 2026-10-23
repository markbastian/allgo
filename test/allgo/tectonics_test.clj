(ns allgo.tectonics-test
  (:require [allgo.procedural.tectonics :as tc]
            [clojure.test :refer [deftest is testing]]))

(def ^:private world (tc/plates {:resolution 4 :plates 10 :seed 7}))

(defn- connected?
  "Whether the cells in `members` form one piece over `neighbors`."
  [neighbors members]
  (let [members (set members)
        start (first members)]
    (loop [seen #{start} frontier [start]]
      (if-let [c (peek frontier)]
        (let [fresh (remove seen (filter members (nth neighbors c)))]
          (recur (into seen fresh) (into (pop frontier) fresh)))
        (= seen members)))))

(deftest plates-test
  (testing "the same seed gives the same world"
    (is (= (:elevation world) (:elevation (tc/plates {:resolution 4 :plates 10 :seed 7})))))
  (testing "every cell is on one of the plates, and each plate is one piece"
    (is (every? #(< -1 % 10) (:plate world)))
    (doseq [[_ cells] (group-by (:plate world) (range (count (:plate world))))]
      (is (connected? (:neighbors world) cells))))
  (testing "the share of continental crust is what was asked"
    (let [c (:continental world)]
      (is (< (abs (- 0.35 (/ (count (filter true? c)) (double (count c))))) 0.01)))))

(deftest boundary-kind-test
  (let [a [1.0 0.0 0.0] b [0.99 0.1 0.0]]
    (testing "plates closing on each other converge, and how hard is how fast"
      (let [[kind s] (tc/boundary-kind a b [0.0 1.0 0.0] [0.0 -1.0 0.0])]
        (is (= :converging kind))
        (is (< 1.9 s 2.0))))
    (testing "plates parting diverge"
      (is (= :diverging (first (tc/boundary-kind a b [0.0 -1.0 0.0] [0.0 1.0 0.0])))))
    (testing "plates sliding past build nothing"
      (is (= :transform (first (tc/boundary-kind a b [0.0 0.0 1.0] [0.0 0.0 -1.0])))))))

(deftest ground-test
  (let [{:keys [elevation continental neighbors plate]} world
        n (count elevation)
        mean (fn [xs] (/ (reduce + xs) (double (count xs))))
        boundary? (fn [c] (some #(not= (plate c) (plate %)) (nth neighbors c)))]
    (testing "continents stand above the sea floor"
      (is (> (mean (map elevation (filter continental (range n))))
             (+ 0.4 (mean (map elevation (remove continental (range n))))))))
    (testing "the highest ground is at plate boundaries, not in the middle of plates"
      (let [top (take (quot n 100) (sort-by #(- (elevation %)) (range n)))
            near (fn [c] (or (boundary? c) (some boundary? (nth neighbors c))))]
        (is (> (count (filter near top)) (* 0.8 (count top))))))))

(deftest terrain-test
  (let [locate (tc/locator world)
        bare (tc/terrain world {:detail 0.0})
        {:keys [points elevation]} world]
    (testing "barycentric weights are a partition of one, none negative"
      (doseq [i (range 0 (count points) 37)]
        (let [p (mapv #(* 0.5 %) (map + (points i) (points (first ((:neighbors world) i)))))
              l (Math/sqrt (reduce + (map * p p)))
              [_ ws] (locate (mapv #(/ % l) p))]
          (is (< (abs (- 1.0 (reduce + ws))) 1e-9))
          (is (every? #(> % -1e-9) ws)))))
    (testing "without detail, the terrain at a cell is that cell's elevation"
      (doseq [i (range 0 (count points) 53)]
        (let [[x y z] (points i)]
          (is (< (abs (- (elevation i) (bare x y z))) 1e-9)))))))
