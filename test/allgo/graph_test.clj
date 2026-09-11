(ns allgo.graph-test
  (:require [allgo.graph :as g]
            [clojure.test :refer [deftest is testing]]))

(deftest disjoint-set-test
  (testing "every node starts alone"
    (let [f (g/disjoint-sets [:a :b :c])]
      (is (not (g/connected? f :a :b)))
      (is (g/connected? f :a :a))))

  (testing "union is transitive"
    (let [f (-> (g/disjoint-sets [:a :b :c :d])
                (g/union :a :b)
                (g/union :b :c))]
      (is (g/connected? f :a :c))
      (is (not (g/connected? f :a :d)))))

  (testing "unioning an already-joined pair changes nothing"
    (let [f (-> (g/disjoint-sets [:a :b]) (g/union :a :b))]
      (is (= f (g/union f :a :b))))))

(deftest components-test
  (is (= #{#{:a :b} #{:c}}
         (g/components [:a :b :c] #{#{:a :b}})))
  (is (= #{#{:a} #{:b} #{:c}}
         (g/components [:a :b :c] #{})))
  (is (g/connected-graph? [:a :b :c] #{#{:a :b} #{:b :c}}))
  (is (not (g/connected-graph? [:a :b :c] #{#{:a :b}}))))

(deftest adjacency-test
  (is (= {:a #{:b} :b #{:a :c} :c #{:b}}
         (g/adjacency #{#{:a :b} #{:b :c}})))
  (testing "isolated nodes appear only when named"
    (is (= #{:a :b :z} (set (keys (g/adjacency [:a :b :z] #{#{:a :b}})))))
    (is (= #{:a :b} (set (keys (g/adjacency #{#{:a :b}})))))))

(deftest complete-graph-test
  (is (= #{#{:a :b} #{:a :c} #{:b :c}} (g/complete-graph [:a :b :c])))
  (is (= 0 (count (g/complete-graph [:a]))))
  ;; n(n-1)/2
  (is (= 45 (count (g/complete-graph (range 10))))))

(deftest minimum-spanning-tree-test
  (testing "a square with one cheap diagonal takes the three cheapest edges"
    ;; a-b 1, b-c 2, c-d 3, d-a 4, a-c 10
    (let [w      {#{:a :b} 1 #{:b :c} 2 #{:c :d} 3 #{:d :a} 4 #{:a :c} 10}
          tree   (g/minimum-spanning-tree [:a :b :c :d] (keys w) w)]
      (is (= 3 (count tree)) "a spanning tree of n nodes has n-1 edges")
      (is (= #{#{:a :b} #{:b :c} #{:c :d}} tree))
      (is (g/connected-graph? [:a :b :c :d] tree))))

  (testing "the tree is minimal, not merely spanning"
    (let [w    {#{:a :b} 1 #{:b :c} 1 #{:a :c} 5}
          tree (g/minimum-spanning-tree [:a :b :c] (keys w) w)]
      (is (= 2 (reduce + (map w tree))))))

  (testing "a disconnected graph yields a spanning forest rather than failing"
    (let [w    {#{:a :b} 1 #{:c :d} 1}
          tree (g/minimum-spanning-tree [:a :b :c :d] (keys w) w)]
      (is (= 2 (count tree)))
      (is (= 2 (count (g/components [:a :b :c :d] tree))))))

  (testing "no edges at all"
    (is (= #{} (g/minimum-spanning-tree [:a :b] #{} (constantly 1)))))

  (testing "a tree over a complete graph still has n-1 edges and is connected"
    (let [nodes (vec (range 12))
          pts   (zipmap nodes (map (fn [i] [(mod i 4) (quot i 4)]) nodes))
          edges (g/complete-graph nodes)
          w     (fn [e] (let [[u v] (vec e)
                              [ux uy] (pts u) [vx vy] (pts v)]
                          (+ (* (- vx ux) (- vx ux)) (* (- vy uy) (- vy uy)))))
          tree  (g/minimum-spanning-tree nodes edges w)]
      (is (= 11 (count tree)))
      (is (g/connected-graph? nodes tree))
      ;; On a 4x3 unit lattice every spanning-tree edge is a unit step.
      (is (= 11 (reduce + (map w tree)))))))

(deftest edges-not-in-test
  (let [edges #{#{:a :b} #{:b :c} #{:a :c}}
        tree  #{#{:a :b} #{:b :c}}]
    (is (= #{#{:a :c}} (g/edges-not-in edges tree)))
    (is (= edges (g/edges-not-in edges #{})))))
