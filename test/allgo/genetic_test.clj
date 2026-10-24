(ns allgo.genetic-test
  (:require [allgo.procedural.genetic :as g]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

(defn- valid? [expr]
  (cond
    (number? expr) true
    (keyword? expr) (contains? #{:x :y :z} expr)
    :else (let [{:keys [arity]} (g/ops (first expr))]
            (and arity (= arity (dec (count expr))) (every? valid? (rest expr))))))

(defn- finite-everywhere? [genome]
  (let [t (g/texture genome)]
    (every? (fn [[x y]] (every? #(and (Double/isFinite %) (<= 0.0 % 1.0)) (t x y 0.3)))
            (for [i (range 7) j (range 7)] [(- (* 0.6 i) 1.8) (- (* 0.6 j) 1.8)]))))

(deftest expression-test
  (testing "an expression evaluates as it reads"
    (is (= 7.0 ((g/compile-expr [:add [:mul 2 :x] :y]) 3 1 0)))
    (is (= 1.0 ((g/compile-expr [:max :z [:abs :x]]) -1 0 0.5))))
  (testing "it prints as it reads"
    (is (= "(add (mul 2.5 x) (noise y z 0.33))"
           (g/expr->string [:add [:mul 2.5 :x] [:noise :y :z 0.333]]))))
  (testing "counting nodes"
    (is (= 5 (g/node-count [:add [:mul 2 :x] :y])))))

(deftest population-test
  (let [pop (g/population (random/rng 3) 12)]
    (testing "random genomes are well formed, drawable and varied"
      (is (= 12 (count pop)))
      (doseq [genome pop]
        (is (valid? (:expr genome)))
        (is (<= 5 (g/node-count (:expr genome)) g/max-nodes))
        (is (= 5 (count (:palette genome))))
        (is (finite-everywhere? genome)))
      (is (= 12 (count (distinct pop)))))
    (testing "a seed repeats its population"
      (is (= pop (g/population (random/rng 3) 12))))))

(deftest breeding-test
  (let [rng (random/rng 9)
        [a b] (g/population rng 2)]
    (testing "a mutant is well formed and differs from its parent"
      (dotimes [_ 30]
        (let [m (g/mutate rng a 0.2)]
          (is (valid? (:expr m)))
          (is (not= m a))
          (is (finite-everywhere? m)))))
    (testing "a child's tree is built from its parents' subtrees"
      (let [subtrees (fn subtrees [e] (cons e (when (vector? e) (mapcat subtrees (rest e)))))
            from-b (set (subtrees (:expr b)))]
        (dotimes [_ 30]
          (let [c (g/cross rng a b)]
            (is (valid? (:expr c)))
            ;; Some subtree of the child came whole from b.
            (is (some from-b (subtrees (:expr c))))
            (is (every? (set (concat (:palette a) (:palette b))) (:palette c)))))))
    (testing "the chosen survive into the next generation, first"
      (is (= [a] (take 1 (g/breed rng [a] 12 0.2))))
      (is (= [a b] (take 2 (g/breed rng [a b] 12 0.2))))
      (is (= 12 (count (g/breed rng [a b] 12 0.2))))
      (is (= 12 (count (g/breed rng [] 12 0.2)))))))
