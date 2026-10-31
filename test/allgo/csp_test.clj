(ns allgo.csp-test
  (:require [allgo.random :as random]
            [allgo.search.csp :as csp]
            [allgo.search.problems :as probs]
            [clojure.test :refer [deftest is testing]]))

;; ---------------------------------------------------------------------------
;; The worked problems of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 6

(def ^:private australia
  "The chapter's first example: Australia's states and territories, neighbors in different
  colors."
  (csp/csp {:variables [:wa :nt :q :nsw :v :sa :t]
            :domains [:red :green :blue]
            :edges [[:wa :nt] [:wa :sa] [:nt :sa] [:nt :q] [:sa :q] [:sa :nsw] [:sa :v] [:q :nsw] [:nsw :v]]
            :constraint (fn [_ x _ y] (not= x y))}))

(def ^:private two-two-four
  "The chapter's cryptarithmetic: TWO + TWO = FOUR, a different digit for each letter, and
  the carries C1 C2 C3 as variables of their own -- constraints of three
  and five variables, and an Alldiff written out as binary inequalities."
  (let [letters [:t :w :o :f :u :r]]
    (csp/csp {:variables (into letters [:c1 :c2 :c3])
              :domains (merge (zipmap letters (repeat (vec (range 10))))
                              {:t (vec (range 1 10)) :f (vec (range 1 10))
                               :c1 [0 1] :c2 [0 1] :c3 [0 1]})
              :edges (for [a letters b letters :when (neg? (compare a b))] [a b])
              :constraint (fn [_ x _ y] (not= x y))
              :constraints [{:scope [:o :r :c1] :ok? (fn [[o r c1]] (= (+ o o) (+ r (* 10 c1))))}
                            {:scope [:c1 :w :u :c2] :ok? (fn [[c1 w u c2]] (= (+ c1 w w) (+ u (* 10 c2))))}
                            {:scope [:c2 :t :o :c3] :ok? (fn [[c2 t o c3]] (= (+ c2 t t) (+ o (* 10 c3))))}
                            {:scope [:c3 :f] :ok? (fn [[c3 f]] (= c3 f))}]})))

(deftest australia-test
  (testing "three colors do it, eighteen ways: six for the mainland, three for Tasmania"
    (is (csp/solution? australia (csp/backtracking australia)))
    (is (= 18 (count (csp/all-solutions australia))))
    (is (every? #(csp/solution? australia %) (csp/all-solutions australia))))
  (testing "two do not -- which arc consistency cannot see: the chapter's case for path consistency"
    (let [two (assoc australia :domains (zipmap (:variables australia) (repeat [:red :green])))]
      (is (nil? (csp/backtracking two)))
      (is (= (:domains two) (csp/ac-3 two)))))
  (testing "MRV with the degree heuristic starts with SA, in five constraints"
    (let [first-frame (second (csp/backtracking-seq australia {:select :mrv}))]
      (is (= :sa (key (first (:assignment (peek (:stack first-frame)))))))))
  (testing "the chapter's forward-checking example: after WA = red, Q = green leaves NT and SA only blue, and V = blue empties SA"
    (let [d0 (:domains australia)
          [d1] (csp/inference australia :forward-checking {} d0 :wa :red)
          [d2] (csp/inference australia :forward-checking {:wa :red} d1 :q :green)]
      (is (= [:green :blue] (d1 :nt) (d1 :sa)))
      (is (= [:blue] (d2 :nt) (d2 :sa)))
      (is (= [:red :blue] (d2 :nsw)))
      (is (nil? (csp/inference australia :forward-checking {:wa :red :q :green} d2 :v :blue)))))
  (testing "removing SA leaves a tree: the chapter's cycle cutset"
    (is (= [:sa] (csp/cycle-cutset australia)))
    (is (csp/solution? australia (csp/cutset-conditioning australia)))
    (is (not (csp/tree? australia)))))

(deftest cryptarithmetic-test
  (testing "TWO + TWO = FOUR: every solution backtracking finds adds up, and it finds every one"
    (let [sols (csp/all-solutions two-two-four)
          sums (set (map (fn [{:keys [t w o]}] (+ (* 100 t) (* 10 w) o)) sols))
          brute (set (for [t (range 1 10) w (range 10) o (range 10)
                           :let [two (+ (* 100 t) (* 10 w) o) four (* 2 two)
                                 digits (fn [n k] (mapv #(mod (quot n (long (Math/pow 10 %))) 10) (range (dec k) -1 -1)))
                                 [f o' u r] (digits four 4)]
                           :when (and (<= 1000 four 9999) (= o o')
                                      (= 6 (count (distinct [t w o f u r]))))]
                       two))]
      (is (seq sols))
      (is (= brute sums))
      (is (every? #(csp/solution? two-two-four %) sols)))))

;; ---------------------------------------------------------------------------
;; Sudoku, queens, maps

(deftest sudoku-test
  (let [[grid solution] (probs/sudoku-puzzle (random/rng 4) 55)
        p (probs/sudoku grid)]
    (testing "the generator's puzzle has its solution, and only it"
      (is (= 1 (count (csp/all-solutions p))))
      (is (= solution (mapv (csp/backtracking p {:select :mrv :inference :mac}) (range 81)))))
    (testing "AC-3 removes no value a solution needs, and fixes some cells outright"
      (let [d (csp/ac-3 p)]
        (is (every? (fn [k] (some #{(solution k)} (d k))) (range 81)))
        (is (> (count (filter #(= 1 (count (d %))) (range 81)))
               (count (remove zero? grid))))))
    (testing "every combination of the heuristics solves it, and inference saves assignments"
      ;; Not plain backtracking: on 55 blanks it ran past ten minutes.
      (let [runs (for [opts [{:select :mrv :inference :forward-checking} {:select :mrv :inference :mac}
                             {:select :mrv :order :lcv :inference :mac}]]
                   (let [end (last (csp/backtracking-seq p opts))]
                     [(= solution (mapv (:value end) (range 81))) (:assignments (:context end))]))]
        (is (every? first runs))
        (is (< (second (nth runs 1)) (second (first runs))))))))

(deftest queens-test
  (testing "backtracking: forward checking and MAC cut the assignments tried"
    (let [q (probs/queens-csp 10)
          tries (fn [opts] (:assignments (:context (last (csp/backtracking-seq q opts)))))]
      (is (csp/solution? q (csp/backtracking q {:select :mrv :inference :forward-checking})))
      (is (> (tries {}) (tries {:select :mrv :inference :forward-checking}) (tries {:select :mrv :inference :mac})))))
  (testing "min-conflicts solves 200 queens in a few hundred steps at most"
    (let [q (probs/queens-csp 200)
          end (last (csp/min-conflicts-seq q {:rng (random/rng 2)}))]
      (is (csp/solution? q (:current end)))
      (is (< (:steps end) 400)))))

(deftest tree-and-cutset-test
  (testing "the tree solver agrees with backtracking on random trees, with no backtracking of its own"
    (doseq [k (range 1 9)]
      (let [rng (random/rng k)
            n 30
            edges (for [i (range 1 n)] [(random/below rng i) i])
            ;; A random relation per edge, so that some trees have no solution.
            allowed (into {} (for [e edges] [e (set (for [x (range 3) y (range 3) :when (< (double (rng)) 0.6)] [x y]))]))
            p (csp/csp {:variables (range n) :domains (range 3) :edges edges
                        :constraint (fn [a x b y] (if (< a b) (contains? (allowed [a b]) [x y]) (contains? (allowed [b a]) [y x])))})]
        (is (csp/tree? p))
        (is (= (some? (csp/backtracking p)) (some? (csp/tree-solver p))) k)
        (when-let [s (csp/tree-solver p)] (is (csp/solution? p s))))))
  (testing "cutset conditioning colors a random map of thirty regions in four colors -- its cutset of 11 is cheap; at 60 regions the cutset is 26, and 300,000 steps were not enough"
    (let [{p :csp} (probs/region-map 30 4 2)]
      (is (not (csp/tree? p)))
      (is (csp/solution? p (csp/cutset-conditioning p)))
      (is (csp/solution? p (csp/min-conflicts p))))))
