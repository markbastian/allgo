(ns allgo.decision.games-test
  (:require [allgo.decision.games :as games]
            [clojure.test :refer [deftest is testing]]))

;; Multiagent decision making, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, chapter 18. Every
;; expected answer here is worked out by hand in its comment.

(defn- close? [a b] (every? #(< (abs %) 1e-9) (map - (flatten [a]) (flatten [b]))))

(def ^:private dilemma
  ;; refuse, testify: years in prison, negated
  (games/bimatrix [:refuse :testify] [:refuse :testify] [[[-1 -1] [-10 0]] [[0 -10] [-5 -5]]]))

(deftest dilemma-test
  (testing "testifying dominates, so both testify -- the one equilibrium, and not Pareto optimal"
    (is (= [1 1] (games/dominant-equilibrium dilemma)))
    (is (= [[1 1]] (games/pure-nash dilemma)))
    (is (not (games/pareto-optimal? dilemma [1 1])))
    (is (games/pareto-optimal? dilemma [0 0]))
    (is (= [[0 0]] (games/utilitarian dilemma)) "refusing together is best for the two of them")))

(deftest elimination-test
  (testing "Right is dominated by Middle; then Down by Up; then Left by Middle: (Up, Middle)"
    (let [g (games/bimatrix [:up :down] [:left :middle :right] [[[1 0] [1 2] [0 1]] [[0 3] [0 1] [2 0]]])
          steps (games/elimination-seq g)]
      (is (= [[1 2] [0 1] [1 0]] (keep :removed steps)))
      (is (= [[0] [1]] (:alive (last steps)))))))

(deftest zero-sum-test
  (testing "two-finger Morra: even indifferent when 2p - 3(1 - p) = -3p + 4(1 - p), p = 7/12, value -1/12"
    (let [morra (games/zero-sum [:one :two] [:one :two] [[2 -3] [-3 4]])
          {:keys [strategy value]} (games/maximin morra)]
      (is (close? [7/12 5/12] strategy))
      (is (close? -1/12 value))
      (is (close? [[[7/12 5/12] [7/12 5/12]]] (games/mixed-nash morra)))))
  (testing "rock, paper, scissors: a third each, value 0"
    (let [rps (games/zero-sum [:r :p :s] [:r :p :s] [[0 -1 1] [1 0 -1] [-1 1 0]])]
      (is (close? [1/3 1/3 1/3] (:strategy (games/maximin rps))))
      (is (close? 0 (:value (games/maximin rps)))))))

(deftest mixed-nash-test
  (testing "battle of the sexes: two pure equilibria and a mixed one, each going to their favorite 2/3 of the time"
    (let [bos (games/bimatrix [:opera :football] [:opera :football] [[[2 1] [0 0]] [[0 0] [1 2]]])
          eqs (games/mixed-nash bos)]
      (is (= 3 (count eqs)))
      (is (some #(close? [[2/3 1/3] [1/3 2/3]] %) eqs))
      (testing "and no player gains by deviating from any of them"
        (doseq [[p q] eqs i [0 1] a [0 1]]
          (let [pure [(if (= a 0) 1.0 0.0) (if (= a 0) 0.0 1.0)]
                u (games/expected-utility bos i p q)
                dev (if (= i 0) (games/expected-utility bos 0 pure q) (games/expected-utility bos 1 p pure))]
            (is (<= dev (+ u 1e-9)))))))))

(deftest repeated-test
  (testing "tit-for-tat loses only the first round to always-defect"
    (is (= [-55 -45] (:totals (nth (games/play-seq games/prisoners-dilemma games/tit-for-tat games/always-defect) 10)))))
  (testing "the reciprocators do best in a round robin"
    (let [scores (games/tournament games/prisoners-dilemma
                                   {:alld games/always-defect :allc games/always-cooperate
                                    :tft games/tit-for-tat :grim games/grim} 10)]
      (is (= (scores :tft) (scores :grim)))
      (is (> (scores :tft) (scores :allc)))
      (is (> (scores :tft) (scores :alld))))))

(deftest bargaining-test
  (testing "backward induction alternates toward Rubinstein's split, (1 - d2)/(1 - d1 d2)"
    (is (close? [1.0 0.1 0.91] (take 3 (games/bargaining-seq 0.9 0.9))))
    (is (close? (first (games/rubinstein 0.9 0.8)) (nth (games/bargaining-seq 0.9 0.8) 400)))
    (is (close? [(/ 1 1.9) (/ 0.9 1.9)] (games/rubinstein 0.9 0.9)))))

(deftest shapley-test
  (testing "the glove game, one left glove and two right: a pair is worth 1; left gets 2/3, each right 1/6"
    (let [v (fn [c] (min (count (filter #{:l} c)) (count (filter #{:r1 :r2} c))))]
      (is (= {:l 2/3 :r1 1/6 :r2 1/6} (games/shapley [:l :r1 :r2] v)))))
  (testing "the values always share out the grand coalition's worth"
    (let [v (fn [c] (* (count c) (count c)))]
      (is (== 16 (reduce + (vals (games/shapley [1 2 3 4] v))))))))

(deftest voting-test
  (let [ballots (concat (repeat 4 [:a :b :c]) (repeat 3 [:b :c :a]) (repeat 2 [:c :b :a]))]
    (testing "four rank a b c, three b c a, two c b a: plurality picks a; Borda b (12 to a's 8 and c's 7); b beats both head to head; the runoff drops c and b wins"
      (is (= #{:a} (games/plurality ballots)))
      (is (= {:a 8 :b 12 :c 7} (games/borda ballots)))
      (is (= :b (games/condorcet-winner ballots)))
      (is (= [#{:a :b :c} #{:a :b}] (games/runoff-seq ballots)))
      (is (= :b (games/runoff-winner ballots)))))
  (testing "a cycle of preferences has no Condorcet winner"
    (is (nil? (games/condorcet-winner [[:a :b :c] [:b :c :a] [:c :a :b]])))))

(deftest auctions-test
  (let [values [10 30 25 7]]
    (testing "the English auction ends one raise past the second-highest value"
      (is (= {:price 26 :in #{1}} (last (games/english-seq values 0 1)))))
    (testing "first price: the winner pays its bid; second price: the runner-up's"
      (is (= {:winner 1 :price 30} (games/sealed-bid values :first)))
      (is (= {:winner 1 :price 25} (games/sealed-bid values :second))))
    (testing "in a second-price auction no bid does better than the truth, whatever the others bid"
      (doseq [others [[10 25 7] [10 35 7] [0 0 0]]
              bid (range 0 50 5)]
        (let [truthful (games/auction-utility (into [30] others) 0 30 :second)
              other (games/auction-utility (into [bid] others) 0 30 :second)]
          (is (>= truthful other))))))
  (testing "VCG for two items: the two highest win, each paying the third-highest bid"
    (is (= {:winners [1 2] :prices {1 18 2 18}} (games/vcg [10 30 25 7 18] 2)))))
