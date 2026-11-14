(ns allgo.search.adversarial-test
  (:require [allgo.random :as random]
            [allgo.search.adversarial :as adv]
            [allgo.search.games :as games]
            [clojure.test :refer [deftest is testing]]))

;; ---------------------------------------------------------------------------
;; The worked problems of Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 5

(def ^:private figure-5-2
  "The chapter's two-ply game: MAX moves a1, a2 or a3, then MIN b, c or d
  (and so on), and the leaves are worth to MAX 3 12 8, 2 4 6, 14 5 2."
  (let [leaves {[0 0] 3 [0 1] 12 [0 2] 8 [1 0] 2 [1 1] 4 [1 2] 6 [2 0] 14 [2 1] 5 [2 2] 2}]
    {:initial []
     :to-move (fn [path] (if (even? (count path)) :max :min))
     :actions (constantly [0 1 2])
     :result conj
     :terminal? #(= 2 (count %))
     :utility (fn [path player] (let [v (leaves path)] (if (= player :max) v (- v))))}))

(deftest figure-5-2-test
  (testing "minimax: MIN's values are 3, 2 and 2, so MAX plays a1 for 3"
    (let [r (adv/minimax figure-5-2 [])]
      (is (= 3 (:value r)))
      (is (= 0 (:best r)))
      (is (= {0 3 1 2 2 2} (:scores r)))))
  (testing "alpha-beta: the same move and value, and the two leaves the book shows it never looking at"
    (let [r (adv/alpha-beta figure-5-2 [])]
      (is (= 3 (:value r)))
      (is (= 0 (:best r)))
      (is (= 2 (:pruned r)))
      (is (= (- (:nodes (adv/minimax figure-5-2 [])) 2) (:nodes r)))))
  (testing "a step a visit: a leaf one, an inner node one and then one per child -- 9 + 4 x 4 = 25, after the first state"
    (is (= 13 (:nodes (:context (last (adv/minimax-seq figure-5-2 []))))))
    (is (= 26 (count (adv/minimax-seq figure-5-2 []))))))

(deftest expectiminimax-test
  ;; MAX chooses a1 or a2, then a coin weighted 0.9/0.1 decides between
  ;; two leaves. Doubling the gap at the far end of a2 changes nothing for
  ;; minimax, but everything for an expectation.
  (let [game (fn [leaves]
               {:initial []
                :to-move (fn [path] (case (count path) 0 :max 1 :chance :max))
                :actions (constantly [0 1])
                :chances (constantly [[0 0.9] [1 0.1]])
                :result conj
                :terminal? #(= 2 (count %))
                :utility (fn [path _] (leaves path))})
        small {[0 0] 2 [0 1] 3 [1 0] 1 [1 1] 4}
        large {[0 0] 20 [0 1] 30 [1 0] 1 [1 1] 400}]
    (testing "chance nodes are worth the expectation of their outcomes"
      (let [r (adv/minimax (game small) [])]
        (is (< (abs (- 2.1 (:value r))) 1e-12))
        (is (= 0 (:best r)))))
    (testing "and so an order-preserving change of the utilities can change the move"
      (is (= 1 (:best (adv/minimax (game large) [])))))))

;; ---------------------------------------------------------------------------
;; Tic-tac-toe

(deftest tic-tac-toe-test
  (let [ttt games/tic-tac-toe
        start (:initial ttt)]
    (testing "perfect play draws, by every search"
      (is (zero? (:value (adv/minimax ttt start))))
      (doseq [opts [{} {:order? true} {:order? true :table? true}]]
        (is (zero? (:value (adv/alpha-beta ttt start opts))))))
    (testing "the whole game tree is 549,946 nodes; alpha-beta needs a few percent of it"
      (let [full (:nodes (adv/minimax ttt start))
            pruned (:nodes (adv/alpha-beta ttt start))
            ordered (:nodes (adv/alpha-beta ttt start {:order? true :table? true}))]
        (is (= 549946 full))
        (is (< ordered pruned (* 0.05 full)))))
    (testing "a win in one is taken, and a loss in one blocked"
      (let [win (games/play ttt [0 3 1 4])         ; x to play 2 and win
            block (games/play ttt [0 3 8 4])]      ; x must block at 5
        (is (= 2 (:best (adv/alpha-beta ttt win))))
        (is (= 5 (:best (adv/alpha-beta ttt block))))
        (is (= 2 (:best (adv/monte-carlo ttt win 2000))))
        (is (= 5 (:best (adv/monte-carlo ttt block 2000))))))))

(deftest monte-carlo-against-random-test
  (testing "Monte Carlo, 400 playouts a move, never loses to a random player"
    (let [ttt games/tic-tac-toe]
      (doseq [k (range 1 11)]
        (let [rng (random/rng k)
              end (loop [s (:initial ttt)]
                    (if ((:terminal? ttt) s)
                      s
                      (recur ((:result ttt) s
                                            (if (= :x (:to-move s))
                                              (:best (adv/monte-carlo ttt s 400 {:rng rng}))
                                              (random/pick rng ((:actions ttt) s)))))))]
          (is (not= :o (:winner end)) k))))))

;; ---------------------------------------------------------------------------
;; Connect Four

(defn- random-positions [n moves]
  (for [k (range 1 (inc n))]
    (let [rng (random/rng k)]
      (loop [s (:initial games/connect-four) m 0]
        (if (or (= m moves) ((:terminal? games/connect-four) s))
          s
          (recur ((:result games/connect-four) s (random/pick rng ((:actions games/connect-four) s)))
                 (inc m)))))))

(deftest connect-four-test
  (let [c4 games/connect-four]
    (testing "cut off at the same depth, alpha-beta's value is minimax's, with or without ordering and a table"
      (doseq [s (random-positions 6 8)
              :when (not ((:terminal? c4) s))]
        (let [v (:value (adv/minimax c4 s {:depth 3}))]
          (doseq [opts [{:depth 3} {:depth 3 :order? true} {:depth 3 :order? true :table? true}]]
            (is (< (abs (- v (:value (adv/alpha-beta c4 s opts)))) 1e-12))))))
    (testing "and looks at fewer nodes"
      (let [s (:initial c4)]
        (is (< (:nodes (adv/alpha-beta c4 s {:depth 4})) (:nodes (adv/minimax c4 s {:depth 4}))))))
    (testing "the evaluation stays inside the utilities, so a real win outranks any guess"
      (doseq [s (random-positions 20 12) p [:x :o]]
        (is (< -1.0 ((:evaluate c4) s p) 1.0))))
    (testing "every search takes a win in one and blocks a loss in one"
      (let [win (games/play c4 [0 1 0 1 0 1])
            block (games/play c4 [0 1 0 1 0])]
        (doseq [choose [#(:best (adv/minimax c4 % {:depth 2}))
                        #(:best (adv/alpha-beta c4 % {:depth 4 :order? true :table? true}))
                        #(:best (:context (last (adv/deepening-seq c4 % {:max-depth 4}))))
                        #(:best (adv/monte-carlo c4 % 1500))]]
          (is (= 0 (choose win)))
          (is (= 0 (choose block))))))
    (testing "iterative deepening ends with the decision of the deepest search"
      (let [s (games/play c4 [3 3 2])
            end (last (adv/deepening-seq c4 s {:max-depth 5 :table? true}))
            direct (adv/alpha-beta c4 s {:depth 5 :order? true})]
        (is (= 5 (:depth end)))
        (is (< (abs (- (:value direct) (:value end))) 1e-12))))))

;; ---------------------------------------------------------------------------
;; Monte Carlo by line of play, and under hidden information (ISMCTS)

(deftest monte-carlo-by-path-test
  (testing "keyed by line of play rather than state, it still takes a win and blocks a loss"
    (let [ttt games/tic-tac-toe
          win (games/play ttt [0 3 1 4])
          block (games/play ttt [0 3 8 4])]
      (is (= 2 (:best (adv/monte-carlo ttt win 2000 {:node-key :path}))))
      (is (= 5 (:best (adv/monte-carlo ttt block 2000 {:node-key :path}))))))
  (testing "cut off and evaluated, on Connect Four"
    (let [c4 games/connect-four]
      (is (= 0 (:best (adv/monte-carlo c4 (games/play c4 [0 1 0 1 0 1]) 1500 {:node-key :path :cutoff 6}))))
      (is (= 0 (:best (adv/monte-carlo c4 (games/play c4 [0 1 0 1 0]) 1500 {:node-key :path :cutoff 6}))))))
  (testing "never loses to a random player at tic-tac-toe"
    (let [ttt games/tic-tac-toe]
      (doseq [k (range 1 11)]
        (let [rng (random/rng k)
              end (loop [s (:initial ttt)]
                    (if ((:terminal? ttt) s)
                      s
                      (recur ((:result ttt) s
                                            (if (= :x (:to-move s))
                                              (:best (adv/monte-carlo ttt s 400 {:rng rng :node-key :path}))
                                              (random/pick rng ((:actions ttt) s)))))))]
          (is (not= :o (:winner end)) k)))))
  (testing "the same seed, the same search"
    (let [s (games/play games/connect-four [3 3])
          run #(adv/monte-carlo games/connect-four s 300 {:rng (random/rng 5) :cutoff 4 :node-key :path})]
      (is (= (run) (run))))))

(def ^:private coin-game
  "A coin lies hidden under O's hand; X may guess it, for +1 if right and
  -1 if wrong, or pass, for a sure +0.2. Guessing is worth nothing to a
  player who cannot see the coin, and everything to one who can."
  (let [value {:pass 0.2}]
    {:to-move :to-move
     :actions (fn [{:keys [done]}] (if done [] [:heads :tails :pass]))
     :result (fn [s a] (assoc s :done a :to-move :o))
     :terminal? (comp some? :done)
     :utility (fn [{:keys [coin done]} player]
                (let [u (get value done (if (= done coin) 1.0 -1.0))]
                  (if (= player :x) u (- u))))
     :determinize (fn [s _player rng] (assoc s :coin (random/pick rng [:heads :tails])))}))

(deftest monte-carlo-does-not-peek-test
  (let [s {:coin :tails :to-move :x :done nil}]
    (testing "a game with :determinize is searched by line of play unless told otherwise"
      (is (= :path (:node-key (first (adv/monte-carlo-seq coin-game s)))))
      (is (= :state (:node-key (first (adv/monte-carlo-seq (dissoc coin-game :determinize) s))))))
    (testing "searching what X knows, X passes"
      (is (= :pass (:best (adv/monte-carlo coin-game s 2000)))))
    (testing "searching the true state -- cheating -- X calls the coin"
      (is (= :tails (:best (adv/monte-carlo (dissoc coin-game :determinize) s 2000)))))
    (testing "and the values say why: a guess is worth about nothing blind"
      (let [{:keys [moves]} (adv/monte-carlo coin-game s 4000)]
        (is (< (abs (:q (moves :heads))) 0.2))
        (is (< (abs (:q (moves :tails))) 0.2))))))

(def ^:private sometimes-game
  "X may always :wait (worth 0), and may :strike (worth +0.5) only when
  the hidden die is even -- which X cannot see."
  {:to-move :to-move
   :actions (fn [{:keys [done die]}] (cond done [] (even? die) [:wait :strike] :else [:wait]))
   :result (fn [s a] (assoc s :done a))
   :terminal? (comp some? :done)
   :utility (fn [{:keys [done]} p] (let [u (if (= done :strike) 0.5 0.0)] (if (= p :x) u (- u))))
   :determinize (fn [s _ rng] (assoc s :die (inc (random/below rng 6))))})

(deftest monte-carlo-availability-test
  (let [st (nth (adv/monte-carlo-seq sometimes-game {:die 2 :to-move :x :done nil} {}) 2000)
        {:keys [moves best]} (adv/monte-carlo-choice st)]
    (testing "an action legal in about half the worlds was available about half the time"
      (is (< 800 (:avail (moves :strike)) 1200))
      (is (= 2000 (:avail (moves :wait)))))
    (testing "every count is no more than its availability"
      (doseq [[_ {:keys [n avail]}] moves] (is (<= n avail))))
    (testing "and where it is legal it is the better move, so it is the most played of its chances"
      (is (> (:q (moves :strike)) (:q (moves :wait))))
      (is (= :wait best) "wait is always there, so it is played in every world where strike is not"))))
