(ns allgo.decision.contract-net-test
  (:require [allgo.decision.contract-net :as cn]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; The contract net, after Russell and Norvig, *Artificial Intelligence: A
;; Modern Approach*, 4th edition, chapter 18; Smith's protocol with
;; Sandholm's marginal-cost bids. Checked against invariants of the
;; protocol and, on small fleets, against the best allocation found by
;; trying every one.

(def ^:private fleet
  [{:id :a :home [0 0] :skills #{:truck}}
   {:id :b :home [10 0] :skills #{:truck}}
   {:id :c :home [5 8] :skills #{:truck :crane}}])

(def ^:private jobs
  [{:id 1 :at [1 1] :needs :truck} {:id 2 :at [9 1] :needs :truck} {:id 3 :at [5 7] :needs :crane}
   {:id 4 :at [2 2] :needs :truck} {:id 5 :at [8 3] :needs :truck}
   {:id 6 :at [6 6] :needs :truck :subtasks [{:id 7 :at [4 6] :needs :crane}]}])

(defn- perms [xs] (if (empty? xs) [[]] (for [x xs p (perms (remove #{x} xs))] (cons x p))))

(defn- best-cost
  "The shortest route from `home` through `pts`, by trying every order."
  [home pts]
  (if (empty? pts) 0.0 (apply min (map #(cn/route-cost home %) (perms pts)))))

(defn- subsets [xs] (reduce (fn [ss x] (concat ss (map #(conj % x) ss))) [#{}] xs))

(defn- optimum
  "The least total route length over every assignment of `tasks` to
  `agents` and every order of each route."
  [agents tasks]
  (apply min (for [assign (reduce (fn [acc _] (for [a acc g (range (count agents))] (conj a g))) [[]] tasks)]
               (reduce + (map-indexed (fn [g ag] (best-cost (:home ag) (map :at (keep-indexed (fn [i t] (when (= g (assign i)) t)) tasks))))
                                      agents)))))

(defn- random-jobs [seed n]
  (let [rng (random/rng seed)] (vec (for [i (range n)] {:id i :at [(* 10 (rng)) (* 8 (rng))]}))))

(deftest protocol
  (let [states (vec (cn/contract-net-seq fleet jobs :dispatcher))
        end (peek states)
        msgs (keep :delivered states)]
    (testing "every task placed with an agent able to do it"
      (is (empty? (:unplaced end)))
      (is (= #{1 2 3 4 5 6 7} (set (keys (:awarded end)))))
      (is (= :c ((:awarded end) 3)) "only c has a crane")
      (is (= :c ((:awarded end) 7))))
    (testing "a message a step, every one handled"
      (is (empty? (:queue end)))
      (is (= (:messages end) (count msgs))))
    (testing "per task: an announcement and a reply from each agent, one award, a rejection for each other bidder"
      (doseq [t (keys (:awarded end))]
        (let [of-t (group-by :type (filter #(= t (:task %)) msgs))]
          (is (= 3 (count (:announce of-t))))
          (is (= 3 (+ (count (:bid of-t)) (count (:decline of-t)))))
          (is (= 1 (count (:award of-t))))
          (is (= (dec (count (:bid of-t))) (count (:reject of-t)))))))
    (testing "the cranes decline what they cannot do: a and b decline the crane jobs"
      (is (= #{:a :b} (set (map :from (filter #(and (= :decline (:type %)) (= 3 (:task %))) msgs))))))
    (testing "the subtask is announced by the contractor of its parent, not the dispatcher"
      (is (= #{((:awarded end) 6)} (set (map :from (filter #(and (= :announce (:type %)) (= 7 (:task %))) msgs))))))
    (testing "the award goes to the cheapest bid, and the winning bids add up to the total cost"
      (let [won (for [{:keys [task to]} (filter #(= :award (:type %)) msgs)]
                  (:cost (first (filter #(and (= :bid (:type %)) (= task (:task %)) (= to (:from %))) msgs))))]
        (doseq [t (keys (:awarded end))]
          (let [bids (filter #(and (= :bid (:type %)) (= t (:task %))) msgs)]
            (is (= ((:awarded end) t) (:from (apply min-key :cost bids))))))
        (is (< (abs (- (reduce + won) (cn/total-cost end))) 1e-9))))))

(deftest capacity
  (testing "an agent at capacity declines; a task nobody can take is left unplaced"
    (let [end (cn/allocate [{:id :a :home [0 0] :capacity 1}] [{:id 1 :at [1 0]} {:id 2 :at [2 0]}] :m)]
      (is (= {:a [1]} (:routes end)))
      (is (= [2] (:unplaced end))))))

(deftest best-route
  (testing "Held and Karp's route is the shortest of every order"
    (let [rng (random/rng 3)]
      (doseq [n (range 8) _ (range 4)]
        (let [pts (vec (repeatedly n #(vector (* 10 (rng)) (* 8 (rng)))))
              home [(rng) (rng)]
              [order c] (cn/best-route home pts)]
          (is (= (range n) (sort order)))
          (is (< (abs (- c (cn/route-cost home (map pts order)))) 1e-9))
          (is (< (abs (- c (best-cost home pts))) 1e-9)))))))

(defn- set-cost
  "What agent `a` pays for tasks `ts` in state `st`: its shortest route."
  [st a ts]
  (best-cost (:home ((:agents st) a)) (map #(:at ((:tasks st) %)) ts)))

(defn- final [st & opts] (peek (vec (apply cn/recontract-seq st opts))))

(deftest against-the-optimum
  (let [agents (mapv #(dissoc % :skills) fleet)]
    (testing "never better than the best allocation; recontracting, of any type, never worse"
      (doseq [seed (range 6)
              contracts [#{:o} #{:s} #{:c} #{:o :s :c}]]
        (let [tasks (random-jobs seed 6)
              end (cn/allocate agents tasks :m)
              re (vec (cn/recontract-seq end :contracts contracts))]
          (is (<= (- (optimum agents tasks) 1e-9) (cn/total-cost (peek re)) (cn/total-cost end)) [seed contracts])
          (is (apply >= (map cn/total-cost re)) [seed contracts]))))
    (testing "where one-task contracts stop short, swaps and clusters reach the optimum"
      (doseq [seed [0 1 4 6]]
        (let [end (cn/allocate agents (random-jobs seed 6) :m)
              best (optimum agents (random-jobs seed 6))]
          (is (< (+ best 0.5) (cn/total-cost (final end :contracts #{:o}))) seed)
          (is (< (abs (- best (cn/total-cost (final end)))) 1e-9) seed))))
    (testing "at the end no task, cluster of tasks or swap lowers the cost"
      (doseq [seed (range 6)
              :let [st (final (cn/allocate agents (random-jobs seed 6) :m))
                    {:keys [routes]} st
                    held (fn [a] (set (routes a)))]
              from (keys routes)
              to (keys routes) :when (not= from to)
              :let [before (+ (set-cost st from (held from)) (set-cost st to (held to)))]]
        (doseq [moving (rest (subsets (held from)))]
          (is (<= before (+ 1e-9 (set-cost st from (remove moving (held from))) (set-cost st to (into (held to) moving))))))
        (doseq [t (held from) u (held to)]
          (is (<= before (+ 1e-9 (set-cost st from (conj (disj (held from) t) u)) (set-cost st to (conj (disj (held to) u) t))))))))
    (testing "and a local optimum of all three that only a contract among three agents leaves"
      (let [tasks (random-jobs 10 6)
            end (cn/allocate agents tasks :m)
            re (vec (cn/recontract-seq end))]
        (is (not-any? :moved re))
        (is (= {:a [2] :b [3 1 4 5] :c [0]} (:routes (peek re))))
        ;; b handing 1, 4 and 5 to c as c hands 0 to a
        (let [st (assoc end :routes {:a [0 2] :b [3] :c [1 4 5]})]
          (is (< (+ (optimum agents tasks) -1e-9) (reduce + (map #(set-cost st % ((:routes st) %)) [:a :b :c]))
                 (+ (optimum agents tasks) 1e-9)))
          (is (< (+ (optimum agents tasks) 0.4) (cn/total-cost (peek re)))))))))

(deftest contract-types
  (let [agents [{:id :a :home [0 0] :capacity 1} {:id :b :home [10 0] :capacity 1}]
        tasks [{:id 1 :at [9 1]} {:id 2 :at [1 1]}]
        ;; each holding the task nearer the other
        st (assoc (cn/allocate agents tasks :m) :routes {:a [1] :b [2]})]
    (testing "a swap: full, neither can take the other's task, but they can trade"
      (is (not-any? :moved (cn/recontract-seq st :contracts #{:o :c})))
      (let [re (vec (cn/recontract-seq st :contracts #{:s}))
            moved (some :moved re)]
        (is (= [:s :a #{1} :b #{2}] ((juxt :type :from :tasks :to :back) moved)))
        (is (= {:a [2] :b [1]} (:routes (peek re)))))))
  (let [agents [{:id :a :home [0 0]} {:id :b :home [10 0]}]
        tasks [{:id 1 :at [5 0]} {:id 2 :at [6 0]}]
        st (assoc (cn/allocate agents tasks :m) :routes {:a [1 2] :b []})]
    (testing "a cluster: neither task alone is worth moving, both together are"
      (is (not-any? :moved (cn/recontract-seq st :contracts #{:o :s})))
      (let [re (vec (cn/recontract-seq st :contracts #{:c}))
            moved (some :moved re)]
        (is (= [:c :a #{1 2} :b] ((juxt :type :from :tasks :to) moved)))
        (is (< (abs (- 1.0 (:gain moved))) 1e-9) "a saves 6, b's route grows by 5")
        (is (= {:a [] :b [2 1]} (:routes (peek re))))))))
