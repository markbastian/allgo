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

(defn- optimum
  "The least total route length over every assignment of `tasks` to
  `agents` and every order of each route."
  [agents tasks]
  (let [best-route (fn [home pts] (if (empty? pts) 0.0 (apply min (map #(cn/route-cost home %) (perms pts)))))]
    (apply min (for [assign (reduce (fn [acc _] (for [a acc g (range (count agents))] (conj a g))) [[]] tasks)]
                 (reduce + (map-indexed (fn [g ag] (best-route (:home ag) (map :at (keep-indexed (fn [i t] (when (= g (assign i)) t)) tasks))))
                                        agents))))))

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

(deftest against-the-optimum
  (let [agents (mapv #(dissoc % :skills) fleet)]
    (testing "never better than the best allocation, often close; recontracting never worse"
      (doseq [seed (range 6)]
        (let [tasks (random-jobs seed 6)
              end (cn/allocate agents tasks :m)
              re (vec (cn/recontract-seq end))
              best (optimum agents tasks)]
          (is (<= (- best 1e-9) (cn/total-cost (peek re)) (cn/total-cost end)) seed)
          (is (apply >= (map cn/total-cost re)) seed))))
    (testing "recontracting moves a task when another agent can do it for less than its holder saves"
      (let [tasks (random-jobs 9 6)
            end (cn/allocate agents tasks :m)
            re (vec (cn/recontract-seq end))]
        (is (seq (keep :moved re)))
        (is (< (cn/total-cost (peek re)) (cn/total-cost end)))
        (is (< (- (cn/total-cost (peek re)) (optimum agents tasks)) 0.05))))
    (testing "at the end no single task can be moved at a profit"
      (let [final (peek (vec (cn/recontract-seq (cn/allocate agents (random-jobs 0 6) :m))))]
        (doseq [[holder route] (:routes final) t route
                :let [without (assoc-in final [:routes holder] (vec (remove #{t} route)))
                      saving (- (cn/agent-cost final holder) (cn/agent-cost without holder))]
                other (keys (:agents final)) :when (not= other holder)]
          (is (>= (+ 1e-9 (second (cn/cheapest-insertion (:home ((:agents final) other))
                                                         (map #(:at ((:tasks final) %)) ((:routes without) other))
                                                         (:at ((:tasks final) t)))))
                  saving)))))))
