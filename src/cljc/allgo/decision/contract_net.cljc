(ns allgo.decision.contract-net
  "The contract net: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 18 -- allocating tasks among
  cooperative agents by announcement, bidding and award. The protocol is
  Smith's (\"The Contract Net Protocol\", IEEE Transactions on Computers,
  1980); bidding marginal costs, and recontracting to improve on the
  first allocation, are Sandholm's (\"An Implementation of the Contract
  Net Protocol Based on Marginal Cost Calculations\", AAAI 1993).

    announce   a manager broadcasts a task, with what a contractor must be
               able to do
    bid        each agent able to take it replies with what taking it
               would cost it -- here the length it would add to its
               route, at the cheapest place in it -- or declines: not
               able, or full
    award      when every reply is in, the manager awards the task to the
               cheapest bid and tells the rest they lost
    subtasks   a contractor awarded a task with parts becomes their
               manager, and announces them in turn: the net is recursive
    recontract once everything is placed, each task is offered again by
               its holder, and moved to another agent if that agent's bid
               is less than what dropping it saves the holder -- a
               hill-climb on total cost, since every move lowers it

  Every message is a step: the state's `:queue` holds those in flight,
  `:delivered` the one just handled. An agent is `{:id :home [x y]
  :skills #{...} :capacity n}`, a task `{:id :at [x y] :needs skill
  :subtasks [tasks]}`; an agent's route starts at its home and visits its
  tasks in order, and its cost is the route's length."
  (:require [allgo.geometry.polygon :as polygon])
  #?(:clj (:import (clojure.lang PersistentQueue))))

(def ^:private empty-queue #?(:clj PersistentQueue/EMPTY :cljs #queue []))

;; ---------------------------------------------------------------------------
;; Routes and costs

(defn route-cost
  "The length of a route from `home` through the points `stops`."
  ^double [home stops]
  (reduce + 0.0 (map (fn [[a b]] (polygon/distance a b)) (partition 2 1 (cons home stops)))))

(defn cheapest-insertion
  "Where in `stops` -- after `home` -- point `p` adds least to the route,
  and how much: `[index added]`."
  [home stops p]
  (let [stops (vec stops)
        base (route-cost home stops)]
    (apply min-key second
           (for [i (range (inc (count stops)))]
             [i (- (route-cost home (concat (subvec stops 0 i) [p] (subvec stops i))) base)]))))

(defn- points [tasks route] (map #(:at (tasks %)) route))

(defn agent-cost
  "What agent `id`'s route costs it in state `st`."
  ^double [{:keys [agents tasks routes]} id]
  (route-cost (:home (agents id)) (points tasks (routes id))))

(defn total-cost
  "The sum of every agent's route."
  ^double [st]
  (reduce + 0.0 (map #(agent-cost st %) (keys (:agents st)))))

(defn- eligible? [{:keys [agents routes]} id task]
  (let [{:keys [skills capacity]} (agents id)]
    (and (or (nil? (:needs task)) (contains? skills (:needs task)))
         (or (nil? capacity) (< (count (routes id)) (long capacity))))))

(defn- bid
  "Agent `id`'s reply to an announcement of `task`: its marginal cost, or
  nil to decline."
  [{:keys [agents tasks routes] :as st} id task]
  (when (eligible? st id task)
    (second (cheapest-insertion (:home (agents id)) (points tasks (routes id)) (:at task)))))

(defn- take-task [{:keys [agents tasks routes] :as st} id task-id]
  (let [stops (vec (routes id))
        [i] (cheapest-insertion (:home (agents id)) (points tasks stops) (:at (tasks task-id)))]
    (assoc-in st [:routes id] (into (conj (subvec stops 0 i) task-id) (subvec stops i)))))

;; ---------------------------------------------------------------------------
;; The protocol

(defn- announce
  "Manager `from` broadcasts `task`: an announcement queued to every
  other agent (and the manager too, if it is an agent: it may bid on its
  own subtasks)."
  [st from task]
  (let [to (keys (:agents st))]
    (-> st
        (assoc-in [:open (:id task)] {:manager from :awaiting (set to) :bids {}})
        (update :queue into (map (fn [a] {:type :announce :from from :to a :task (:id task)}) to)))))

(defn- next-announcement [{:keys [pending] :as st}]
  (if-let [[from task] (first pending)]
    (announce (assoc st :pending (subvec pending 1)) from task)
    st))

(defn- handle
  [{:keys [tasks] :as st} {:keys [type from to task] :as msg}]
  (case type
    :announce (let [c (bid st to (tasks task))]
                (update st :queue conj (if c
                                         {:type :bid :from to :to from :task task :cost c}
                                         {:type :decline :from to :to from :task task})))
    (:bid :decline)
    (let [st (-> st
                 (update-in [:open task :awaiting] disj from)
                 (cond-> (= type :bid) (assoc-in [:open task :bids from] (:cost msg))))
          {:keys [awaiting bids manager]} (get-in st [:open task])]
      (if (seq awaiting)
        st
        ;; Every reply is in: award to the cheapest, or leave it unplaced.
        (let [st (update st :open dissoc task)]
          (if (empty? bids)
            (-> st (update :unplaced conj task) next-announcement)
            (let [winner (key (apply min-key val (sort-by key bids)))]
              (update st :queue into (cons {:type :award :from manager :to winner :task task}
                                           (for [a (keys bids) :when (not= a winner)]
                                             {:type :reject :from manager :to a :task task}))))))))
    :award (let [st (-> st (take-task to task) (assoc-in [:awarded task] to))
                 subs (:subtasks (tasks task))]
             ;; A task with parts makes its contractor their manager.
             (next-announcement (if (seq subs)
                                  (update st :pending #(into (mapv (fn [s] [to s]) subs) %))
                                  st)))
    :reject st))

(defn- all-tasks [tasks] (mapcat (fn [t] (cons t (all-tasks (:subtasks t)))) tasks))

(defn contract-net-seq
  "The contract net allocating `tasks`, announced in order by the manager
  `manager` (not itself an agent unless it is one of `agents`), a message
  a step. Each state: `:queue` the messages in flight, `:delivered` the
  one just handled, `:routes` each agent's tasks in order, `:awarded`
  task to agent, `:unplaced` the tasks no one could take, `:open` the
  announcements still collecting replies. Ends when no message is left."
  [agents tasks manager]
  (let [agents (into {} (map (juxt :id identity)) agents)
        st (-> {:agents agents
                :tasks (into {} (map (juxt :id identity)) (all-tasks tasks))
                :routes (zipmap (keys agents) (repeat []))
                :queue empty-queue :pending (mapv (fn [t] [manager t]) tasks)
                :open {} :awarded {} :unplaced [] :messages 0}
               next-announcement)]
    (cons st
          (->> st
               (iterate (fn [{:keys [queue] :as st}]
                          (when-let [msg (peek queue)]
                            (-> st (assoc :queue (pop queue) :delivered msg)
                                (update :messages inc)
                                (handle msg)))))
               rest
               (take-while some?)))))

(defn allocate
  "The final state of the contract net on `agents` and `tasks`."
  [agents tasks manager]
  (last (contract-net-seq agents tasks manager)))

;; ---------------------------------------------------------------------------
;; Routes planned whole

(def ^:private exact-limit
  "The most stops `best-route` plans exactly: Held and Karp's program
  takes 2^n n^2 steps."
  12)

(defn- held-karp
  "The shortest path from `home` through every point of `pts`, by Held and
  Karp's dynamic program: for each set of points and each one of them,
  the shortest path from home through the set that ends there, built up
  from the sets one smaller. `[order cost]`."
  [home pts]
  (let [pts (vec pts)
        n (count pts)
        full (dec (bit-shift-left 1 n))
        ;; table: mask -> point j -> [cost previous-point], j in the mask
        table (reduce (fn [table mask]
                        (conj table
                              (vec (for [j (range n)]
                                     (cond (not (bit-test mask j)) nil
                                           (= mask (bit-set 0 j)) [(polygon/distance home (pts j)) nil]
                                           :else (let [m (bit-clear mask j)]
                                                   (apply min-key first
                                                          (for [k (range n) :when (bit-test m k)]
                                                            [(+ (double (first ((table m) k))) (polygon/distance (pts k) (pts j))) k]))))))))
                      [[]] (range 1 (inc full)))
        last-j (apply min-key #(first ((table full) %)) (range n))
        order (loop [mask full j last-j order ()]
                (if (nil? j)
                  order
                  (recur (bit-clear mask j) (second ((table mask) j)) (cons j order))))]
    [(vec order) (first ((table full) last-j))]))

(defn best-route
  "The shortest route from `home` through every one of `stops` (points),
  as `[order cost]`, `order` indices into `stops`. Exact for up to a dozen
  stops; past that, the stops inserted one at a time each at its
  cheapest place."
  [home stops]
  (let [stops (vec stops)
        n (count stops)]
    (cond (zero? n) [[] 0.0]
          (<= n exact-limit) (held-karp home stops)
          :else (let [order (reduce (fn [order i]
                                      (let [[at] (cheapest-insertion home (map stops order) (stops i))]
                                        (into (conj (subvec order 0 at) i) (subvec order at))))
                                    [] (range n))]
                  [order (route-cost home (map stops order))]))))

;; ---------------------------------------------------------------------------
;; Recontracting

(defn- subsets
  "Every subset of `xs` with between `lo` and `hi` members, smallest first."
  [xs lo hi]
  (let [xs (vec xs)]
    (for [k (range lo (inc (min hi (count xs))))
          s (letfn [(choose [i k] (cond (zero? k) [[]]
                                        (> k (- (count xs) i)) []
                                        :else (concat (map #(cons (xs i) %) (choose (inc i) (dec k)))
                                                      (choose (inc i) k))))]
              (choose 0 k))]
      (set s))))

(defn- offers
  "A pass's offers, in the order they are made: every one-task offer,
  then every task offered for a swap, then every cluster of two or more
  of an agent's tasks, up to `max-cluster`."
  [{:keys [routes]} contracts max-cluster]
  (let [held (sort-by key routes)]
    (concat (when (contracts :o) (for [[a r] held t r] {:type :o :from a :tasks #{t}}))
            (when (contracts :s) (for [[a r] held t r] {:type :s :from a :tasks #{t}}))
            (when (contracts :c) (for [[a r] held T (subsets r 2 (or max-cluster (count r)))]
                                   {:type :c :from a :tasks T})))))

(defn- replies
  "What each other agent would give to take `offer`'s tasks: `[gain to
  back]`, the fall in the two agents' summed costs and, for a swap, the
  task it hands back. `cost` values an agent's set of tasks."
  [{:keys [agents tasks routes]} cost {type :type from :from moving :tasks}]
  (let [able? (fn [a ts n] (let [{:keys [skills capacity]} (agents a)]
                             (and (every? #(or (nil? (:needs (tasks %))) (contains? skills (:needs (tasks %)))) ts)
                                  (or (nil? capacity) (<= n (long capacity))))))
        mine (set (routes from))
        before-from (cost from mine)
        after-from (cost from (reduce disj mine moving))]
    (for [b (sort (keys agents)) :when (not= b from)
          :let [theirs (set (routes b))]
          back (if (= type :s) (map vector (routes b)) [[]])
          :let [back (set back)]
          :when (and (able? b moving (+ (count theirs) (count moving) (- (count back))))
                     (able? from back (+ (count mine) (count back) (- (count moving)))))]
      [(- (+ before-from (cost b theirs))
          (if (= type :s) (cost from (into (disj mine (first moving)) back)) after-from)
          (cost b (into (reduce disj theirs back) moving)))
       b back])))

(defn- replan
  "Agent `a`'s route made the best through `ts`."
  [{:keys [agents tasks] :as st} a ts]
  (let [ts (vec (sort ts))
        [order] (best-route (:home (agents a)) (map #(:at (tasks %)) ts))]
    (assoc-in st [:routes a] (mapv ts order))))

(defn- recontract-step
  [{:keys [pending routes] :as st} cost]
  (when-let [{:keys [from tasks] :as offer} (first pending)]
    (let [st (assoc st :pending (rest pending) :offer offer :moved nil)]
      (if-not (every? (set (routes from)) tasks)
        st
        (let [best (when-let [rs (seq (replies st cost offer))]
                     (apply max-key first rs))
              [gain to back] best]
          (if (and best (> (double gain) 1e-9))
            (-> st
                (replan from (into (reduce disj (set (routes from)) tasks) back))
                (replan to (into (reduce disj (set (routes to)) back) tasks))
                (assoc :moved (assoc offer :to to :back back :gain gain)))
            st))))))

(defn recontract-seq
  "Recontracting from the allocation `st`, an offer a step, after
  Sandholm's contract types (\"Contract Types for Satisficing Task
  Allocation: I Theoretical Results\", AAAI Spring Symposium, 1998):

    :o  one task, offered by its holder to the others -- the contract
        net's own contract
    :s  a swap: a task offered in exchange for one of the taker's
    :c  a cluster: two or more of the holder's tasks, offered as one

  Each agent first replans its route, the best through its tasks, so
  that what an agent's work costs it is a matter of which tasks it holds
  and not the order they came in. An offer goes to every other agent;
  the one whose reply lowers the two agents' summed cost most takes it,
  if anything lowers it -- the contract individually rational to both
  once the gain is split by a payment -- and both replan. Passes of
  every offer repeat until one moves nothing. `:contracts` (default
  all three) chooses the types, `:max-cluster` caps a cluster's size.
  Each state's `:offer` is `{:type :from :tasks}`, `:moved` that with
  `:to`, `:back` (a swap's return) and `:gain` when it was taken. Total
  cost never rises.

  Sandholm shows that each type alone can stop short of the best
  allocation, and that only a contract free to move any tasks among any
  agents at once -- his OCSM-contract -- always reaches it. These three
  together can stop short too: when it takes three agents trading at
  once."
  [st & {:keys [contracts max-cluster] :or {contracts #{:o :s :c}}}]
  (let [costs (atom {})
        cost (fn [a ts]
               (let [k [a ts]]
                 (or (@costs k)
                     (let [{:keys [agents tasks]} st
                           c (second (best-route (:home (agents a)) (map #(:at (tasks %)) (sort ts))))]
                       (swap! costs assoc k c)
                       c))))
        planned (assoc (reduce (fn [st a] (replan st a (set (get-in st [:routes a])))) st (keys (:agents st)))
                       :replanned true :offer nil :moved nil)]
    (letfn [(pass [st]
              (let [states (->> (assoc st :pending (offers st contracts max-cluster) :replanned false)
                                (iterate #(some-> % (recontract-step cost)))
                                rest
                                (take-while some?))]
                (when (seq states)
                  (lazy-cat states (when (some :moved states) (pass (last states)))))))]
      (cons st (cons planned (pass planned))))))
