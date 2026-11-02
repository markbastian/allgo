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
  the shortest path from home through the set that ends there, each
  extended by one point into the sets one larger. `[order cost]`."
  [home pts]
  (let [pts (vec pts)
        n (count pts)
        size (bit-shift-left 1 n)
        full (dec size)
        dist (double-array (* n n))
        ;; cost and previous point, at (mask * n + last point)
        cost (double-array (* size n) ##Inf)
        prev (long-array (* size n) -1)]
    (dotimes [i n]
      (dotimes [j n] (aset dist (+ (* i n) j) (polygon/distance (pts i) (pts j))))
      (aset cost (+ (* (bit-shift-left 1 i) n) i) (polygon/distance home (pts i))))
    (doseq [mask (range 1 size) j (range n)
            :let [c (aget cost (+ (* (long mask) n) (long j)))]
            :when (< c ##Inf)
            k (range n) :when (not (bit-test mask k))
            :let [at (+ (* (bit-set mask k) n) (long k))
                  v (+ c (aget dist (+ (* (long j) n) (long k))))]
            :when (< v (aget cost at))]
      (aset cost at v)
      (aset prev at (long j)))
    (let [last-j (apply min-key #(aget cost (+ (* full n) (long %))) (range n))]
      [(loop [mask (long full) j (long last-j) order ()]
         (if (neg? j)
           (vec order)
           (recur (bit-clear mask j) (aget prev (+ (* mask n) j)) (cons j order))))
       (aget cost (+ (* full n) (long last-j)))])))

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
  every task offered for a swap, every cluster of two or more of an
  agent's tasks up to `max-cluster`, every task offered to start an
  M-contract, and every set of tasks, smallest first, up to `max-ocsm`,
  put up for an OCSM-contract."
  [{:keys [routes]} contracts max-cluster max-ocsm]
  (let [held (sort-by key routes)
        all (sort (mapcat val held))]
    (concat (when (contracts :o) (for [[a r] held t r] {:type :o :from a :tasks #{t}}))
            (when (contracts :s) (for [[a r] held t r] {:type :s :from a :tasks #{t}}))
            (when (contracts :c) (for [[a r] held T (subsets r 2 (or max-cluster (count r)))]
                                   {:type :c :from a :tasks T}))
            (when (contracts :m) (for [[a r] held t r] {:type :m :from a :tasks #{t}}))
            (when (contracts :ocsm) (for [T (subsets all 1 (or max-ocsm (count all)))]
                                      {:type :ocsm :tasks T})))))

(defn- candidates
  "The contracts `offer` admits, each a map of task to the agent it
  moves to. `holder` maps each task to the agent holding it."
  [{:keys [agents routes]} holder {:keys [type from tasks]}]
  (let [others (fn [a] (remove #{a} (sort (keys agents))))
        t (first tasks)]
    (case type
      :o (for [b (others from)] {t b})
      :c (for [b (others from)] (zipmap tasks (repeat b)))
      :s (for [b (others from) u (routes b)] {t b u from})
      ;; Three or more agents, each handing on at most one task: `from`
      ;; hands on `t`, each other agent one of its own or none.
      :m (for [b (others from)
               moves (reduce (fn [ms x]
                               (for [m ms
                                     choice (cons nil (for [u (routes x) y (others x)] [u y]))]
                                 (if choice (assoc m (first choice) (second choice)) m)))
                             [{t b}] (others from))
               :when (<= 3 (count (into (set (vals moves)) (map holder (keys moves)))))]
           moves)
      ;; Every reallocation of the tasks among every agent.
      :ocsm (->> (reduce (fn [ms u] (for [m ms y (sort (keys agents))] (assoc m u y))) [{}] (sort tasks))
                 (remove (fn [m] (every? (fn [[u y]] (= y (holder u))) m)))))))

(defn- outcome
  "What making the contract `moves` does: the sets of tasks the agents
  party to it end with, and the fall in their summed cost; nil if one of
  them cannot do, or hold, what it would end with. `owned` maps each
  agent to the set of tasks it holds."
  [{:keys [agents tasks]} holder owned cost moves]
  ;; Whether each can do and hold what it is given, from the moves alone,
  ;; before any route is planned.
  (when (and (every? (fn [[u y]] (let [needs (:needs (tasks u))]
                                   (or (nil? needs) (contains? (:skills (agents y)) needs))))
                     moves)
             (every? (fn [[a n]] (let [capacity (:capacity (agents a))]
                                   (or (nil? capacity) (<= (+ (count (owned a)) (long n)) (long capacity)))))
                     (reduce (fn [d [u y]] (-> d (update (holder u) (fnil dec 0)) (update y (fnil inc 0)))) {} moves)))
    (let [sets (reduce (fn [sets [u y]]
                         (let [x (holder u)]
                           (-> sets
                               (update x (fnil disj (owned x)) u)
                               (update y (fnil conj (owned y)) u))))
                       {} moves)]
      {:sets sets
       :gain (- (reduce + (map #(cost % (owned %)) (keys sets)))
                (reduce + (map (fn [[a ts]] (cost a ts)) sets)))})))

(defn- replan
  "Agent `a`'s route made the best through `ts`."
  [{:keys [agents tasks] :as st} a ts]
  (let [ts (vec (sort ts))
        [order] (best-route (:home (agents a)) (map #(:at (tasks %)) ts))]
    (assoc-in st [:routes a] (mapv ts order))))

(defn- weigh
  "Weighs up to `chunk` more of the candidates of the offer in hand; when
  they run out, makes the best contract, if any lowers the cost."
  [{:keys [weighing] :as st} cost chunk]
  (let [{:keys [offer holder owned left best]} weighing
        now (take chunk left)
        later (drop chunk left)
        best (reduce (fn [best moves]
                       (let [o (outcome st holder owned cost moves)]
                         (if (and o (> (double (:gain o)) (double (:gain best 1e-9))))
                           (assoc o :moves moves)
                           best)))
                     best now)
        st (update st :weighed + (count now))]
    (cond (seq later) (assoc st :weighing (assoc weighing :left later :best best))
          best (-> (reduce (fn [st [a ts]] (replan st a ts)) st (:sets best))
                   (assoc :weighing nil :concluded true
                          :moved (assoc offer :gain (:gain best)
                                        :moves (into (sorted-map) (for [[u y] (:moves best) :when (not= y (holder u))]
                                                                    [u [(holder u) y]])))))
          :else (assoc st :weighing nil :concluded true))))

(defn- recontract-step
  [{:keys [pending routes weighing] :as st} cost chunk]
  (let [st (assoc st :moved nil :concluded false)]
    (if weighing
      (weigh st cost chunk)
      (when-let [{:keys [type from tasks] :as offer} (first pending)]
        (let [st (assoc st :pending (rest pending) :offer offer :weighed 0)
              holder (into {} (for [[a r] routes t r] [t a]))
              owned (update-vals routes set)]
          (if-not (if (= type :ocsm) (every? holder tasks) (every? (owned from) tasks))
            (assoc st :concluded true)
            (weigh (assoc st :weighing {:offer offer :holder holder :owned owned
                                        :left (candidates st holder offer) :best nil})
                   cost chunk)))))))

(defn recontract-seq
  "Recontracting from the allocation `st`, after
  Sandholm's contract types (\"Contract Types for Satisficing Task
  Allocation: I Theoretical Results\", AAAI Spring Symposium, 1998):

    :o     one task, offered by its holder to the others -- the contract
           net's own contract
    :s     a swap: a task offered in exchange for one of the taker's
    :c     a cluster: two or more of the holder's tasks, offered as one
    :m     multiagent: three or more agents at once, each handing on at
           most one task -- offered by the holder of one of them
    :ocsm  any reallocation of a set of tasks among any agents: all of
           the above at once

  Each agent first replans its route, the best through its tasks, so
  that what an agent's work costs it is a matter of which tasks it holds
  and not the order they came in. Of the contracts an offer admits, the
  one that lowers the summed cost of the agents party to it most is
  made, if any lowers it -- individually rational to each once the gain
  is shared out by payments -- and they replan. Passes of every offer
  repeat until one moves nothing. A step weighs at most `:chunk` (default
  256) of an offer's contracts, so an offer of many -- an M- or
  OCSM-contract's -- takes several. Each state's `:offer` is `{:type
  :from :tasks}`, `:weighed` how many of its contracts have been weighed,
  and `:concluded` whether it was settled in this step; `:moved`, when a
  contract was made, is the offer with `:moves` `{task [from to]}` and
  `:gain`. Total cost never rises.

  `:contracts` chooses the types, by default one-task, swap and cluster.
  `:max-cluster` caps a cluster's size and `:max-ocsm` the tasks an
  OCSM-contract moves. Sandholm shows that each of O, C, S and M alone
  can stop short of the best allocation, and that OCSM-contracts cannot:
  unbounded, they reach it whatever the start. That is because an
  OCSM-contract may be the whole problem -- with every task put up, the
  search is over every allocation, n^k of them for k tasks among n
  agents -- so bounded, as it must be past a handful, it can stop short
  too."
  [st & {:keys [contracts max-cluster max-ocsm chunk] :or {contracts #{:o :s :c} chunk 256}}]
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
              (let [states (->> (assoc st :pending (offers st contracts max-cluster max-ocsm) :replanned false)
                                (iterate #(some-> % (recontract-step cost chunk)))
                                rest
                                (take-while some?))]
                (when (seq states)
                  (lazy-cat states (when (some :moved states) (pass (last states)))))))]
      (cons st (cons planned (pass planned))))))
