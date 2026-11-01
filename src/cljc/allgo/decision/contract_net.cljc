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
;; Recontracting

(defn- recontract-step
  [{:keys [offers agents tasks routes] :as st}]
  (when-let [[holder task-id] (first offers)]
    (let [st (assoc st :offers (rest offers) :offer [holder task-id] :moved nil)]
      (if-not (some #{task-id} (routes holder))
        st
        (let [without (assoc-in st [:routes holder] (vec (remove #{task-id} (routes holder))))
              saving (- (agent-cost st holder) (agent-cost without holder))
              bids (into {} (for [a (keys agents)
                                  :when (not= a holder)
                                  :let [c (bid without a (tasks task-id))]
                                  :when c]
                              [a c]))
              [best c] (when (seq bids) (apply min-key val (sort-by key bids)))]
          (if (and best (< (double c) (- saving 1e-9)))
            (-> without (take-task best task-id) (assoc-in [:awarded task-id] best)
                (assoc :moved [task-id holder best (- saving c)]))
            st))))))

(defn recontract-seq
  "Recontracting from the allocation `st`, an offer a step: each agent in
  turn offers each task it holds, and it moves to the cheapest other
  bidder when that bid is below what dropping it saves the holder.
  Passes repeat until one moves nothing. Each state's `:offer` is the
  `[holder task]` just offered and `:moved` `[task from to gain]` when it
  moved. Total cost never rises."
  [st]
  (letfn [(offers-of [st] (for [[a route] (sort-by key (:routes st)) t route] [a t]))
          (pass [st]
            (let [states (take-while some? (iterate #(some-> % recontract-step) (assoc st :offers (offers-of st) :moved nil)))
                  states (rest states)]
              (when (seq states)
                (lazy-cat states
                          (when (some :moved states)
                            (pass (last states)))))))]
    (cons st (pass st))))
