(ns allgo.search.online
  "Online search: Russell and Norvig's chapter 4, its last section. An online agent does
  not know the state space in advance; it learns what an action does only
  by doing it, so it must act to explore, and it can only expand the
  state it is standing in.

    online DFS   figure 4.20's ONLINE-DFS-AGENT: depth-first, physically
                 -- it walks back to backtrack, which needs every action
                 to be undoable
    LRTA*        figure 4.23's LRTA*-AGENT: moves toward the neighbor
                 that looks cheapest by its cost plus the estimate from
                 there, and raises its estimate of each state it leaves
                 to what it learned, so that a dead end, once tried, stops
                 looking promising

  An agent is a pure function `(fn [memory percept]) -> [action memory']`
  -- the book's persistent tables are its memory, passed back in -- and
  `run-seq` runs it in an environment, a step being one percept and one
  action, by `iterate`. A problem is a map:

    :start      where the agent begins
    :actions    (fn [state]) -> the actions it may take there
    :result     (fn [state action]) -> where the action really leads; the
                environment knows this, and the agent learns it
    :goal?      (fn [state]) -> whether it has arrived
    :cost       (fn [state action state']) -> the cost of a move
                (default 1), for LRTA*
    :h          (fn [state]) -> an estimate of the cost to a goal
                (default 0), for LRTA*")

;; ---------------------------------------------------------------------------
;; Online depth-first search

(defn online-dfs-agent
  "The online DFS agent for `problem`'s `:actions` and `:goal?`, and its
  empty memory: `[agent memory]`.

  One change to the figure: a move that leaves the agent where it was --
  into a wall -- is not recorded as a way back to where it was, which the
  figure would do, and later walk into the wall again to 'backtrack'."
  [{:keys [actions goal?]}]
  [(fn [{:keys [result untried unbacktracked s a] :as memory} s']
     (if (goal? s')
       [:stop memory]
       (let [untried (if (contains? untried s') untried (assoc untried s' (vec (actions s'))))
             new? (and (some? s) (not= s' (get result [s a])))
             result (if new? (assoc result [s a] s') result)
             unbacktracked (if (and new? (not= s s'))
                             (update unbacktracked s' #(cons s %))
                             unbacktracked)
             memory (assoc memory :result result :untried untried :unbacktracked unbacktracked)]
         (if (empty? (untried s'))
           (if (empty? (unbacktracked s'))
             [:stop memory]
             ;; Back the way it came: an action already found to lead to
             ;; the state it backtracks to.
             (let [back (first (unbacktracked s'))
                   b (first (filter #(= back (get result [s' %])) (actions s')))]
               [b (assoc memory :unbacktracked (update unbacktracked s' rest) :s s' :a b)]))
           (let [b (peek (untried s'))]
             [b (assoc memory :untried (update untried s' pop) :s s' :a b)])))))
   {:result {} :untried {} :unbacktracked {} :s nil :a nil}])

;; ---------------------------------------------------------------------------
;; LRTA*

(defn lrta-star-agent
  "The LRTA* agent for `problem`, and its empty memory: `[agent memory]`.
  Its memory's `:H` is what it has learned of each state's cost to go."
  [{:keys [actions goal? cost h]}]
  (let [cost (or cost (constantly 1))
        h (or h (constantly 0))
        ;; LRTA*-COST: an action not yet tried is assumed to lead
        ;; somewhere as good as the estimate for here -- optimism under
        ;; uncertainty.
        lrta-cost (fn [result H s b]
                    (if-let [s' (get result [s b])]
                      (+ (cost s b s') (get H s' (h s')))
                      (h s)))]
    [(fn [{:keys [result H s a] :as memory} s']
       (if (goal? s')
         [:stop memory]
         (let [H (if (contains? H s') H (assoc H s' (h s')))
               result (if (some? s) (assoc result [s a] s') result)
               H (if (some? s)
                   (assoc H s (apply min (map #(lrta-cost result H s %) (actions s))))
                   H)
               b (apply min-key #(lrta-cost result H s' %) (reverse (actions s')))]
           [b (assoc memory :result result :H H :s s' :a b)])))
     {:result {} :H {} :s nil :a nil}]))

;; ---------------------------------------------------------------------------
;; Running an agent

(defn run-seq
  "The agent made by `make-agent` (`online-dfs-agent` or `lrta-star-agent`)
  acting in `problem`, a step at a time: each state `{:state :action
  :memory :steps :cost}` is where it stands, the action it chose there,
  its memory after choosing, and the moves and cost so far. It ends when
  the agent stops -- at a goal, or, for online DFS, having exhausted
  everything it could reach -- or after `limit` steps."
  [make-agent {:keys [start result cost] :as problem} limit]
  (let [[agent memory] (make-agent problem)
        cost (or cost (constantly 1))
        [a m] (agent memory start)]
    (take-while some?
                (iterate (fn [{:keys [state action memory steps] :as st}]
                           (when (and st (not= :stop action) (< (long steps) (long limit)))
                             (let [s' (result state action)
                                   [a m] (agent memory s')]
                               {:state s' :action a :memory m :steps (inc (long steps))
                                :cost (+ (:cost st) (cost state action s'))})))
                         {:state start :action a :memory m :steps 0 :cost 0}))))
