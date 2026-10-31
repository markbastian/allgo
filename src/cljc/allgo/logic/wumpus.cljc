(ns allgo.logic.wumpus
  "The wumpus world, and the hybrid agent that explores it: Russell and
  Norvig, *Artificial Intelligence: A Modern Approach*, 4th edition,
  chapter 7 (HYBRID-WUMPUS-AGENT, figure 7.20).

  The world is the chapter's: a 4 by 4 cave, squares `[x y]` from `[1 1]`,
  where the agent starts facing east. Each square but `[1 1]` holds a pit
  with probability 0.2; the wumpus and the gold are each in one square
  other than `[1 1]`, chosen at random. The agent perceives a stench next
  to the wumpus, a breeze next to a pit, glitter where the gold is, a bump
  when it walks into a wall and a scream when its one arrow kills the
  wumpus. Scores are as aima-python (MIT) has them: -1 an action, +1000
  for climbing out of `[1 1]` with the gold, -1000 for walking into a pit
  or a live wumpus. The book also charges the arrow extra; aima-python
  leaves the arrow out and aima-java keeps no score, so that charge is
  `:arrow-cost` here, default -1 like any other action.

  The agent is the figure's, with one simplification of my own: it knows
  where it is and which way it faces by keeping count, so the knowledge
  base holds only what it cannot see -- the pits and the wumpus. It
  starts as the atemporal 'wumpus physics' (a breeze in a square iff a
  pit next to it; a stench iff the wumpus next to it; exactly one wumpus;
  none in `[1 1]`), and each percept adds what was felt where. Whether a
  square is safe is asked of it -- does it entail no pit and no live
  wumpus there? -- by DPLL; routes are planned by A* (`allgo.search`) over
  the squares known safe. A shot that brings no scream tells the agent the
  squares the arrow crossed hold no wumpus, which the book's temporal
  axioms would infer and which is told here directly."
  (:require [allgo.logic.propositional :as pl]
            [allgo.random :as random]
            [allgo.search :as search]))

;; ---------------------------------------------------------------------------
;; The world

(def size 4)

(def ^:private squares (vec (for [y (range 1 (inc size)) x (range 1 (inc size))] [x y])))

(defn neighbors [[x y]]
  (filterv (fn [[a b]] (and (<= 1 a size) (<= 1 b size)))
           [[(inc x) y] [x (inc y)] [(dec x) y] [x (dec y)]]))

(defn world
  "A random cave from `seed`."
  [seed]
  (let [rng (random/rng seed 29 0 0)
        others (vec (remove #{[1 1]} squares))]
    {:pits (set (filter (fn [_] (< (double (rng)) 0.2)) others))
     :wumpus (random/pick rng others)
     :gold (random/pick rng others)}))

(def ^:private headings [:east :north :west :south])
(def ^:private step-of {:east [1 0] :north [0 1] :west [-1 0] :south [0 -1]})

(defn- ahead [[x y] dir] (let [[dx dy] (step-of dir)] [(+ x dx) (+ y dy)]))

(def ^:private heading-index (zipmap headings (range)))

(defn- turn [dir by] (headings (mod (+ (long (heading-index dir)) (long by)) 4)))

(defn start
  "The state of play in `world`: the agent at `[1 1]` facing east, with
  its arrow."
  ([w] (start w {}))
  ([w {:keys [arrow-cost] :or {arrow-cost -1}}]
   {:world w :pos [1 1] :dir :east :arrow? true :gold? false
    :wumpus-alive? true :alive? true :out? false :score 0
    :bump? false :scream? false :arrow-cost arrow-cost}))

(defn percept
  "What the agent senses in `state`: `{:stench :breeze :glitter :bump
  :scream}`."
  [{:keys [world pos bump? scream? gold?]}]
  (let [{:keys [pits wumpus gold]} world
        near (set (neighbors pos))]
    ;; Next to the wumpus only, as the physics says -- not in its own
    ;; square, which the agent can stand in once it is dead.
    {:stench (contains? near wumpus)
     :breeze (boolean (some pits near))
     :glitter (and (= pos gold) (not gold?))
     :bump bump?
     :scream scream?}))

(defn execute
  "The state after the agent does `action` in `state`."
  [{:keys [world pos dir arrow? wumpus-alive?] :as state} action]
  (let [state (assoc state :bump? false :scream? false)
        state (update state :score + (if (= action :shoot) (:arrow-cost state) -1))]
    (case action
      :left (assoc state :dir (turn dir 1))
      :right (assoc state :dir (turn dir -1))
      :forward (let [p (ahead pos dir)
                     [x y] p]
                 (if (and (<= 1 x size) (<= 1 y size))
                   (let [dead? (or (contains? (:pits world) p) (and wumpus-alive? (= p (:wumpus world))))]
                     (cond-> (assoc state :pos p)
                       dead? (-> (assoc :alive? false) (update :score - 1000))))
                   (assoc state :bump? true)))
      :grab (if (= pos (:gold world)) (assoc state :gold? true) state)
      :shoot (if arrow?
               (let [line (take-while (fn [[x y]] (and (<= 1 x size) (<= 1 y size)))
                                      (rest (iterate #(ahead % dir) pos)))
                     hit? (and wumpus-alive? (some #{(:wumpus world)} line))]
                 (cond-> (assoc state :arrow? false)
                   hit? (assoc :wumpus-alive? false :scream? true)))
               state)
      :climb (if (= pos [1 1])
               (cond-> (assoc state :out? true)
                 (:gold? state) (update :score + 1000))
               state))))

(defn over? [{:keys [alive? out?]}] (or (not alive?) out?))

;; ---------------------------------------------------------------------------
;; The knowledge base

(defn- sym [prefix [x y]] (keyword (str prefix x "," y)))
(defn pit [sq] (sym "P" sq))
(defn wumpus [sq] (sym "W" sq))
(defn breeze [sq] (sym "B" sq))
(defn stench [sq] (sym "S" sq))

(def physics
  "The atemporal wumpus physics, as clauses."
  (pl/cnf (into [:and [:not (pit [1 1])] [:not (wumpus [1 1])]
                 (into [:or] (map wumpus squares))]
                (concat
                 (for [sq squares]
                   [:iff (breeze sq) (into [:or] (map pit (neighbors sq)))])
                 (for [sq squares]
                   [:iff (stench sq) (into [:or] (map wumpus (neighbors sq)))])
                 ;; At most one wumpus.
                 (for [i (range (count squares)) j (range (inc i) (count squares))]
                   [:or [:not (wumpus (squares i))] [:not (wumpus (squares j))]])))))

(defn- ask
  "Whether `clauses` entail the literal or sentence `q`, counting queries."
  [memory q]
  (let [clauses (into (:kb memory) (pl/cnf [:not q]))]
    [(nil? (pl/dpll clauses)) (update memory :queries inc)]))

(defn- ok-sentence [memory sq]
  (if (:wumpus-dead? memory)
    [:not (pit sq)]
    [:and [:not (pit sq)] [:not (wumpus sq)]]))

;; ---------------------------------------------------------------------------
;; Routes

(defn- route
  "PLAN-ROUTE: the actions from `pos` facing `dir` to any square
  satisfying `goal?` -- or, with `aim`, to a pose where `(aim sq dir)`
  holds -- through squares in `allowed`, by A*; nil if there is none."
  [pos dir allowed goal? aim]
  (let [goal-pose? (fn [[sq d]] (if aim (aim sq d) (goal? sq)))
        moves (fn [[sq d]]
                (let [f (ahead sq d)]
                  (cond-> [[sq (turn d 1)] [sq (turn d -1)]]
                    (contains? allowed f) (conj [f d]))))
        path (search/a-star {:start [pos dir] :goal nil :goal? goal-pose?
                             :neighbors moves
                             :heuristic (constantly 0)})]
    (when path
      (vec (for [[[a da] [b db]] (partition 2 1 path)]
             (cond (not= a b) :forward
                   (= db (turn da 1)) :left
                   :else :right))))))

(defn- in-line?
  "Whether an arrow shot from `sq` facing `dir` would cross `target`."
  [sq dir target]
  (some #{target} (take size (rest (iterate #(ahead % dir) sq)))))

;; ---------------------------------------------------------------------------
;; The agent

(defn- tell [memory sentence] (update memory :kb into (pl/cnf sentence)))

(defn- update-knowledge
  "What the agent can now prove about the squares it has not proved
  anything about: safe, a pit, the wumpus. Proofs are kept -- the KB only
  grows, so what it entails stays entailed -- and only squares next to a
  visited one are asked about, the only ones a percept can bear on."
  [memory]
  (let [frontier (set (mapcat neighbors (:visited memory)))
        open (remove (fn [sq] (or ((:safe memory) sq)
                                  (and ((:pits memory) sq) ((:no-wumpus memory) sq))))
                     frontier)]
    (reduce (fn [m sq]
              (let [[safe? m] (ask m (ok-sentence m sq))
                    m (if safe? (update m :safe conj sq) m)
                    [pit? m] (if (or safe? ((:pits m) sq)) [false m] (ask m (pit sq)))
                    m (if pit? (update m :pits conj sq) m)
                    [no-w? m] (if ((:no-wumpus m) sq) [true m] (ask m [:not (wumpus sq)]))
                    m (if no-w? (update m :no-wumpus conj sq) m)]
                m))
            memory open)))

(defn- after-shot
  "The memory told what a shot without a scream proves: no wumpus on the
  arrow's line."
  [memory scream?]
  (if-let [[sq d] (and (not scream?) (:shot-from memory))]
    (let [line (take-while (fn [[x y]] (and (<= 1 x size) (<= 1 y size))) (rest (iterate #(ahead % d) sq)))]
      (-> (reduce (fn [m s] (tell m [:not (wumpus s)])) memory line)
          (dissoc :shot-from)))
    (dissoc memory :shot-from)))

(defn hybrid-agent
  "The hybrid wumpus agent as `(fn [memory percept]) -> [action memory']`,
  and its starting memory: `[agent memory]`. The memory's `:safe`,
  `:pits`, `:no-wumpus` and `:visited` are what it knows, `:plan` what it
  means to do, and `:queries` how many questions it has put to DPLL."
  []
  [(fn [memory {:keys [glitter scream] felt-breeze :breeze felt-stench :stench}]
     (let [{:keys [pos dir t] :as memory} (after-shot memory scream)
           memory (-> memory
                      (update :visited conj pos)
                      (update :safe conj pos)
                      (tell (if felt-breeze (breeze pos) [:not (breeze pos)]))
                      (tell (if felt-stench (stench pos) [:not (stench pos)]))
                      (cond-> scream (assoc :wumpus-dead? true)))
           memory (update-knowledge memory)
           {:keys [safe visited]} memory
           plan (:plan memory)
           plan (if glitter
                  (into [:grab] (conj (or (route pos dir safe #{[1 1]} nil) []) :climb))
                  plan)
           unvisited (fn [sq] (not (visited sq)))
           plan (or (seq plan)
                    (route pos dir safe #(and (safe %) (unvisited %)) nil))
           ;; Where the wumpus might be: no proof it is not there.
           possible-wumpus (when (and (empty? plan) (:arrow? memory) (not (:wumpus-dead? memory)))
                             (set (remove (:no-wumpus memory) (remove #{[1 1]} squares))))
           plan (or (seq plan)
                    (when (seq possible-wumpus)
                      (when-let [r (route pos dir safe nil
                                          (fn [sq d] (and (safe sq) (some #(in-line? sq d %) possible-wumpus))))]
                        (conj r :shoot))))
           ;; No choice but a risk: a square not proved unsafe.
           plan (or (seq plan)
                    (let [frontier (set (filter unvisited (mapcat neighbors visited)))
                          risky (remove (:pits memory) frontier)
                          ;; Not proved to hold the live wumpus either.
                          risky (filter (fn [sq] (or (:wumpus-dead? memory)
                                                     (not= #{sq} (set (remove (:no-wumpus memory) squares)))))
                                        risky)]
                      (route pos dir (into safe risky) (set risky) nil)))
           plan (or (seq plan)
                    (conj (or (route pos dir safe #{[1 1]} nil) []) :climb))
           [action & more] plan
           memory (cond-> (assoc memory :plan (vec more) :t (inc (long t)) :last action)
                    (= action :shoot) (assoc :arrow? false :shot-from [pos dir])
                    (= action :left) (assoc :dir (turn dir 1))
                    (= action :right) (assoc :dir (turn dir -1))
                    (= action :forward) (assoc :pos (ahead pos dir)))]
       [action memory]))
   {:kb physics :t 0 :pos [1 1] :dir :east :arrow? true :plan []
    :safe #{[1 1]} :pits #{} :no-wumpus #{[1 1]} :visited #{} :queries 0}])

(defn run-seq
  "The hybrid agent playing `world`, a step at a time: each state
  `{:state :memory :action :percept}` is the world after the action, what
  the agent knows, what it chose and what it felt. Ends when the agent
  climbs out or dies, or after `limit` actions (default 300)."
  ([w] (run-seq w {}))
  ([w {:keys [limit] :or {limit 300} :as opts}]
   (let [[agent memory] (hybrid-agent)
         s0 (start w opts)
         p0 (percept s0)
         [a m] (agent memory p0)]
     (take-while some?
                 (iterate (fn [{:keys [state action memory steps] :as st}]
                            (when (and st (not (over? state)) (< (long steps) (long limit)))
                              (let [s (execute state action)
                                    p (percept s)]
                                (if (over? s)
                                  (assoc st :state s :action nil :percept p :steps (inc (long steps)))
                                  (let [[a m] (agent memory p)]
                                    {:state s :memory m :action a :percept p :steps (inc (long steps))})))))
                          {:state s0 :memory m :action a :percept p0 :steps 0})))))
