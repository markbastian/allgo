(ns allgo.agents
  "Agent programs: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, chapter 2, in the vacuum world.

    table-driven       TABLE-DRIVEN-AGENT (figure 2.7): the whole percept
                       sequence looked up in a table of actions
    reflex vacuum      REFLEX-VACUUM-AGENT (figure 2.8): suck if dirty,
                       else move to the other square
    simple reflex      SIMPLE-REFLEX-AGENT (figure 2.10): the first
                       condition-action rule that matches the percept
    model-based        MODEL-BASED-REFLEX-AGENT (figure 2.12): rules
                       matched against a state kept up to date from
                       percepts, the last action and a model of the world

  An agent program is `(fn [memory percept]) -> [action memory']`, its
  persistent variables the memory, as for the online agents of
  `allgo.search.online`; `run-seq` runs one in the vacuum world by
  `iterate`. The world is aima-python's (MIT) TrivialVacuumEnvironment:
  squares A and B, a percept `[location status]`, and a score of +10
  for each dirty square sucked clean and -1 for each move.")

;; ---------------------------------------------------------------------------
;; The vacuum world

(defn vacuum-world
  "The agent at `location` (`:A` or `:B`), each square `:dirty` or
  `:clean` as `status` says."
  [location status]
  {:location location :status status :score 0})

(defn percept [{:keys [location status]}] [location (status location)])

(defn execute
  "The world after `action`: `:left`, `:right`, `:suck` or `:no-op`."
  [{:keys [location status] :as w} action]
  (case action
    :right (-> w (assoc :location :B) (update :score dec))
    :left (-> w (assoc :location :A) (update :score dec))
    :suck (cond-> (assoc-in w [:status location] :clean)
            (= :dirty (status location)) (update :score + 10))
    :no-op w))

(defn- run-step
  "The action carried out, and the program's next choice."
  [program {:keys [world action memory]}]
  (let [w (execute world action)
        [a m] (program memory (percept w))]
    {:world w :action a :memory m}))

(defn run-seq
  "`agent`, a `[program memory]` pair, acting in `world` for `steps` steps:
  each state `{:world :action :memory}`."
  [[program memory] world steps]
  (take (inc steps)
        (iterate (partial run-step program)
                 (let [[a m] (program memory (percept world))]
                   {:world world :action a :memory m}))))

;; ---------------------------------------------------------------------------
;; Agent programs

(defn table-driven-agent
  "TABLE-DRIVEN-AGENT: the action `table` gives for the percepts so far,
  `:no-op` when it has none."
  [table]
  [(fn [percepts p] (let [ps (conj percepts p)] [(get table ps :no-op) ps]))
   []])

(def vacuum-table
  "The table aima-python gives for the vacuum world: percept sequences of
  up to three steps."
  {[[:A :clean]] :right
   [[:A :dirty]] :suck
   [[:B :clean]] :left
   [[:B :dirty]] :suck
   [[:A :dirty] [:A :clean]] :right
   [[:A :clean] [:B :dirty]] :suck
   [[:B :clean] [:A :dirty]] :suck
   [[:B :dirty] [:B :clean]] :left
   [[:A :dirty] [:A :clean] [:B :dirty]] :suck
   [[:B :dirty] [:B :clean] [:A :dirty]] :suck})

(def reflex-vacuum-agent
  "REFLEX-VACUUM-AGENT."
  [(fn [m [location status]]
     [(cond (= status :dirty) :suck (= location :A) :right :else :left) m])
   nil])

(defn simple-reflex-agent
  "SIMPLE-REFLEX-AGENT: `interpret` turns a percept into a state, and the
  first of `rules`, `[condition action]` pairs, whose condition holds of
  it gives the action."
  [interpret rules]
  [(fn [m percept]
     (let [state (interpret percept)]
       [(some (fn [[condition action]] (when (condition state) action)) rules) m]))
   nil])

(defn model-based-reflex-agent
  "MODEL-BASED-REFLEX-AGENT: `update-state` folds the last action and a
  percept into the agent's idea of the world, and the first of `rules`
  whose condition holds of that state gives the action."
  [update-state rules initial-state]
  [(fn [{:keys [state action]} percept]
     (let [state (update-state state action percept)
           a (some (fn [[condition act]] (when (condition state) act)) rules)]
       [a {:state state :action a}]))
   {:state initial-state :action nil}])

(def model-based-vacuum-agent
  "The model-based reflex agent for the vacuum world: it remembers what it
  has seen of each square, so once both are clean it stops, where the
  simple reflex agent shuttles back and forth for ever."
  (model-based-reflex-agent
   (fn [state _ [location status]] (assoc state location status :at location))
   [[(fn [s] (= :dirty (s (:at s)))) :suck]
    [(fn [s] (= :clean (s :A) (s :B))) :no-op]
    [(fn [s] (= :A (:at s))) :right]
    [(constantly true) :left]]
   {}))
