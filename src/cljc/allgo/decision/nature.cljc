(ns allgo.decision.nature
  "Games against nature: LaValle, *Planning Algorithms*, Cambridge
  University Press, 2006, section 9.2 -- choosing an action when the
  outcome also depends on something nature chooses, which the decision
  maker cannot see.

  A problem is a cost table: `costs` maps each action to a map from each
  of nature's actions to the cost. How to choose depends on what is
  assumed about nature:

    nondeterministic  nothing: the worst case of each action, minimized --
                      the minimax decision
    probabilistic     a distribution over nature's actions: the expected
                      cost, minimized
    regret            the cost of an action against the best that could
                      have been done had nature's choice been known;
                      Savage's minimax regret takes the worst regret, and
                      is less timid than minimax when one bad case
                      dominates
    Hurwicz           a weighted mean of each action's best and worst
                      cases, the weight `alpha` an index of pessimism

  `observe` adds what section 9.2.3 adds: an observation that narrows
  nature's possible choices before the decision.")

(defn worst-case
  "The minimax action -- least worst-case cost -- and that cost: `[action
  cost]`."
  [costs]
  (apply min-key second (map (fn [[a row]] [a (reduce max (vals row))]) costs)))

(defn expected
  "The action of least expected cost under `p`, nature's distribution (a
  map from its actions to probabilities), and that cost."
  [costs p]
  (apply min-key second (map (fn [[a row]] [a (reduce + (map (fn [[n c]] (* c (p n))) row))]) costs)))

(defn regrets
  "Each action's regret against each of nature's: its cost less the best
  cost any action has against that choice of nature's."
  [costs]
  (let [natures (keys (val (first costs)))
        best (into {} (for [n natures] [n (reduce min (map #(get % n) (vals costs)))]))]
    (into {} (for [[a row] costs] [a (into {} (for [[n c] row] [n (- c (best n))]))]))))

(defn minimax-regret
  "Savage's choice: the action whose worst regret is least, and that
  regret."
  [costs]
  (worst-case (regrets costs)))

(defn hurwicz
  "Hurwicz's choice: the action minimizing alpha times its worst cost
  plus 1 - alpha times its best, `alpha` in [0, 1]; 1 is minimax, 0 the
  optimist's."
  [costs alpha]
  (apply min-key second (map (fn [[a row]] (let [cs (vals row)] [a (+ (* alpha (reduce max cs)) (* (- 1 alpha) (reduce min cs)))])) costs)))

(defn observe
  "The problem left once an observation rules out all but the `possible`
  choices of nature's."
  [costs possible]
  (into {} (for [[a row] costs] [a (select-keys row possible)])))
