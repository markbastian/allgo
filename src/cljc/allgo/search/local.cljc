(ns allgo.search.local
  "Local search: Russell and Norvig's chapter 4, as steps.

  Local search keeps no paths, only a current state (or a few) and a way
  to judge it, and moves to better ones. Each algorithm here is a step
  function from one search state to the next, run by `iterate` as
  `allgo.search` runs the systematic ones; the `-seq` forms are the
  states, one per move or generation, for watching, and the plain forms
  run to the end and return the best state found.

    hill climbing         steepest ascent (figure 4.2), and the section's
                          variants: sideways moves on a plateau, stochastic
                          (a random uphill move) and first-choice (the
                          first random successor that is uphill)
    random restarts       any of them, again from random states, until one
                          reaches a goal -- a composition of seqs
    simulated annealing   SIMULATED-ANNEALING: a random move, taken if uphill, and
                          if downhill only with probability e^(dE/T)
    local beam            k states, each step the best k of all their
                          successors
    genetic algorithm     figure 4.7: a population bred by tournament
                          selection, one-point recombination and mutation

  Continuous local search -- the chapter's gradient ascent, Newton-
  Raphson, the simplex method of Nelder and Mead -- is
  `allgo.numerics.optimize`, minimizing where the book maximizes.

  A problem is a map:

    :initial           a state to start from
    :value             (fn [state]) -> how good it is; higher is better
    :neighbors         (fn [state]) -> every successor
    :random-neighbor   (fn [rng state]) -> one successor at random
    :random-state      (fn [rng]) -> a fresh random state, for restarts
                       and beams
    :goal?             (fn [state]) -> whether it is a solution (optional)

  Everything random takes `:rng`, a function of no arguments returning
  numbers in [0, 1) -- `(allgo.random/rng seed)` for a repeatable run."
  (:require [allgo.random :as random]
            [allgo.search :as search]
            [clojure.math :as math]))

(defn- start [{:keys [initial value] :as problem}]
  (assoc problem :current initial :current-value (value initial) :steps 0 :done? false))

(defn- move [{:keys [value] :as search} state]
  (-> search
      (assoc :current state :current-value (value state))
      (update :steps inc)))

(defn- finish [search] (assoc search :done? true))

(defn- step-unless-done
  "`step`, unless the search is over or stands on a goal -- from which a
  sideways move or a downhill one could only walk away."
  [step]
  (fn [{:keys [goal? current] :as search}]
    (when-not (or (:done? search) (and goal? (goal? current)))
      (step search))))

;; ---------------------------------------------------------------------------
;; Hill climbing

(defn- steepest-step
  "Figure 4.2: to the best successor if it is better, else stop. With
  `:sideways`, as many as that many moves in a row to a successor no
  worse, which walks a plateau instead of stopping on it."
  [{:keys [neighbors value current current-value sideways rng flat] :as search}]
  (let [scored (map (fn [s] [s (value s)]) (neighbors current))
        best (if (seq scored) (apply max (map second scored)) ##-Inf)
        ;; Ties broken at random, so a plateau is not walked in a line.
        choice (when (seq scored)
                 (random/pick rng (vec (keep (fn [[s v]] (when (= v best) s)) scored))))
        flat (long (or flat 0))]
    (cond
      (> best current-value) (assoc (move search choice) :flat 0)
      (and (= best current-value) (< flat (long (or sideways 0))))
      (assoc (move search choice) :flat (inc flat))
      :else (finish search))))

(defn- stochastic-step
  "A move to an uphill successor chosen at random, not the steepest."
  [{:keys [neighbors value current current-value rng] :as search}]
  (let [up (vec (filter #(> (value %) current-value) (neighbors current)))]
    (if (seq up) (move search (random/pick rng up)) (finish search))))

(defn- first-choice-step
  "Successors drawn at random until one is uphill, up to `:tries` (default
  100) of them: for when a state has too many successors to list."
  [{:keys [random-neighbor value current current-value rng tries] :as search}]
  (if-let [s (some (fn [_] (let [s (random-neighbor rng current)]
                             (when (> (value s) current-value) s)))
                   (range (or tries 100)))]
    (move search s)
    (finish search)))

(defn hill-climbing-seq
  "Every state of a hill climb, one per move, ending where it stops.
  `:variant` is `:steepest` (the default), `:stochastic` or
  `:first-choice`; `:sideways` allows that many moves in a row across a
  plateau."
  [problem]
  (let [step (case (:variant problem :steepest)
               :steepest steepest-step
               :stochastic stochastic-step
               :first-choice first-choice-step)]
    (search/step-seq (step-unless-done step) (start (update problem :rng #(or % (random/rng 1)))))))

(defn random-restart-seq
  "`search-seq` -- any of these -- from the problem's initial state, then
  again from a random state each time it ends short of a goal, up to
  `:restarts` times (default 100). Each state carries `:restart`, its
  run's number. A composition of seqs: the book's random-restart hill
  climbing is this with `hill-climbing-seq`."
  [search-seq {:keys [rng random-state goal? restarts] :as problem}]
  (let [rng (or rng (random/rng 1))
        problem (assoc problem :rng rng)]
    (letfn [(runs [k p]
              (lazy-seq
               (let [states (map #(assoc % :restart k) (search-seq p))
                     end (last states)]
                 (if (or (goal? (:current end)) (>= k (long (or restarts 100))))
                   states
                   (concat states (runs (inc k) (assoc p :initial (random-state rng))))))))]
      (runs 0 problem))))

;; ---------------------------------------------------------------------------
;; Simulated annealing

(defn exponential-schedule
  "A temperature falling from `t0` by the factor `alpha` a step, zero
  once it is below `floor`: T(t) = t0 alpha^t."
  [t0 alpha floor]
  (fn [t] (let [temp (* t0 (math/pow alpha t))] (if (< temp floor) 0.0 temp))))

(defn- annealing-step
  "SIMULATED-ANNEALING: a random successor, taken if it is better, and if worse
  only with probability e^(dE/T); done when the schedule reaches zero.
  It keeps the best state seen as `:best`, which the book's version, which
  returns where it is when the temperature runs out, does not need."
  [{:keys [random-neighbor value schedule current current-value rng t best-value] :as search}]
  (let [temp (schedule t)]
    (if (zero? temp)
      (finish search)
      (let [nxt (random-neighbor rng current)
            v (value nxt)
            de (- v current-value)
            search (assoc search :t (inc (long t)) :temperature temp)]
        (if (or (pos? de) (< (double (rng)) (math/exp (/ de temp))))
          (cond-> (move search nxt)
            (> v best-value) (assoc :best nxt :best-value v))
          search)))))

(defn simulated-annealing-seq
  "Every state of a simulated annealing run, one per proposed move.
  `:schedule` maps the step number to a temperature (default
  `(exponential-schedule 20 0.995 0.001)`)."
  [problem]
  (let [{:keys [current current-value] :as s}
        (start (-> problem
                   (update :rng #(or % (random/rng 1)))
                   (update :schedule #(or % (exponential-schedule 20.0 0.995 0.001)))))]
    (search/step-seq (step-unless-done annealing-step)
                     (assoc s :t 1 :best current :best-value current-value))))

;; ---------------------------------------------------------------------------
;; Local beam search

(defn- beam-step
  "All the successors of all k states, and the best k of them; or, with
  `:stochastic?`, k of them drawn with probability rising with their
  value. Done when no successor beats the best state held."
  [{:keys [neighbors value beam k stochastic? rng] :as search}]
  (let [best-held (apply max (map second beam))
        pool (vec (distinct (for [[s] beam n (neighbors s)] [n (value n)])))
        top (apply max (map second pool))]
    (if (<= top best-held)
      (finish search)
      (let [chosen (if stochastic?
                     ;; Weighted by value, shifted to be positive.
                     (let [lo (apply min (map second pool))
                           weights (mapv (fn [[_ v]] (+ 1.0 (- v lo))) pool)
                           total (reduce + weights)]
                       (vec (repeatedly k (fn []
                                            (let [r (* total (double (rng)))]
                                              (loop [i 0 acc 0.0]
                                                (let [acc (+ acc (weights i))]
                                                  (if (or (< r acc) (= i (dec (count pool))))
                                                    (pool i)
                                                    (recur (inc i) acc)))))))))
                     (vec (take k (sort-by (comp - second) pool))))
            [s v] (apply max-key second chosen)]
        (-> search
            (assoc :beam chosen :current s :current-value v)
            (update :steps inc))))))

(defn local-beam-seq
  "Every state of a local beam search: `:k` states (default 8) started at
  random, each step replaced by the best k of all their successors --
  which, unlike k separate climbs, passes the effort to wherever the
  climbing is good. Each state holds the `:beam` of `[state value]`."
  [{:keys [random-state value k] :as problem}]
  (let [rng (or (:rng problem) (random/rng 1))
        k (or k 8)
        beam (vec (for [_ (range k) :let [s (random-state rng)]] [s (value s)]))
        [s v] (apply max-key second beam)]
    (search/step-seq (step-unless-done beam-step)
                     (assoc problem :rng rng :k k :beam beam :current s :current-value v :steps 0 :done? false))))

;; ---------------------------------------------------------------------------
;; The genetic algorithm

(defn- select
  "Figure 4.7's SELECT: the best `rho` of a uniform sample of 2 rho."
  [rng population fitness rho]
  (take rho (sort-by (comp - fitness) (repeatedly (* 2 rho) #(random/pick rng population)))))

(defn- recombine
  "One-point crossover: the first part of `x` and the rest of `y`."
  [rng x y]
  (let [c (random/below rng (inc (count x)))]
    (into (subvec x 0 c) (subvec y c))))

(defn- mutate
  "Each gene replaced, with probability `rate`, by a letter of `alphabet`."
  [rng alphabet rate individual]
  (mapv (fn [g] (if (< (double (rng)) rate) (random/pick rng alphabet) g)) individual))

(defn- generation-step
  [{:keys [population fitness alphabet mutation rng goal? generations] :as search}]
  (let [best (apply max-key fitness population)]
    (if (or (and goal? (goal? best)) (>= (long (:generation search)) (long generations)))
      (finish search)
      (let [nxt (vec (repeatedly (count population)
                                 #(let [[x y] (select rng population fitness 2)]
                                    (mutate rng alphabet mutation (recombine rng x y)))))
            b (apply max-key fitness nxt)]
        (-> search
            (assoc :population nxt :current b :current-value (fitness b))
            (update :generation inc)
            (update :steps inc))))))

(defn genetic-seq
  "Every generation of a genetic algorithm: `:population`, a vector of
  equal-length vectors over `:alphabet`, bred by `:fitness` (higher is
  better) until `:goal?` holds of its best or `:generations` (default
  1000) have passed. Each child has one-point crossover of two parents,
  each the fitter of a pair drawn at random, and then each gene mutated
  with probability `:mutation` (default 1/length)."
  [{:keys [population fitness] :as problem}]
  (let [best (apply max-key fitness population)]
    (search/step-seq (step-unless-done generation-step)
                     (-> problem
                         (update :rng #(or % (random/rng 1)))
                         (update :mutation #(or % (/ 1.0 (count (first population)))))
                         (update :generations #(or % 1000))
                         (assoc :generation 0 :steps 0 :done? false
                                :current best :current-value (fitness best))))))

(defn random-population
  "`n` individuals of `length` genes drawn from `alphabet`."
  [rng n length alphabet]
  (vec (repeatedly n (fn [] (vec (repeatedly length #(random/pick rng alphabet)))))))

;; ---------------------------------------------------------------------------
;; Answers

(defn result
  "The best state a local search's seq found: `{:state :value :steps}`,
  `:steps` counting moves, or generations, to the end."
  [states]
  (let [end (last states)
        [s v] (if (:best end) [(:best end) (:best-value end)] [(:current end) (:current-value end)])]
    {:state s :value v :steps (:steps end) :restarts (:restart end)}))

(def hill-climbing (comp result hill-climbing-seq))
(def simulated-annealing (comp result simulated-annealing-seq))
(def local-beam (comp result local-beam-seq))
(def genetic (comp result genetic-seq))
