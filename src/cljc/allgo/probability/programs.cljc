(ns allgo.probability.programs
  "Probabilistic programs: the idea of Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, chapter 15 -- a model
  written as a program that makes random choices and conditions on what
  is observed, inference running the program.

  The book's own programs for the chapter (its citation, NET-VISA, radar,
  image and Markov-letter figures) are linked from the aima-pseudocode
  index but missing from the repository, so the examples here are this
  library's: a Markov chain of letters read through noise, and an
  open-universe count -- how many aircraft are there, when the radar
  misses some and invents others? -- each checked against an exact
  answer.

  A program is `(fn [ctx])` that calls `sample!` for each random choice
  and `observe!` for each observation; `likelihood-weighting-seq` runs it
  again and again, drawing every choice from its prior and weighting each
  run by the likelihood of the observations, and estimates the
  distribution of whatever the program returns. A distribution is a map
  `{value p}`."
  (:require [allgo.random :as random]
            [clojure.math :as math]))

(defn sample!
  "A random choice from `dist`."
  [{:keys [rng]} dist]
  (random/weighted rng dist))

(defn observe!
  "Condition the run on `value` having been drawn from `dist`: its
  probability multiplies the run's weight."
  [{:keys [weight]} dist value]
  (vswap! weight * (get dist value 0.0)))

(defn likelihood-weighting-seq
  "The estimate of the distribution of `program`'s result after each run,
  `{:estimate :runs}`."
  [program rng]
  (rest (reductions (fn [{:keys [totals runs]} _]
                      (let [ctx {:rng rng :weight (volatile! 1.0)}
                            v (program ctx)
                            w @(:weight ctx)
                            totals (update totals v (fnil + 0.0) w)
                            z (reduce + (vals totals))]
                        {:totals totals :runs (inc (long runs))
                         :estimate (if (pos? z) (into {} (map (fn [[k t]] [k (/ t z)]) totals)) {})}))
                    {:totals {} :runs 0}
                    (range))))

;; ---------------------------------------------------------------------------
;; Markov letters

(defn markov-letters
  "A program writing `n` letters from a Markov chain -- `start` the first
  letter's distribution, `next` each letter's given the one before --
  each read correctly with probability `accuracy` and otherwise as any
  other letter, and conditioned on reading `seen`. It returns the letter
  at `position`."
  [{:keys [letters start next accuracy]} seen position]
  (fn [ctx]
    (let [noise (fn [l] (into {} (for [m letters] [m (if (= m l) accuracy (/ (- 1.0 accuracy) (dec (count letters))))])))
          word (reduce (fn [word i]
                         (let [l (sample! ctx (if (zero? i) start (next (peek word))))]
                           (observe! ctx (noise l) (nth seen i))
                           (conj word l)))
                       [] (range (count seen)))]
      (word position))))

;; ---------------------------------------------------------------------------
;; An open universe: how many aircraft?

(defn poisson [lambda n-max]
  (let [ps (mapv (fn [k] (* (math/exp (- lambda)) (/ (math/pow lambda k) (reduce * 1.0 (range 1 (inc k)))))) (range (inc n-max)))
        z (reduce + ps)]
    (into {} (map-indexed (fn [k p] [k (/ p z)]) ps))))

(defn aircraft-count
  "A program for an open universe: the number of aircraft is unknown,
  Poisson with mean `lambda` (up to `n-max`); each is detected with
  probability `detect`; false alarms are Poisson with mean `clutter`; and
  the radar shows `blips`. It returns the number of aircraft."
  [{:keys [lambda n-max detect clutter]} blips]
  (fn [ctx]
    (let [n (sample! ctx (poisson lambda n-max))
          detected (count (filter true? (repeatedly n #(sample! ctx {true detect false (- 1.0 detect)}))))]
      ;; The false alarms are whatever the blips the aircraft do not
      ;; explain: observed, not drawn.
      (observe! ctx (poisson clutter (+ n-max blips)) (- blips detected))
      n)))

(defn aircraft-count-exact
  "The exact posterior of the same model, summing over how many aircraft
  were detected: P(n | k) proportional to P(n) times the sum over d of
  Binomial(d; n, detect) Poisson(k - d; clutter)."
  [{:keys [lambda n-max detect clutter]} blips]
  (let [prior (poisson lambda n-max)
        fa (poisson clutter (+ n-max blips))
        choose (fn [n d] (/ (reduce * 1.0 (range (inc (- n d)) (inc n))) (reduce * 1.0 (range 1 (inc d)))))
        unnorm (into {} (for [n (range (inc n-max))]
                          [n (* (prior n)
                                (reduce + (for [d (range (inc (min n blips)))]
                                            (* (choose n d) (math/pow detect d) (math/pow (- 1.0 detect) (- n d))
                                               (get fa (- blips d) 0.0)))))]))
        z (reduce + (vals unnorm))]
    (into {} (map (fn [[n p]] [n (/ p z)]) unnorm))))
