(ns allgo.decision.assistance
  "Assistance games: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, sections 16.7.2 and 18.2.5 -- a robot that is
  to serve a human's preferences without knowing them.

  The human, Harriet, knows her preferences theta; the robot, Robbie, has
  only a prior P(theta); and both play for the same payoff, Harriet's.
  Robbie's uncertainty is what makes him useful and safe: Harriet's
  choices tell him about theta, and she, knowing he watches, can choose so
  as to tell him. (Introduced as cooperative inverse reinforcement
  learning by Hadfield-Menell et al., NeurIPS 2016.)

    off-switch game   Robbie can act now, switch himself off, or defer:
                      say what he plans and let Harriet switch him off.
                      Deferring is worth the expected payoff with its
                      negative part cut away -- never less than acting,
                      and more whenever he is unsure. A cost for her time,
                      or a chance that she errs, can make acting better
    paperclip game    Harriet makes a few paperclips or staples, then
                      Robbie makes many; her choice is a signal. Myopic
                      best response, a step a player, finds the
                      equilibrium: a code in which she makes one of each
                      not only when indifferent but whenever close to it,
                      and he reads that as 'make 50 of each'

  A prior over a number is given as samples on a grid, `[[value
  weight] ...]`, the weights summing to 1.")

(defn uniform-prior
  "The uniform distribution over [`lo`, `hi`] as `n` equally weighted
  points, the ends included -- so an odd `n` has the midpoint, where a
  choice between two things worth the same can be exactly indifferent."
  [lo hi n]
  (vec (for [i (range n)] [(+ lo (* (- hi lo) (/ i (double (dec n))))) (/ 1.0 n)])))

;; ---------------------------------------------------------------------------
;; The off-switch game

(defn off-switch
  "The values to Harriet of Robbie's three choices, given his `prior` over
  her utility u for his proposed action: `:act`, E[u]; `:off`, 0; and
  `:defer`, letting her decide -- she lets him go ahead when u > 0 and
  switches him off otherwise, so E[max(u, 0)], less the `cost` of her
  time, and, with probability `error`, she decides the wrong way (default
  0, 0). `{:act :off :defer :best}`."
  ([prior] (off-switch prior 0.0 0.0))
  ([prior cost error]
   (let [act (reduce + (map (fn [[u w]] (* w u)) prior))
         defer (- (reduce + (map (fn [[u w]] (* w (if (pos? u) (* (- 1.0 error) u) (* error u)))) prior)) cost)
         values {:act act :off 0.0 :defer defer}]
     (assoc values :best (key (apply max-key val values))))))

;; ---------------------------------------------------------------------------
;; The paperclip game

(def paperclip-game
  "Section 18.2.5's game: Harriet makes two paperclips, one of each, or
  two staples; then Robbie makes 90 paperclips, 50 of each, or 90
  staples. With theta the value of a paperclip and 1 - theta of a staple,
  p paperclips and s staples are worth p theta + s (1 - theta)."
  {:human [[2 0] [1 1] [0 2]]
   :robot [[90 0] [50 50] [0 90]]})

(defn payoff
  "What `[p s]` paperclips and staples are worth to a Harriet who values
  a paperclip at `theta`."
  [theta [p s]]
  (+ (* p theta) (* s (- 1.0 theta))))

(defn- best-choice
  "The choice of `choices` that `value` rates highest; among those within
  1e-12 of it, the most balanced -- fewest more of one thing than the
  other -- as the section's greedy Harriet makes one of each when she is
  exactly indifferent."
  [value choices]
  (let [top (reduce max (map value choices))]
    (apply min-key (fn [[p s]] (abs (- p s))) (filter #(>= (value %) (- top 1e-12)) choices))))

(defn greedy-human
  "Harriet's strategy were she on her own: at each theta of the `prior`,
  the choice worth most to her by itself. A map from theta to choice."
  [{:keys [human]} prior]
  (into {} (for [[theta] prior] [theta (best-choice #(payoff theta %) human)])))

(defn robot-response
  "Robbie's best reply to Harriet's `strategy` (theta to choice): for each
  choice she might make, his belief about theta is the prior restricted
  to the thetas at which she makes it, and he makes what is worth most
  under that belief. A map from her choice to his."
  [{:keys [human robot]} prior strategy]
  (into {} (for [h human
                 :let [region (filter (fn [[theta]] (= h (strategy theta))) prior)]
                 :when (seq region)]
             (let [z (reduce + (map second region))
                   mean (/ (reduce + (map (fn [[theta w]] (* theta w)) region)) z)]
               ;; the payoff is linear in theta, so the posterior mean decides
               [h (best-choice #(payoff mean %) robot)]))))

(defn human-response
  "Harriet's best reply to Robbie's `response` (her choice to his): at
  each theta, the choice whose value with his reply added is greatest.
  A choice he has no reply to yet counts as followed by nothing."
  [{:keys [human]} prior response]
  (into {} (for [[theta] prior]
             [theta (best-choice (fn [h] (+ (payoff theta h) (if-let [r (response h)] (payoff theta r) 0.0))) human)])))

(defn- best-response-step
  "Robbie's best reply to Harriet's strategy, then hers to his."
  [game prior {:keys [human]}]
  (let [robot (robot-response game prior human)]
    {:human (human-response game prior robot) :robot robot}))

(defn best-response-seq
  "Myopic best response, as the section solves the game: from Harriet's
  greedy strategy, Robbie's best reply to it, then hers to his, and so
  on. Each state `{:human :robot}`, the strategies; it has reached an
  equilibrium when a step changes neither."
  [game prior]
  (iterate (partial best-response-step game prior)
           (let [human (greedy-human game prior)]
             {:human human :robot (robot-response game prior human)})))

(defn equilibrium
  "The first state of `best-response-seq` that its next repeats."
  [game prior]
  (some (fn [[a b]] (when (= a b) a)) (partition 2 1 (best-response-seq game prior))))

(defn signal-range
  "The thetas at which Harriet's `strategy` makes `choice`: `[lo hi]`."
  [strategy choice]
  (let [ts (sort (keep (fn [[theta h]] (when (= h choice) theta)) strategy))]
    (when (seq ts) [(first ts) (last ts)])))

(defn threshold
  "Where Harriet is indifferent between two plans, each `[[p s] ...]`
  summed: the theta at which their payoffs are equal, the payoff being
  linear in theta. Nil if they are parallel."
  [plan-a plan-b]
  (let [[pa sa] (apply map + plan-a) [pb sb] (apply map + plan-b)
        ;; p theta + s (1 - theta) = s + (p - s) theta
        slope (- (- pa sa) (- pb sb))]
    (when-not (zero? slope) (/ (- sb sa) (double slope)))))
