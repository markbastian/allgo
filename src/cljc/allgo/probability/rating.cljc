(ns allgo.probability.rating
  "Bayesian ratings of players' skill: the relational model of Russell and
  Norvig, *Artificial Intelligence: A Modern Approach*, 4th edition,
  section 15.1.2.

    Skill(i)            ~ N(mu, sigma^2)
    Performance(i, g)   ~ N(Skill(i), beta^2)
    Win(i, j, g)        = Performance(i, g) > Performance(j, g)

  and for teams, a team's performance the sum of its players'. A rating
  is the posterior over a player's skill, `{:mu :sigma}`.

  The book leaves inference to an engine. Here it is the one Microsoft's
  TrueSkill uses for this model (Herbrich, Minka and Graepel, \"TrueSkill:
  A Bayesian skill rating system\", NIPS 2006): assumed density filtering,
  each game folded in as it is played and the posterior kept Gaussian by
  matching its first two moments. The difference of two teams'
  performances is Gaussian, mean the difference of their skills' sums and
  variance c^2 = (players) beta^2 + the sum of their skill variances;
  learning that it was positive truncates it, and with t its mean over c,

    v(t) = phi(t) / Phi(t),     w(t) = v(t) (v(t) + t)

  are the truncated Gaussian's shift and shrinkage. Each winner's mean
  moves up by sigma_i^2 v / c, each loser's down as far, and each
  variance shrinks by the factor 1 - sigma_i^2 w / c^2. For one game
  between players with Gaussian ratings this is the exact posterior's
  mean and variance; over many it is an approximation, which TrueSkill
  then improves by passes of expectation propagation, left out here.

  `:tau` adds skill drift before each game -- players change -- so that a
  rating never grows so certain it stops learning."
  (:require [allgo.numerics.special :as special]
            [clojure.math :as math]))

(defn rating
  "A rating: mean `mu`, standard deviation `sigma`."
  [mu sigma]
  {:mu (double mu) :sigma (double sigma)})

(defn- team-moments
  "A team's skill sum, as mean and variance."
  [ratings team]
  [(reduce + (map #(:mu (ratings %)) team))
   (reduce + (map #(math/pow (:sigma (ratings %)) 2) team))])

(defn win-probability
  "The probability that team `a` (a collection of player ids) outperforms
  team `b`, given `ratings` and the performance spread `beta`: Phi of the
  difference of their means over c."
  [ratings a b beta]
  (let [[ma va] (team-moments ratings a)
        [mb vb] (team-moments ratings b)
        c (math/sqrt (+ (* (+ (count a) (count b)) beta beta) va vb))]
    (special/normal-cdf (/ (- ma mb) c))))

(defn- truncation
  "v(t) and w(t), the shift and shrinkage of a Gaussian truncated below
  zero, at its mean over its spread `t`. Far in the losing tail,
  phi/Phi is computed as its limit, -t."
  [t]
  (let [p (special/normal-cdf t)
        v (if (< p 1e-300) (- t) (/ (special/normal-pdf t) p))]
    [v (* v (+ v t))]))

(defn update-ratings
  "`ratings` after a game won by team `winner` over team `loser`: each
  player's rating moved by the truncated Gaussian's moments. `:beta` the
  performance spread, `:tau` the skill drift added first (default 0)."
  [ratings {:keys [winner loser]} {:keys [beta tau] :or {tau 0.0}}]
  (let [drift (fn [r] (update r :sigma #(math/sqrt (+ (* % %) (* tau tau)))))
        ratings (reduce (fn [rs id] (update rs id drift)) ratings (concat winner loser))
        [mw vw] (team-moments ratings winner)
        [ml vl] (team-moments ratings loser)
        c2 (+ (* (+ (count winner) (count loser)) beta beta) vw vl)
        c (math/sqrt c2)
        [v w] (truncation (/ (- mw ml) c))
        move (fn [sign] (fn [{:keys [mu sigma]}]
                          (let [s2 (* sigma sigma)]
                            {:mu (+ mu (* sign (/ s2 c) v))
                             :sigma (math/sqrt (* s2 (max 0.0 (- 1.0 (* (/ s2 c2) w)))))})))]
    (as-> ratings rs
      (reduce (fn [rs id] (update rs id (move 1.0))) rs winner)
      (reduce (fn [rs id] (update rs id (move -1.0))) rs loser))))

(defn- rating-step [opts ratings game] (update-ratings ratings game opts))

(defn rating-seq
  "The ratings after each of `games` in turn, from `prior`, a map of
  player to rating, the prior first. Each game `{:winner [ids] :loser
  [ids]}`; options as `update-ratings`'s."
  [prior games opts]
  (reductions (partial rating-step opts) prior games))

(defn ranking
  "The players of `ratings` from most to least skilled by a conservative
  estimate, mu - 3 sigma: a player with few games ranks below one equally
  good on record but better known."
  [ratings]
  (map key (sort-by (fn [[_ {:keys [mu sigma]}]] (- (- mu (* 3.0 sigma)))) ratings)))
