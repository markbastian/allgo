(ns allgo.robotics.randomization
  "Domain randomization: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, section 26.8's sim-to-real transfer --
  \"train policies that will work with a variety of models by sampling
  different parameters in the simulations\".

  The robot is a mass to be brought to rest at the origin, pushed through
  an actuator that answers a few steps late and can push only so hard, by
  a policy u = -k1 x - k2 v. Policy search (section 22.5) picks the
  gains: Nelder-Mead on the cost of rolling the policy out, averaged over
  the simulations it is given. Searched in one simulation of the nominal
  robot, the gains come out as aggressive as that robot allows -- and a
  real robot that is lighter than modeled, the delay the same, overshoots
  and shakes itself apart. Searched over simulations whose mass and
  friction are drawn at random, the gains are gentler: a little worse on
  the nominal robot, stable on all of them.

  A model is `{:mass :friction :delay}`, the delay in steps of `dt`."
  (:require [allgo.numerics.optimize :as opt]
            [allgo.search :as search]))

(def ^:private dt 0.05)
(def ^:private force-limit 20.0)
(def ^:private steps 200)

(defn- rollout-step
  "One step of the policy `[k1 k2]` on `model`: the force it asks for
  queued, the one asked `delay` steps ago applied, and the cost of the
  position and the effort added. nil after the last step, or once the
  mass has run off."
  [{:keys [mass friction]} [k1 k2] {:keys [t x v queue cost] :as st}]
  (when (and (< t steps) (<= (abs x) 20.0))
    (let [u (max (- force-limit) (min force-limit (- (- (* k1 x)) (* k2 v))))
          queue (conj queue u)
          a (/ (- (first queue) (* friction v)) mass)
          v' (+ v (* dt a))]
      (assoc st :t (inc t) :x (+ x (* dt v')) :v v' :queue (subvec queue 1)
             :cost (+ cost (* dt (+ (* x x) (* 0.001 u u))))))))

(defn rollout
  "The states of the policy `gains` on `model` from rest at x = 1: each
  `{:t :x :v :cost}`."
  [{:keys [delay] :as model} gains]
  (search/step-seq (partial rollout-step model gains) {:t 0 :x 1.0 :v 0.0 :queue (vec (repeat delay 0.0)) :cost 0.0}))

(defn cost
  "The cost of the policy `gains` on `model`: the integral of x^2 and a
  thousandth of u^2, and 100 more if the mass ran off."
  [model gains]
  (let [{:keys [x cost]} (last (rollout model gains))]
    (+ cost (if (> (abs x) 20.0) 100.0 0.0))))

(defn search-policy
  "The gains minimizing the mean cost over `models`, by Nelder-Mead from
  [1 1]."
  [models]
  (:x (opt/nelder-mead (fn [gains] (/ (reduce + (map #(cost % gains) models)) (count models)))
                       [1.0 1.0] {:max-iter 800})))

(defn randomized
  "`n` models drawn at random: mass uniform on `[lo hi]`, friction on [0
  `friction`], the delay as given."
  [rng n [lo hi] friction delay]
  (vec (repeatedly n (fn [] {:mass (+ lo (* (- hi lo) (rng))) :friction (* friction (rng)) :delay delay}))))
