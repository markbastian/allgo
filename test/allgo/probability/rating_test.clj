(ns allgo.probability.rating-test
  "Bayesian skill ratings, after Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, section 15.1.2. The book
  gives the model and no numbers, so the update is checked against the
  exact posterior, by integration, and the ratings against skills that
  generated the games."
  (:require [allgo.numerics.special :as special]
            [allgo.probability.rating :as rating]
            [allgo.random :as random]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private beta 4.0)

(defn- exact-one-game
  "The posterior over two skills with Gaussian priors after the first
  beats the second, by integration on a grid: the evidence (the win
  probability), and each skill's posterior mean and standard deviation.
  The likelihood of the win given the skills is Phi((s1 - s2) / (sqrt 2
  beta)), the two performances' difference being N(s1 - s2, 2 beta^2)."
  [{m1 :mu s1 :sigma} {m2 :mu s2 :sigma}]
  (let [n 500
        grid (fn [m s] (mapv #(+ (- m (* 8.0 s)) (* % (/ (* 16.0 s) n))) (range (inc n))))
        g1 (grid m1 s1) g2 (grid m2 s2)
        h1 (/ (* 16.0 s1) n) h2 (/ (* 16.0 s2) n)
        pdf (fn [x m s] (/ (special/normal-pdf (/ (- x m) s)) s))
        cells (for [x g1 y g2]
                [x y (* (pdf x m1 s1) (pdf y m2 s2) (special/normal-cdf (/ (- x y) (* (math/sqrt 2.0) beta))) h1 h2)])
        z (reduce + (map #(nth % 2) cells))
        moment (fn [f] (/ (reduce + (map (fn [[x y w]] (* w (f x y))) cells)) z))
        e1 (moment (fn [x _] x)) e2 (moment (fn [_ y] y))]
    {:evidence z
     :mu1 e1 :sigma1 (math/sqrt (- (moment (fn [x _] (* x x))) (* e1 e1)))
     :mu2 e2 :sigma2 (math/sqrt (- (moment (fn [_ y] (* y y))) (* e2 e2)))}))

(deftest one-game-is-exact
  (testing "for a single game between Gaussian ratings, moment matching gives the exact posterior's means and spreads, and Phi(.../c) the evidence"
    (doseq [[a b] [[(rating/rating 25 8) (rating/rating 25 8)]
                   [(rating/rating 30 3) (rating/rating 20 6)]
                   [(rating/rating 18 5) (rating/rating 27 2)]]]
      (let [exact (exact-one-game a b)
            after (rating/update-ratings {:a a :b b} {:winner [:a] :loser [:b]} {:beta beta})]
        (is (< (abs (- (:evidence exact) (rating/win-probability {:a a :b b} [:a] [:b] beta))) 1e-6))
        (is (< (abs (- (:mu1 exact) (get-in after [:a :mu]))) 1e-4))
        (is (< (abs (- (:sigma1 exact) (get-in after [:a :sigma]))) 1e-4))
        (is (< (abs (- (:mu2 exact) (get-in after [:b :mu]))) 1e-4))
        (is (< (abs (- (:sigma2 exact) (get-in after [:b :sigma]))) 1e-4))))))

(deftest surprises-teach-more
  (let [prior {:strong (rating/rating 30 4) :weak (rating/rating 20 4)}
        moved (fn [game] (let [r (rating/update-ratings prior game {:beta beta})]
                           (- (get-in r [:strong :mu]) 30.0)))]
    (testing "an upset moves the ratings further than the expected result does"
      (is (> (abs (moved {:winner [:weak] :loser [:strong]}))
             (* 3 (abs (moved {:winner [:strong] :loser [:weak]}))))))
    (testing "and every game, won or lost, makes a rating more certain"
      (is (every? #(< (:sigma %) 4.0) (vals (rating/update-ratings prior {:winner [:weak] :loser [:strong]} {:beta beta})))))))

(deftest team-games
  (testing "skills recovered from team games alone, the teams changing from game to game, as the section says they must"
    (let [rng (random/rng 15)
          skills {:p1 10.0 :p2 14.0 :p3 18.0 :p4 22.0 :p5 26.0 :p6 30.0 :p7 34.0 :p8 38.0}
          players (vec (keys skills))
          perform (fn [team] (reduce + (map #(random/gaussian rng (skills %) beta) team)))
          games (vec (for [_ (range 3000)
                           :let [[a b c d] (take 4 (random/shuffle rng players))
                                 [t1 t2] [[a b] [c d]]]]
                       (if (> (perform t1) (perform t2)) {:winner t1 :loser t2} {:winner t2 :loser t1})))
          prior (zipmap players (repeat (rating/rating 25 8)))
          final (last (rating/rating-seq prior games {:beta beta :tau 0.05}))
          by-mean (sort-by #(:mu (final %)) players)]
      (is (= (sort-by skills players) by-mean) "ranked in the order of their true skills")
      (is (= (reverse (sort-by skills players)) (rating/ranking final)))
      ;; only differences of skill are seen, so the means are pinned only
      ;; relative to each other: compare them centered
      (let [center (fn [m] (let [avg (/ (reduce + (vals m)) (count m))] (update-vals m #(- % avg))))
            truth (center skills)
            estimate (center (update-vals final :mu))]
        (is (every? #(< (abs (- (truth %) (estimate %))) 2.5) players) "each within 2.5 of the truth")))))
