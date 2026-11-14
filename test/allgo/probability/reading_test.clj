(ns allgo.probability.reading-test
  "Reading degraded text with a generative program, after Russell and
  Norvig, *Artificial Intelligence: A Modern Approach*, 4th edition,
  section 15.4: figure 15.13's word, read by MCMC in 25 iterations, and
  figure 15.14's lesson, that a bigram model of English reads noisy text
  better than independent letters."
  (:require [allgo.probability.reading :as reading]
            [allgo.random :as random]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- edit-distance [a b]
  (let [n (count b)]
    (peek (reduce (fn [prev ca]
                    (reduce (fn [row j]
                              (conj row (min (inc (peek row)) (inc (prev (inc j)))
                                             (+ (prev j) (if (= ca (nth b j)) 0 1)))))
                            [(inc (first prev))] (range n)))
                  (vec (range (inc n))) a))))

(deftest the-font
  (testing "every two glyphs differ in at least two pixels"
    (is (every? (fn [[a b]] (>= (reduce + (map (fn [ra rb] (reduce + (map #(abs (- %1 %2)) ra rb)))
                                               (reading/glyphs a) (reading/glyphs b)))
                                2.0))
                (for [a reading/letters b reading/letters :when (pos? (compare a b))] [a b])))))

(deftest figure-15-13
  (testing "\"uncertainty\" read from a noisy image, by three runs of 25 MCMC iterations, each from an empty trace"
    (let [m (reading/model)
          image (reading/degrade (random/rng 1) (reading/render "uncertainty" 14 false) 0.1)]
      (doseq [seed [2 3 4]]
        (is (= "uncertainty" (:text (nth (reading/mcmc-seq m image (random/rng seed)) 25))))))))

(deftest exact-against-every-text
  (testing "the exact posterior over lengths, cell by cell and by the forward algorithm, is the sum over every text and noise of the joint probability"
    (let [image (reading/degrade (random/rng 5) (reading/render "an" 2 false) 0.6)
          texts (concat [""] (map str reading/letters) (for [a reading/letters b reading/letters] (str a b)))
          noises (mapv #(+ 0.1 (* 0.05 %)) (range 19))]
      (doseq [m [(reading/model :max-letters 2)
                 (reading/model :max-letters 2 :bigrams (reading/english-bigrams reading/common-words))]]
        (let [lj (for [t texts noise noises] [(count t) (reading/log-joint m image t noise)])
              top (reduce max (map second lj))
              z (reduce + (map #(math/exp (- (second %) top)) lj))
              brute (reduce (fn [acc [n w]] (update acc n (fnil + 0.0) (/ (math/exp (- w top)) z))) {} lj)
              exact (:length (reading/exact-posterior m image))]
          (doseq [n [0 1 2]]
            (is (< (abs (- (brute n) (exact n))) 1e-9))))))))

(deftest mcmc-against-exact
  (testing "MCMC samples the exact posterior: the length, and each letter"
    (let [m (reading/model :max-letters 6)
          image (reading/degrade (random/rng 6) (reading/render "cat" 6 false) 0.7)
          exact (reading/exact-posterior m image)
          samples (take 400 (drop 20 (reading/mcmc-seq m image (random/rng 7))))
          length (frequencies (map (comp count :text) samples))]
      (is (< (abs (- (get-in exact [:length 3]) (/ (get length 3 0) 400.0))) 0.05))
      (doseq [[i ch] (map-indexed vector "cat")]
        (let [of-length-3 (filter #(= 3 (count (:text %))) samples)
              p (/ (count (filter #(= ch (nth (:text %) i)) of-length-3)) (double (count of-length-3)))]
          (is (< (abs (- p (get-in exact [:letters i ch]))) 0.05) (str ch)))))))

(deftest figure-15-14
  (testing "from noisy images of English words it was not estimated on, a bigram model of English reads far closer to the truth than independent letters"
    (let [held-out ["reading" "letters" "noisy" "program" "monitor" "station" "citation" "pattern"
                    "planet" "garden" "silver" "window" "corner" "button" "winter"]
          bigrams (reading/english-bigrams reading/common-words)
          errors (fn [m] (reduce + (map-indexed (fn [k w]
                                                  (let [image (reading/degrade (random/rng (+ 100 k)) (reading/render w 14 false) 0.8)]
                                                    (edit-distance w (:text (reading/exact-posterior m image)))))
                                                held-out)))]
      (is (< (errors (reading/model :bigrams bigrams)) (* 0.75 (errors (reading/model))))))))
