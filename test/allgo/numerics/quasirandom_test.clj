(ns allgo.numerics.quasirandom-test
  (:require [allgo.numerics.quasirandom :as qr]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

;; Low-dispersion and low-discrepancy sampling, after LaValle, *Planning
;; Algorithms*, Cambridge University Press, 2006, sections 5.2.3-5.2.4.

(deftest sequences-test
  (testing "van der Corput in base 2: 1/2, 1/4, 3/4, 1/8, 5/8, 3/8, 7/8 -- each in the largest gap left"
    (is (= [0.5 0.25 0.75 0.125 0.625 0.375 0.875] (mapv #(qr/van-der-corput % 2) (range 1 8)))))
  (testing "base 3: 1/3, 2/3, 1/9, 4/9, 7/9"
    (is (every? #(< (abs %) 1e-15) (map - [(/ 1.0 3) (/ 2.0 3) (/ 1.0 9) (/ 4.0 9) (/ 7.0 9)] (mapv #(qr/van-der-corput % 3) (range 1 6))))))
  (testing "a Sukharev grid of 4 a side in two dimensions: the sixteen cell centers"
    (is (= 16 (count (qr/sukharev 4 2))))
    (is (= [0.125 0.125] (first (qr/sukharev 4 2))))))

(deftest coverage-test
  (let [n 256
        probes (qr/sukharev 40 2)
        corners (qr/sukharev 20 2)
        random-points (qr/uniform n 2 (random/rng 1))]
    (testing "the Sukharev grid has the least dispersion: no point farther than half a cell, 1/32, from it"
      (is (<= (qr/dispersion (qr/sukharev 16 2) probes) (+ (/ 1.0 32) 1e-12))))
    (testing "Halton and Hammersley leave smaller gaps than as many random points"
      (is (< (qr/dispersion (qr/halton n 2) probes) (qr/dispersion random-points probes)))
      (is (< (qr/dispersion (qr/hammersley n 2) probes) (qr/dispersion random-points probes))))
    (testing "and fill boxes nearer their share: a lower star discrepancy than random points"
      (is (< (qr/star-discrepancy (qr/halton n 2) corners) (qr/star-discrepancy random-points corners)))
      (is (< (qr/star-discrepancy (qr/hammersley n 2) corners) (qr/star-discrepancy random-points corners))))))
