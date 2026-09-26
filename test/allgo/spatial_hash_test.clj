(ns allgo.spatial-hash-test
  (:require [allgo.spatial.hash :as sh]
            [clojure.test :refer [deftest is testing]]))

(defn- positions [coords] (double-array (map double coords)))

(defn- brute-force-within
  "The answer the hash has to match, computed the slow obvious way."
  [^doubles p n i max-dist]
  (let [a  (* 3 i)
        d2 (* max-dist max-dist)]
    (set (for [j (range n)
               :when (not= j i)
               :let  [b  (* 3 j)
                      dx (- (aget p a) (aget p b))
                      dy (- (aget p (+ a 1)) (aget p (+ b 1)))
                      dz (- (aget p (+ a 2)) (aget p (+ b 2)))]
               :when (<= (+ (* dx dx) (* dy dy) (* dz dz)) d2)]
           j))))

(defn- brute-force-pairs
  "Direct double loop. Building it out of `brute-force-within` would
  recompute a set per object and turn a slow reference into an unusably
  slow one."
  [^doubles p n max-dist]
  (let [d2 (* max-dist max-dist)]
    (set (for [i (range n) j (range (inc i) n)
               :let  [a  (* 3 i) b (* 3 j)
                      dx (- (aget p a) (aget p b))
                      dy (- (aget p (+ a 1)) (aget p (+ b 1)))
                      dz (- (aget p (+ a 2)) (aget p (+ b 2)))]
               :when (<= (+ (* dx dx) (* dy dy) (* dz dz)) d2)]
           [i j]))))

(defn- random-positions [n seed spread]
  (let [r (java.util.Random. seed)]
    (positions (repeatedly (* 3 n) #(* spread (- (.nextDouble r) 0.5))))))

(defn- lattice
  "A cube of points one unit apart; every neighbor count is known."
  [side]
  (positions (for [x (range side) y (range side) z (range side) c [x y z]] c)))

;; ---------------------------------------------------------------------------

(deftest lattice-test
  (let [side 6
        n    (* side side side)
        p    (lattice side)
        h    (sh/spatial-hash 1.0 n)]
    (sh/rebuild! h p n)

    (testing "a corner point has three neighbors, an interior one six"
      (is (= 3 (count (sh/within h p 0 1.01))))
      ;; index = x*36 + y*6 + z, so 43 is (1,1,1)
      (is (= 6 (count (sh/within h p 43 1.01)))))

    (testing "widening the radius picks up the diagonals"
      ;; The twelve edge-diagonal neighbors sit at sqrt(2).
      (is (= 18 (count (sh/within h p 43 1.5)))))

    (testing "a radius covering everything returns everything but self"
      (is (= (dec n) (count (sh/within h p 0 100.0)))))

    (testing "nothing is reported twice"
      (doseq [i [0 43 100 (dec n)]]
        (let [found (sh/within h p i 2.0)]
          (is (= (count found) (count (distinct found))) (str "object " i)))))))

(deftest agrees-with-brute-force-test
  (testing "random clouds, at radii either side of the cell size"
    (doseq [n        [1 2 50 400]
            spread   [1.0 10.0]
            max-dist [0.05 0.3 1.0]]
      (let [p (random-positions n 1234 spread)
            h (sh/spatial-hash 0.3 n)]
        (sh/rebuild! h p n)
        (doseq [i (range n)]
          (is (= (brute-force-within p n i max-dist)
                 (set (sh/within h p i max-dist)))
              (str "n=" n " spread=" spread " dist=" max-dist " i=" i))))))

  (testing "overlapping-pairs matches too, with no pair repeated"
    (doseq [n [2 60 300]]
      (let [p     (random-positions n 99 4.0)
            h     (sh/spatial-hash 0.25 n)
            pairs (sh/overlapping-pairs h p n 0.25)]
        (is (= (brute-force-pairs p n 0.25) (set pairs)))
        (is (= (count pairs) (count (distinct pairs))))
        (is (every? (fn [[i j]] (< i j)) pairs) "each pair is reported once, ordered")))))

(deftest awkward-input-test
  (testing "negative coordinates are handled, not folded onto positive ones"
    ;; The hash multiplies cell indices and masks the result; if that were
    ;; done carelessly, mirrored points would collide or go missing.
    (let [p (positions [-5.0 -5.0 -5.0,  5.0 5.0 5.0,  -5.1 -5.0 -5.0])
          h (sh/spatial-hash 0.5 3)]
      (sh/rebuild! h p 3)
      (is (= #{2} (set (sh/within h p 0 0.2))))
      (is (= #{} (set (sh/within h p 1 0.2))) "the far point sees nobody")))

  (testing "objects piled at the same point all find each other"
    (let [n 20
          p (positions (repeat (* 3 n) 1.0))
          h (sh/spatial-hash 0.5 n)]
      (sh/rebuild! h p n)
      (is (= (dec n) (count (sh/within h p 0 1e-9))))))

  (testing "a cell size far from the query radius still gives the right answer"
    ;; Only the cost changes: too small and a query sweeps many cells, too
    ;; large and each holds objects that are nowhere near.
    (doseq [spacing [0.01 0.1 1.0 25.0]]
      (let [p (random-positions 120 7 5.0)
            h (sh/spatial-hash spacing 120)]
        (sh/rebuild! h p 120)
        (is (= (brute-force-within p 120 0 0.5)
               (set (sh/within h p 0 0.5)))
            (str "spacing " spacing)))))

  (testing "an empty table answers nothing rather than failing"
    (let [h (sh/spatial-hash 1.0 10)]
      (sh/rebuild! h (positions []) 0)
      (is (= 0 (sh/query-point! h 0.0 0.0 0.0 1.0)))))

  (testing "rebuilding replaces the contents, it does not accumulate"
    (let [h (sh/spatial-hash 1.0 8)
          near (positions [0.0 0.0 0.0, 0.1 0.0 0.0])
          far  (positions [0.0 0.0 0.0, 50.0 0.0 0.0])]
      (sh/rebuild! h near 2)
      (is (= 1 (count (sh/within h near 0 0.5))))
      (sh/rebuild! h far 2)
      (is (= 0 (count (sh/within h far 0 0.5))) "the stale neighbor is gone"))))

(deftest query-mechanics-test
  (let [p (lattice 4)
        n 64
        h (sh/spatial-hash 1.0 n)]
    (sh/rebuild! h p n)

    (testing "query! reports a count and neighbor reads that many out"
      (let [found (sh/query! h p 0 1.01)]
        (is (pos? found))
        (is (= found (count (sh/neighbors h found))))
        (is (every? #(< -1 % n) (sh/neighbors h found)) "all are real ids")))

    (testing "candidates include the object itself; within excludes it"
      (let [found (sh/query! h p 5 1.01)]
        (is (contains? (set (sh/neighbors h found)) 5))
        (is (not (contains? (set (sh/within h p 5 1.01)) 5)))))

    (testing "candidates are a superset of the true neighbors"
      ;; The hash is allowed to over-report -- that is the bargain that
      ;; lets it be approximate and cheap -- but never to under-report.
      (doseq [i (range n)]
        (let [found      (sh/query! h p i 1.01)
              candidates (set (sh/neighbors h found))]
          (is (every? candidates (brute-force-within p n i 1.01))
              (str "object " i)))))))

(deftest scale-test
  (testing "it stays exact at a few thousand objects"
    (let [n 1200
          p (random-positions n 2024 6.0)
          h (sh/spatial-hash 0.2 n)
          pairs (sh/overlapping-pairs h p n 0.2)]
      (is (= (brute-force-pairs p n 0.2) (set pairs)))
      (is (pos? (count pairs)) "and the test is not vacuous"))))
