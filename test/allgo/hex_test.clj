(ns allgo.hex-test
  (:require [allgo.geometry.hex :as h]
            [allgo.search :as search]
            [clojure.test :refer [deftest is testing]]))

(def ^:private all-kinds [:odd-r :even-r :odd-q :even-q])

(defn- patch
  "A block of hexes covering both signs on both axes."
  ([] (patch 8))
  ([n] (for [q (range (- n) (inc n)) r (range (- n) (inc n))] [q r])))

;; ---------------------------------------------------------------------------

(deftest coordinate-test
  (testing "cube adds the derived third coordinate"
    (is (= [1 2 -3] (h/cube [1 2])))
    (is (= [1 2 -3] (h/cube [1 2 -3]))) "idempotent")

  (testing "axial drops it"
    (is (= [1 2] (h/axial [1 2 -3])))
    (is (= [1 2] (h/axial [1 2]))))

  (testing "every cube coordinate sums to zero"
    (is (every? h/cube? (map h/cube (patch 5))))
    (is (not (h/cube? [1 1 1])))))

(deftest arithmetic-test
  (is (= [3 1] (h/add [1 2] [2 -1])))
  (is (= [-1 3] (h/subtract [1 2] [2 -1])))
  (is (= [3 6] (h/scale [1 2] 3)))

  (testing "length counts steps from the origin"
    (is (= 0 (h/hex-length [0 0])))
    (is (every? #(= 1 (h/hex-length %)) h/directions))
    (is (= 3 (h/hex-length [3 0])))
    (is (= 3 (h/hex-length [0 -3])))
    (is (= 3 (h/hex-length [3 -3]))))

  (testing "distance is a metric"
    (let [hexes (vec (patch 3))]
      (doseq [a (take 20 hexes) b (take 20 (drop 7 hexes))]
        (is (= (h/distance a b) (h/distance b a)) "symmetric")
        (is (= (zero? (h/distance a b)) (= a b)) "zero only on itself"))
      (doseq [a (take 8 hexes) b (take 8 (drop 5 hexes)) c (take 8 (drop 11 hexes))]
        (is (<= (h/distance a c) (+ (h/distance a b) (h/distance b c)))
            "triangle inequality")))))

(deftest direction-test
  (testing "the six neighbours are each one step away and all distinct"
    (let [ns (h/neighbors [3 -1])]
      (is (= 6 (count (distinct ns))))
      (is (every? #(= 1 (h/distance [3 -1] %)) ns))))

  (testing "diagonals are two steps away and off-axis"
    ;; Not the directions doubled: [2 0] is two steps along one axis,
    ;; whereas the diagonal [2 -1] shares only a corner with the centre.
    (is (every? #(= 2 (h/distance [0 0] %)) h/diagonals))
    (is (every? (fn [d] (not-any? zero? (h/cube d))) h/diagonals))
    (is (= 6 (count (distinct (h/diagonal-neighbors [2 2]))))))

  (testing "direction indices wrap"
    (is (= (h/direction 0) (h/direction 6) (h/direction -6)))
    (is (= (h/neighbor [0 0] 1) (h/neighbor [0 0] 7))))

  (testing "walking a direction moves one step per element"
    (is (= [[0 0] [1 0] [2 0] [3 0]] (take 4 (h/walk [0 0] 0))))
    (is (every? (fn [[i hex]] (= i (h/distance [0 0] hex)))
                (map-indexed vector (take 7 (h/walk [0 0] 3)))))))

(deftest rotation-test
  (testing "six turns is the identity, and the two directions undo each other"
    (doseq [hex [[2 -1] [0 0] [3 1] [-2 5]]]
      (is (= hex (nth (iterate h/rotate-left hex) 6)))
      (is (= hex (nth (iterate h/rotate-right hex) 6)))
      (is (= hex (h/rotate-left (h/rotate-right hex))))))

  (testing "rotation preserves distance from the point it turns about"
    (doseq [hex (patch 3)]
      (is (= (h/hex-length hex) (h/hex-length (h/rotate-left hex))))
      (is (= (h/distance [2 2] hex) (h/distance [2 2] (h/rotate-left hex [2 2]))))))

  (testing "rotating about a centre fixes that centre"
    (is (= [2 2] (h/rotate-left [2 2] [2 2])))
    (is (= [2 2] (h/rotate-right [2 2] [2 2]))))

  (testing "a ring rotates onto itself"
    (is (= (set (h/ring [0 0] 3))
           (set (map h/rotate-left (h/ring [0 0] 3)))))))

(deftest reflection-test
  (testing "reflecting twice is the identity"
    (doseq [hex (patch 3)
            f   [h/reflect-q h/reflect-r h/reflect-s]]
      (is (= (h/axial hex) (f (f hex))))))

  (testing "each keeps its own axis fixed"
    (is (every? #(= (first (h/cube %)) (first (h/cube (h/reflect-q %)))) (patch 3)))
    (is (every? #(= (second (h/cube %)) (second (h/cube (h/reflect-r %)))) (patch 3)))
    (is (every? #(= (nth (h/cube %) 2) (nth (h/cube (h/reflect-s %)) 2)) (patch 3))))

  (testing "reflection preserves distance from the origin, and from a centre"
    (doseq [hex (patch 3)]
      (is (= (h/hex-length hex) (h/hex-length (h/reflect-q hex))))
      (is (= (h/distance [1 -2] hex) (h/distance [1 -2] (h/reflect-r hex [1 -2])))))))

(deftest round-test
  (testing "an exact hex rounds to itself"
    (is (every? #(= % (h/round %)) (patch 4))))

  (testing "the result is always a legal hex, however it is nudged"
    (doseq [[q r] (patch 3)
            dq [-0.4 -0.1 0.0 0.1 0.4]
            dr [-0.4 -0.1 0.0 0.1 0.4]]
      (is (h/cube? (h/cube (h/round [(+ q dq) (+ r dr)]))))))

  (testing "it picks a hex adjacent to where it started"
    (doseq [[q r] (patch 2) dq [-0.3 0.3] dr [-0.3 0.3]]
      (is (<= (h/distance [q r] (h/round [(+ q dq) (+ r dr)])) 1)))))

(deftest line-test
  (testing "a line runs from a to b, one step at a time"
    (doseq [target (h/hexes-within [0 0] 4)]
      (let [l (h/line [0 0] target)]
        (is (= [0 0] (first l)))
        (is (= target (last l)))
        (is (= (inc (h/distance [0 0] target)) (count l)) "no gaps or repeats")
        (is (every? (fn [[a b]] (= 1 (h/distance a b))) (partition 2 1 l))
            "each step is to a neighbour"))))

  (testing "a line to itself is just itself"
    (is (= [[2 2]] (h/line [2 2] [2 2]))))

  (testing "the line is reversible"
    (is (= (h/line [0 0] [3 -1]) (reverse (h/line [3 -1] [0 0]))))))

(deftest region-test
  (testing "hexes-within has 1 + 3n(n+1) hexes, all in range"
    (doseq [n (range 6)]
      (let [hs (h/hexes-within [2 -1] n)]
        (is (= (+ 1 (* 3 n (inc n))) (count hs)))
        (is (= (count hs) (count (distinct hs))))
        (is (every? #(<= (h/distance [2 -1] %) n) hs)))))

  (testing "a ring has 6n hexes, every one exactly n away"
    (is (= [[0 0]] (h/ring [0 0] 0)))
    (doseq [n (range 1 6)]
      (let [r (h/ring [1 1] n)]
        (is (= (* 6 n) (count r)))
        (is (= (count r) (count (distinct r))))
        (is (every? #(= n (h/distance [1 1] %)) r))
        (is (every? (fn [[a b]] (= 1 (h/distance a b)))
                    (partition 2 1 (concat r [(first r)])))
            "the ring is a closed walk"))))

  (testing "a spiral is the rings stacked, and covers the same ground"
    (doseq [n (range 5)]
      (is (= (set (h/hexes-within [0 0] n)) (set (h/spiral [0 0] n))))
      (is (= (count (h/hexes-within [0 0] n)) (count (h/spiral [0 0] n))))
      (is (= [0 0] (first (h/spiral [0 0] n))))
      (is (apply <= (map #(h/distance [0 0] %) (h/spiral [0 0] n)))
          "ordered outward")))

  (testing "intersecting ranges is exactly the set intersection"
    (doseq [[c1 n1 c2 n2] [[[0 0] 2 [2 0] 2]
                           [[0 0] 3 [0 0] 1]
                           [[0 0] 1 [9 9] 1]
                           [[-2 3] 4 [3 -1] 3]]]
      (is (= (into #{} (filter (set (h/hexes-within c2 n2))) (h/hexes-within c1 n1))
             (set (h/intersecting-ranges c1 n1 c2 n2)))))))

(deftest map-shape-test
  (is (= 12 (count (h/parallelogram 0 2 0 3))))
  (is (every? #(h/cube? (h/cube %)) (h/parallelogram -2 2 -2 2)))

  (testing "a triangle of side n has n(n+1)/2 hexes"
    (doseq [n (range 1 6)]
      (is (= (/ (* (inc n) (+ n 2)) 2) (count (h/triangle n))))))

  (is (= (set (h/hexes-within [0 0] 3)) (set (h/hexagon 3))))

  (testing "a rectangle holds width*height hexes with no repeats"
    (doseq [o [:pointy :flat]]
      (let [rect (h/rectangle o 5 4)]
        (is (= 20 (count rect)))
        (is (= 20 (count (distinct rect))))))))

(deftest offset-test
  (testing "every offset system round-trips, negatives included"
    (doseq [k all-kinds]
      (is (every? #(= % (h/offset-> k (h/->offset k %))) (patch))
          (str k))))

  (testing "and is a bijection -- no two hexes share a cell"
    (doseq [k all-kinds]
      (let [hexes (vec (patch))]
        (is (= (count hexes) (count (set (map #(h/->offset k %) hexes)))) (str k)))))

  (testing "the shifted parity is the one the name says"
    ;; odd-r pushes odd rows right; even-r pushes even rows right.
    (is (= [0 1] (h/->offset :odd-r [0 1])))
    (is (= [1 1] (h/->offset :even-r [0 1])))
    (is (= [0 0] (h/->offset :odd-r [0 0])))
    (is (= [0 0] (h/->offset :even-r [0 0]))))

  (testing "neighbours stay adjacent through the conversion"
    (doseq [k all-kinds]
      (is (every? #(= 1 (h/distance [2 -3] (h/offset-> k (h/->offset k %))))
                  (h/neighbors [2 -3]))
          (str k)))))

(deftest doubled-test
  (testing "both doubled systems round-trip"
    (is (every? #(= % (h/doublewidth-> (h/->doublewidth %))) (patch)))
    (is (every? #(= % (h/doubleheight-> (h/->doubleheight %))) (patch))))

  (testing "their own distance formulas agree with hex distance"
    ;; This is what the doubled systems are for: you can measure without
    ;; converting back. If the conversion were wrong, this would diverge.
    (let [hexes (h/hexagon 3)]
      (doseq [a hexes b hexes]
        (is (= (h/distance a b)
               (h/doublewidth-distance (h/->doublewidth a) (h/->doublewidth b))))
        (is (= (h/distance a b)
               (h/doubleheight-distance (h/->doubleheight a) (h/->doubleheight b)))))))

  (testing "one axis really is doubled"
    (is (= [2 0] (h/->doublewidth [1 0])))
    (is (= [0 2] (h/->doubleheight [0 1])))))

(deftest layout-test
  (doseq [o [:pointy :flat]]
    (let [l (h/layout o 12.0 [50.0 40.0])]
      (testing (str o ": the origin hex sits at the layout origin")
        (is (= [50.0 40.0] (h/->pixel l [0 0]))))

      (testing (str o ": pixel conversion round-trips")
        (is (every? #(= % (h/pixel->hex l (h/->pixel l %))) (h/hexagon 6))))

      (testing (str o ": a point near a centre lands in that hex")
        (doseq [hex (h/hexagon 3)]
          (let [[x y] (h/->pixel l hex)]
            (is (= hex (h/pixel->hex l [(+ x 2.0) (- y 1.5)]))))))

      (testing (str o ": corners sit one radius out, six of them")
        (let [[cx cy] (h/->pixel l [0 0])
              cs      (h/corners l [0 0])]
          (is (= 6 (count cs)))
          (is (every? (fn [[x y]]
                        (< (abs (- 12.0 (Math/hypot (- x cx) (- y cy)))) 1e-9))
                      cs))))

      (testing (str o ": neighbouring centres are all the same distance apart")
        (let [[cx cy] (h/->pixel l [0 0])
              ds      (map (fn [n] (let [[x y] (h/->pixel l n)]
                                     (Math/hypot (- x cx) (- y cy))))
                           (h/neighbors [0 0]))]
          (is (< (- (apply max ds) (apply min ds)) 1e-9))))

      (testing (str o ": size scales the layout linearly")
        (let [small (h/layout o 1.0)
              big   (h/layout o 10.0)]
          (is (every? (fn [hex]
                        (let [[sx sy] (h/->pixel small hex)
                              [bx by] (h/->pixel big hex)]
                          (and (< (abs (- (* 10 sx) bx)) 1e-9)
                               (< (abs (- (* 10 sy) by)) 1e-9))))
                      (h/hexagon 3)))))))

  (testing "the two orientations differ"
    (is (not= (h/->pixel (h/layout :pointy 10.0) [1 1])
              (h/->pixel (h/layout :flat 10.0) [1 1])))))

(deftest obstacle-test
  (testing "with nothing in the way, reachable is just the range"
    (doseq [n (range 4)]
      (is (= (set (h/hexes-within [0 0] n)) (h/reachable [0 0] n (constantly false))))))

  (testing "a wall all the way around leaves nowhere to go"
    (is (= #{[0 0]} (h/reachable [0 0] 5 (set (h/ring [0 0] 1))))))

  (testing "reachable flows around a wall rather than through it"
    (let [wall (set (for [r (range -4 4)] [1 r]))
          got  (h/reachable [0 0] 3 wall)]
      (is (not-any? wall got))
      (is (contains? got [0 0]))
      (is (< (count got) (count (h/hexes-within [0 0] 3)))
          "the wall costs it some hexes")))

  (testing "visibility stops at a wall but includes it"
    (let [wall #{[2 0]}
          seen (h/visible [0 0] 4 wall)]
      (is (contains? seen [2 0]) "you can see the wall")
      (is (not (contains? seen [4 0])) "but not straight through it")
      (is (contains? seen [0 2]) "other directions are unaffected")))

  (testing "with nothing blocking, everything in range is visible"
    (is (= (set (h/hexes-within [0 0] 3)) (h/visible [0 0] 3 (constantly false))))))

(deftest pathfinding-test
  (let [path-between (fn [start goal blocked]
                       (search/a-star {:start      start
                                       :goal       goal
                                       :neighbours #(remove blocked (h/neighbors %))
                                       ;; hex/distance is exactly the right
                                       ;; heuristic and has the right shape
                                       ;; to be handed over as-is.
                                       :heuristic  h/distance}))]
    (testing "in the open, the path is as short as the distance allows"
      (doseq [goal (h/ring [0 0] 4)]
        (let [p (path-between [0 0] goal #{})]
          (is (= (inc (h/distance [0 0] goal)) (count p)))
          (is (= [0 0] (first p)))
          (is (= goal (last p))))))

    (testing "every step of a path is to a neighbour"
      (let [wall (set (for [r (range -3 4) :when (not= r 3)] [2 r]))
            p    (path-between [0 0] [4 0] wall)]
        (is (every? (fn [[a b]] (= 1 (h/distance a b))) (partition 2 1 p)))
        (is (not-any? wall p) "and never enters a wall")
        (is (> (count p) (inc (h/distance [0 0] [4 0]))) "the detour costs it")))

    (testing "an unreachable goal gives nil, not a wrong path"
      (is (nil? (path-between [0 0] [5 0] (set (h/ring [0 0] 1))))))

    (testing "the heuristic does not change the answer, only the work"
      (doseq [goal [[4 -2] [-3 1] [0 5]]]
        (is (= (count (path-between [0 0] goal #{}))
               (count (search/dijkstra {:start [0 0] :goal goal
                                        :neighbours h/neighbors}))))))))

(deftest wraparound-test
  (testing "a hex already inside is left alone"
    (is (every? #(= % (h/wrap % 3)) (h/hexagon 3))))

  (testing "anything outside comes back inside"
    (doseq [hex (h/hexes-within [0 0] 7)]
      (is (<= (h/distance [0 0] (h/wrap hex 3)) 3))))

  (testing "there are seven centres and they are distinct"
    (let [cs (h/mirror-centers 3)]
      (is (= 7 (count cs)))
      (is (= 7 (count (distinct cs))))
      (is (some #{[0 0]} cs))))

  (testing "wrapping is idempotent"
    (is (every? #(= (h/wrap % 3) (h/wrap (h/wrap % 3) 3)) (h/hexes-within [0 0] 9))))

  (testing "hexes a whole map apart wrap to the same place"
    (doseq [center (rest (h/mirror-centers 3))
            hex    (h/hexagon 3)]
      (is (= hex (h/wrap (h/add hex center) 3)))))

  (testing "a long walk stays on the map"
    ;; Several map widths out. One subtraction only pulls a hex back from
    ;; an adjacent copy, so this is what catches a wrap that does not repeat.
    (let [walk (map #(h/wrap % 3) (take 60 (h/walk [0 0] 0)))]
      (is (every? some? walk) "nothing falls off the edge")
      (is (every? #(<= (h/hex-length %) 3) walk)))))

(deftest depth-first-needs-a-bounded-space-test
  ;; Depth-first follows one branch to its end before trying another, so on
  ;; the open hex grid it walks away from the goal indefinitely -- measured,
  ;; it exhausts the heap. Bounded, it behaves.
  (let [board (set (h/hexagon 4))
        p     (search/depth-first {:start      [0 0]
                                   :goal       [3 -1]
                                   :neighbours #(filterv board (h/neighbors %))})]
    (is (some? p) "it finds a route")
    (is (= [0 0] (first p)))
    (is (= [3 -1] (last p)))
    (is (every? (fn [[a b]] (= 1 (h/distance a b))) (partition 2 1 p))
        "and it is a real walk, even if a roundabout one")
    (is (every? board p))
    (is (>= (count p) (inc (h/distance [0 0] [3 -1])))
        "no shorter than the shortest path, and usually longer")))
