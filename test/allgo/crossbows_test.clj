(ns allgo.crossbows-test
  (:require [allgo.geometry.vec3 :as v]
            [allgo.physics.rigid :as rigid]
            [allgo.simulation.crossbows :as c]
            [clojure.test :refer [deftest is testing]]))

(defn- play-out
  "Step `game` until the shot in flight has been resolved."
  [game]
  (loop [g game] (if (= :flying (:phase g)) (recur (c/advance g)) g)))

(defn- tower-of [game side]
  (some #(when (and (= :tower (:kind %)) (= side (:side %))) %) (:bodies (:world game))))

(deftest board-test
  (testing "each side has the set's pieces"
    (doseq [side c/sides
            kind [:tower :wall :flag :warrior]]
      (is (= (c/pieces kind)
             (count (filter #(and (= side (:side %)) (= kind (:kind %))) (c/board))))
          (str (name side) " " (name kind)))))

  (testing "a new game stands, asleep, with nothing down"
    (let [g (c/new-game)]
      (is (= {:vikings {} :barbarians {}} (c/tally (:world g))))
      (is (every? rigid/inert? (:bodies (:world g))))
      (is (= :vikings (:turn g)))
      (is (= {:vikings 10 :barbarians 10} (:discs g))))))

(deftest aim-test
  (testing "the catapult's solved shot passes through its target"
    ;; In a vacuum, which is what the solution assumes: the parabola from
    ;; the launch point, sampled finely, comes within a few centimeters.
    (doseq [side c/sides]
      (let [g (c/new-game)
            target (c/tower-face (:world g) side)
            {:keys [weapon yaw pitch speed]} (c/aim-at side :catapult target)
            p0 (c/launch-point side weapon)
            v0 (c/launch-velocity side weapon yaw pitch speed)
            miss (apply min (for [i (range 4000)
                                  :let [t (* i 0.001)]]
                              (v/distance target (v/add (v/add p0 (v/scale v0 t))
                                                        (v/scale c/gravity (* 0.5 t t))))))]
        (is (< miss 0.05) (str (name side) " missed by " miss))))))

(deftest shot-test
  (testing "a clean catapult hit high on the tower's face knocks it over"
    (let [g (c/new-game)
          {:keys [weapon yaw pitch speed]} (c/aim-at :vikings :catapult (c/tower-face (:world g) :vikings))
          after (play-out (c/fire g weapon yaw pitch speed))]
      (is (c/down? (tower-of after :barbarians)))
      (is (= :vikings (:winner after)))
      (is (= :over (:phase after)))))

  (testing "a crossbow shot slides through the gate to the tower's foot"
    ;; The gate is the only way in along the ground. The shot reaching the
    ;; tower at all is the test: it moves it, and a wall in the way would
    ;; have stopped the disc short.
    (let [g (c/new-game)
          [x _ z :as before] (:pos (tower-of g :barbarians))
          {:keys [weapon yaw pitch speed]} (c/aim-at :vikings :crossbow [x 0.1 z])
          after (play-out (c/fire g weapon yaw pitch speed))]
      (is (> (v/distance before (:pos (tower-of after :barbarians))) 0.01))
      (is (not (c/down? (tower-of after :barbarians))) "a crossbow shot is not meant to be enough")
      (is (nil? (:winner after))))))

(deftest turn-test
  (testing "a shot uses a disc and passes the turn"
    (let [g (c/new-game)
          after (play-out (c/fire g :catapult 30.0 35.0 8.0))]
      (is (= 9 (get-in after [:discs :vikings])))
      (is (= :barbarians (:turn after)))
      (is (= :aiming (:phase after)))))

  (testing "nothing happens out of turn, or out of discs"
    (let [g (c/new-game)]
      (let [busy (assoc g :phase :flying)]
        (is (= busy (c/fire busy :catapult 0.0 35.0 10.0))))
      (is (= :aiming (:phase (c/fire (assoc-in g [:discs :vikings] 0) :catapult 0.0 35.0 10.0))))))

  (testing "both out of discs with both towers up is a draw"
    (let [g (-> (c/new-game) (assoc :discs {:vikings 1 :barbarians 0}))
          after (play-out (c/fire g :catapult 30.0 35.0 8.0))]
      (is (= :draw (:winner after))))))
