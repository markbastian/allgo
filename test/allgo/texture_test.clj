(ns allgo.texture-test
  (:require [allgo.procedural.texture :as tx]
            [clojure.test :refer [deftest is testing]]))

(defn- points [n]
  (for [i (range n)] [(* 0.1379 i) (- (* 0.2713 i) 3.0) (* 0.0317 i)]))

(defn- close? [a b] (every? #(< (abs (double %)) 1e-9) (map - a b)))

(deftest catalog-test
  (testing "every texture gives a color in [0, 1], a finite height and a shine in [0, 1]"
    (doseq [{:keys [name make]} tx/catalog]
      (let [tex (make {})]
        (doseq [[x y z] (points 60)]
          (let [[r g b h shine :as t] (tex x y z)]
            (is (= 5 (count t)) name)
            (is (every? #(<= 0.0 % 1.0) [r g b shine]) name)
            (is (Double/isFinite h) name))))))

  (testing "a seed is a texture: the same one twice, the same everywhere; another, different"
    (doseq [{:keys [name make]} tx/catalog]
      (let [a (make {:seed 3}) b (make {:seed 3}) c (make {:seed 4})
            ps (points 40)]
        (is (every? (fn [[x y z]] (= (a x y z) (b x y z))) ps) name)
        (is (some (fn [[x y z]] (not= (a x y z) (c x y z))) ps) name)))))

(deftest brick-test
  (let [tex (tx/brick {:width 1.0 :height 0.36 :mortar 0.05})
        [_ _ _ face] (tex 0.5 0.18 0.0)
        [_ _ _ bed] (tex 0.5 0.36 0.0)
        [_ _ _ head] (tex 1.0 0.18 0.0)
        [_ _ _ shifted] (tex 1.0 0.54 0.0)]
    (testing "the joints are sunk below the faces, along the courses and between bricks"
      (is (< bed (- face 0.3)))
      (is (< head (- face 0.3))))
    (testing "and each course is offset half a brick, so a head joint below is a face above"
      (is (> shifted (- face 0.2))))))

(deftest flagstone-test
  (testing "mortar runs between stones: some points are sunk, most are not"
    (let [tex (tx/flagstone {})
          hs (for [i (range 80) j (range 80)] (nth (tex (* 0.05 i) (* 0.05 j) 0.0) 3))
          low (count (filter #(< % 0.1) hs))]
      (is (< 0.03 (/ low (count hs)) 0.4)))))

(deftest color-map-test
  (testing "the ramp passes through its interior knots"
    (let [ramp (tx/color-map [[0 0 0] [0 0 0] [1 0.5 0] [1 1 1] [1 1 1]])]
      (is (close? [0.0 0.0 0.0] (ramp 0.0)))
      (is (close? [1.0 0.5 0.0] (ramp 0.5)))
      (is (close? [1.0 1.0 1.0] (ramp 1.0))))))

(deftest shade-test
  (let [tex (tx/brick {})
        flat (fn [x y] (tx/plane-pixel tex x y 0.01 {:bump 0.0}))]
    (testing "lit colors stay in [0, 1]"
      (doseq [{:keys [make]} tx/catalog
              :let [t (make {})]
              [x y] [[0.1 0.2] [0.97 0.36] [2.3 -1.1]]]
        (is (every? #(<= 0.0 % 1.0) (tx/plane-pixel t x y 0.01 {})))))
    (testing "with no bump, a matte surface is its color times the light it gets"
      ;; Clouds have no shine, so no highlight of another color.
      (let [sky (tx/clouds {})
            [r g b] (sky 0.5 0.18 0.0)
            [lr lg lb] (tx/plane-pixel sky 0.5 0.18 0.01 {:bump 0.0})]
        (is (< (abs (- (/ lr r) (/ lg g))) 1e-9))
        (is (< (abs (- (/ lr r) (/ lb b))) 1e-9))))
    (testing "the bumps light a joint's two sides differently, and flat they are lit alike"
      ;; The bed joint at y = 0.36: the brick below it rounds over
      ;; toward the light, the one above rounds away from it.
      (let [lum (fn [[r g b]] (+ r g b))
            up (tx/plane-pixel tex 0.5 0.33 0.005 {})
            down (tx/plane-pixel tex 0.5 0.39 0.005 {})]
        (is (> (lum up) (+ (lum down) 0.1)))
        (is (< (abs (- (lum (flat 0.5 0.33)) (lum (flat 0.5 0.39)))) 0.25)))))

  (testing "a sphere is drawn inside the unit disk and nowhere else"
    (let [tex (tx/marble {})]
      (is (nil? (tx/sphere-pixel tex 0.8 0.7 1.0 0.01 {})))
      (is (some? (tx/sphere-pixel tex 0.5 0.5 1.0 0.01 {})))))

  (testing "facing the viewer, a sphere's surface frame is the plane's, bumps and all"
    ;; Its nearest point is z = 1, moved here onto the plane z = 0.
    (let [tex (tx/cobblestone {})]
      (is (close? (tx/plane-pixel tex 0.0 0.0 0.01 {})
                  (tx/sphere-pixel (fn [x y z] (tex x y (- z 1.0))) 0.0 0.0 1.0 0.01 {}))))))
