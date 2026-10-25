(ns allgo.gesture-test
  (:require [allgo.geometry.vec3 :as v]
            [allgo.procedural.gesture :as gesture]
            [allgo.procedural.noise :as noise]
            [clojure.test :refer [deftest is testing]]))

(defn- max-change
  "The largest difference in any angle between two poses."
  [a b]
  (apply max (for [[j x] a] (apply max (map #(abs (- (double %1) (double %2))) x (get b j))))))

(def ^:private rest-lengths
  (let [f (gesture/figure {})]
    (into {} (for [[a b] gesture/bones] [[a b] (v/distance (f a) (f b))]))))

(deftest figure-test
  (testing "at rest the figure stands: feet on the floor, head above the hips"
    (let [f (gesture/figure {})]
      (is (< (abs (- 0.08 (second (f :l-ankle)))) 1e-9))
      (is (> (second (f :crown)) 1.7))
      (is (> (second (f :l-knee)) (second (f :l-ankle))))))
  (testing "in any action, bones keep their lengths and feet stay put"
    (doseq [[layer names] gesture/actions
            name (keys names)
            t [0.3 1.7 4.2]]
      (let [perf (gesture/start 0.0 {layer name})
            f (gesture/figure (gesture/pose perf t {}))]
        (doseq [[[a b] len] rest-lengths]
          (is (< (abs (- len (v/distance (f a) (f b)))) 1e-6) [name a b]))
        (is (< (v/distance (f :l-ankle) [0.1 0.08 0.0]) 1e-9) name)
        (is (< (v/distance (f :r-ankle) [-0.1 0.08 0.0]) 1e-9) name))))
  (testing "lowering the hips bends the knees forward"
    (let [up (gesture/figure {})
          down (gesture/figure {:root [0.0 -0.1 0.0]})]
      (is (> (nth (down :l-knee) 2) (nth (up :l-knee) 2) 0.0)))))

(deftest action-test
  (testing "every angle stays in its range"
    (let [basis (noise/gradient-basis {:seed 3})]
      (doseq [[_ names] gesture/actions
              [name action] names
              t (range 0.0 20.0 0.37)
              [j values] (gesture/evaluate action t basis)
              [dof value] (map vector (get-in action [:joints j]) values)]
        (let [[lo hi] (if (number? dof) [dof dof] dof)
              k (if (= j :root) 1.0 (/ Math/PI 180.0))]
          (is (<= (- (* k (min lo hi)) 1e-9) value (+ (* k (max lo hi)) 1e-9)) [name j])))))
  (testing "a sine repeats with its rate; noise does not"
    (let [basis (noise/gradient-basis {:seed 0})]
      (is (< (abs (- (gesture/signal [:sin 0.2] 0.3 2.0 basis) (gesture/signal [:sin 0.2] 0.8 2.0 basis))) 1e-12))
      (is (> (abs (- (gesture/signal [:noise 4] 0.3 2.0 basis) (gesture/signal [:noise 4] 0.8 2.0 basis))) 1e-3))))
  (testing "the same seed gives the same performance, another seed another"
    (let [perf (gesture/start 0.0 {:body "fidget"})]
      (is (= (gesture/pose perf 3.3 {:seed 5}) (gesture/pose perf 3.3 {:seed 5})))
      (is (pos? (max-change (gesture/pose perf 3.3 {:seed 5}) (gesture/pose perf 3.3 {:seed 6})))))))

(deftest blend-test
  (let [opts {:transition 1.0}
        perf (gesture/start 0.0 {:body "rest"})
        switched (gesture/change perf :body "dance" 2.0 1.0)]
    (testing "a change starts where the old action was, and ends at the new one"
      (is (< (max-change (gesture/pose perf 2.0 opts) (gesture/pose switched 2.0 opts)) 1e-12))
      (is (< (max-change (gesture/pose (gesture/start 0.0 {:body "dance"}) 3.5 opts)
                         (gesture/pose switched 3.5 opts))
             1e-12)))
    (testing "and gets there smoothly"
      (is (every? #(< (max-change (gesture/pose switched % opts) (gesture/pose switched (+ % 0.01) opts)) 0.06)
                  (range 1.9 3.2 0.01))))
    (testing "a change in the middle of a change does not jump either"
      (let [again (gesture/change switched :body "tired" 2.4 1.0)]
        (is (< (max-change (gesture/pose switched 2.4 opts) (gesture/pose again 2.4 opts)) 1e-12))))
    (testing "an upper layer takes only its own joints, and the body shows through"
      (let [waving (gesture/change perf :arms "wave" 0.0 1.0)
            a (gesture/pose perf 5.0 opts)
            b (gesture/pose waving 5.0 opts)]
        (is (= (:chest a) (:chest b)))
        (is (= (:l-shoulder a) (:l-shoulder b)))
        (is (not= (:r-shoulder a) (:r-shoulder b)))))))
