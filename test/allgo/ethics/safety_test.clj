(ns allgo.ethics.safety-test
  (:require [allgo.ethics.safety :as safety]
            [clojure.test :refer [deftest is testing]]))

;; AI safety, after Russell and Norvig, *Artificial Intelligence: A Modern
;; Approach*, 4th edition, section 27.3.7.

(def ^:private service
  "An online service, down if both servers fail, or the network does, or
  the power and then the backup generator."
  {:or [{:and [{:event :server-a} {:event :server-b}]}
        {:event :network}
        {:and [{:event :power} {:event :generator}]}]})

(deftest fault-tree-test
  (testing "the chance the service goes down: 1 - (1 - 0.01 x 0.02)(1 - 0.001)(1 - 0.05 x 0.1)"
    (is (< (abs (- (- 1.0 (* (- 1.0 0.0002) (- 1.0 0.001) (- 1.0 0.005)))
                   (safety/failure-probability service {:server-a 0.01 :server-b 0.02 :network 0.001 :power 0.05 :generator 0.1})))
           1e-15)))
  (testing "its minimal cut sets"
    (is (= #{#{:server-a :server-b} #{:network} #{:power :generator}} (safety/cut-sets service))))
  (testing "a cause shared by both branches of an AND is a cut set on its own, and the larger set holding it is not minimal"
    (is (= #{#{:a} #{:b :c}} (safety/cut-sets {:and [{:or [{:event :a} {:event :b}]} {:or [{:event :a} {:event :c}]}]})))))

(deftest fmea-test
  (testing "failure modes ranked by severity x occurrence x detection"
    (is (= [:sensor :bolt :fan]
           (map :mode (safety/fmea [{:mode :bolt :severity 9 :occurrence 2 :detection 3}
                                    {:mode :fan :severity 3 :occurrence 4 :detection 2}
                                    {:mode :sensor :severity 7 :occurrence 3 :detection 8}]))))))

(deftest low-impact-test
  (let [rows ["#######" "#.....#" "#S.V.C#" "#.....#" "#######"]]
    (testing "told only to fetch the coffee, the agent walks straight through the vase in its way"
      (let [{:keys [cells broken]} (safety/route (safety/coffee-world rows -0.04 0.0 0.95))]
        (is (= 5 (count cells)))
        (is (= #{[3 2]} broken))))
    (testing "charged for whatever it changes, it goes round -- two steps longer, nothing broken"
      (let [{:keys [cells broken]} (safety/route (safety/coffee-world rows -0.04 0.5 0.95))]
        (is (= 7 (count cells)))
        (is (empty? broken))
        (is (= [5 2] (peek cells)))))))
