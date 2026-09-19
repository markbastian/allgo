(ns allgo.settlement-test
  (:require [allgo.procedural.island :as island]
            [allgo.procedural.settlement :as st]
            [clojure.test :refer [deftest is testing]]))

(defn- seeded [seed]
  (let [r (java.util.Random. seed)]
    (fn [] (.nextDouble r))))

(defn- world
  ([] (world {}))
  ([opts]
   (st/populate (island/generate {:points 500 :seed 5 :rng (seeded 5)}) opts)))

(def ^:private m (delay (world)))

(defn- hops
  "Cells within `n` land steps of `start`."
  [{:keys [centers]} start n]
  (loop [frontier #{start} seen #{start} i 0]
    (if (or (= i n) (empty? frontier))
      seen
      (let [nf (into #{} (comp (mapcat #(:neighbors (centers %)))
                               (remove seen)
                               (remove #(:water? (centers %))))
                     frontier)]
        (recur nf (into seen nf) (inc i))))))

(deftest towns-test
  (let [{:keys [centers towns] :as w} @m]

    (testing "there are towns, and not more than were asked for"
      (is (seq towns))
      (is (<= (count towns) (:towns st/defaults))))

    (testing "nobody lives in the sea"
      (is (every? #(not (:water? (centers %))) towns))
      (is (every? #(:town? (centers %)) towns)))

    (testing "and not on top of each other"
      ;; The refusal rule. Without it every town lands in the best river
      ;; valley, because the second best site in it still beats anywhere
      ;; else on the map.
      (let [spacing (:spacing st/defaults)]
        (is (every? (fn [t]
                      (let [near (disj (hops w t (dec spacing)) t)]
                        (not-any? (set towns) near)))
                    towns))))

    (testing "they choose fresh water and low ground over the rest"
      ;; What `appeal` is for. Not every town will be on a river, but the
      ;; towns should be markedly lower than the land at large.
      (let [land (remove :water? centers)
            mean-elev (fn [xs] (/ (reduce + (map :elevation xs)) (max 1 (count xs))))]
        (is (< (mean-elev (map centers towns)) (mean-elev land)))))))

(deftest roads-test
  (let [{:keys [centers edges roads towns]} @m]

    (testing "a spanning tree, not a mesh"
      ;; One link per town less one, when the towns are all on one
      ;; landmass. More would mean parallel roads down the same valley.
      (is (<= (count roads) (dec (count towns)))))

    (testing "roads are on land and join cells that really are neighbours"
      (let [road-edges (filter #(pos? (long (:road % 0))) edges)]
        (is (seq road-edges))
        (is (every? (fn [e]
                      (let [[a b] (:centers e)]
                        (and a b
                             (not (:water? (centers a)))
                             (not (:water? (centers b)))
                             (some #{b} (:neighbors (centers a))))))
                    road-edges))))

    (testing "every link is a connected route between two towns"
      ;; The road edges of a link have to form a path. Checked by walking
      ;; the road subgraph from one town and requiring the other to be in
      ;; the same component.
      (let [road-adj (reduce (fn [m e]
                               (if (pos? (long (:road e 0)))
                                 (let [[a b] (:centers e)]
                                   (-> m (update a (fnil conj #{}) b)
                                       (update b (fnil conj #{}) a)))
                                 m))
                             {} edges)
            reach (fn [s] (loop [f #{s} seen #{s}]
                            (if (empty? f)
                              seen
                              (let [nf (into #{} (comp (mapcat road-adj) (remove seen)) f)]
                                (recur nf (into seen nf))))))]
        (is (every? (fn [link]
                      (let [[a b] (vec link)]
                        (contains? (reach a) b)))
                    roads))))))

(deftest territory-test
  (let [{:keys [centers towns]} @m
        land (remove :water? centers)]

    (testing "a town holds itself"
      (is (every? (fn [t] (= t (:territory (centers t)))) towns)))

    (testing "every claim names a real town"
      (is (every? (set towns) (keep :territory land))))

    (testing "the land is claimed, bar what cannot be walked to"
      ;; An islet with no town on it has no owner, which is right.
      (is (> (count (filter :territory land)) (* 0.9 (count land)))))

    (testing "water is claimed by nobody"
      (is (every? #(nil? (:territory %)) (filter :water? centers))))

    (testing "reach rises away from the town that holds you"
      (is (every? #(zero? (:reach (centers %))) towns))
      (is (every? #(>= (:reach %) 0.0) (filter :territory land))))

    (testing "borders end up on the high ground"
      ;; The claim the namespace is built on, and the difference between
      ;; assigning territory by travel cost and by distance. A ridge is
      ;; expensive from both sides, so neither town reaches over it, and
      ;; the frontier settles there without anything looking for one.
      ;; Assigned by distance this would be a Voronoi diagram of the
      ;; towns and the two figures would match.
      (let [owned (filter :territory land)
            frontier? (fn [c]
                        (some (fn [n]
                                (let [o (centers n)]
                                  (and (:territory o) (not= (:territory o) (:territory c)))))
                              (:neighbors c)))
            mean (fn [xs] (/ (reduce + (map :elevation xs)) (max 1 (count xs))))
            on-border (filter frontier? owned)
            inside (remove frontier? owned)]
        (is (seq on-border))
        (is (> (mean on-border) (* 1.2 (mean inside)))
            (str "borders " (mean on-border) " against interiors " (mean inside)))))))

(deftest options-test
  (testing "asking for fewer towns gives fewer"
    (is (<= (count (:towns (world {:towns 5}))) 5)))

  (testing "asking for more spacing spreads them further"
    (let [tight (world {:towns 30 :spacing 1})
          loose (world {:towns 30 :spacing 5})]
      (is (> (count (:towns tight)) (count (:towns loose)))))))

(deftest determinism-test
  (testing "the same island gives the same towns, roads and borders"
    (let [a (world) b (world)]
      (is (= (:towns a) (:towns b)))
      (is (= (:roads a) (:roads b)))
      (is (= (map :territory (:centers a)) (map :territory (:centers b)))))))
