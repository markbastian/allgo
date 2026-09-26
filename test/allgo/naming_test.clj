(ns allgo.naming-test
  (:require [allgo.procedural.naming :as nm]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]))

(defn- names-of [lang n] (mapv #(nm/name-for lang %) (range n)))

(deftest language-test
  (testing "a language is a small inventory, not the whole pool"
    ;; The narrowness is the point. A language holding every consonant is
    ;; the same bag everyone else draws from, and its names will not
    ;; sound like they belong together.
    (doseq [seed (range 12)]
      (let [{:keys [onsets nuclei codas]} (nm/language seed)]
        (is (<= 6 (count onsets) 18) (str "seed " seed))
        (is (<= 3 (count nuclei) 8) (str "seed " seed))
        (is (<= 2 (count codas) 6) (str "seed " seed))
        (is (apply distinct? onsets))
        (is (apply distinct? nuclei))
        (is (apply distinct? codas)))))

  (testing "different seeds are different languages"
    (is (not= (nm/language 1) (nm/language 2)))
    (is (= (nm/language 1) (nm/language 1)))))

(deftest word-test
  (let [lang (nm/language 5)
        ns (names-of lang 300)]

    (testing "every name is a capitalized, pronounceable-length word"
      (is (every? seq ns))
      (is (every? #(re-matches #"[A-Z][a-z]+" %) ns))
      (is (every? #(<= 2 (count %) 18) ns)))

    (testing "no name has a pile of consonants in it"
      ;; Four in a row is always an artifact of gluing syllables, never a
      ;; cluster any of these languages has. Three is allowed, because
      ;; `thr` and `str` are real onsets in the pool.
      (is (not-any? #(re-find #"[^aeiouy]{4}" (str/lower-case %)) ns)))

    (testing "and no doubled vowel, which is only ever a join"
      (is (not-any? #(re-find #"([aeiouy])\1" (str/lower-case %)) ns)))))

(deftest family-test
  (testing "a language's names are built only from its own sounds"
    ;; The whole claim of the namespace, stated as an invariant: if a
    ;; language never drew `k`, none of its places have a `k` in them.
    ;; This is what makes two towns in one culture sound related and two
    ;; across a border sound foreign.
    (doseq [seed (range 8)]
      (let [{:keys [onsets nuclei codas] :as lang} (nm/language seed)
            allowed (set (mapcat seq (concat onsets nuclei codas)))]
        (is (every? (fn [w] (every? allowed (str/lower-case w)))
                    (names-of lang 200))
            (str "seed " seed " used a letter it does not have")))))

  (testing "two languages really do sound different"
    ;; Not a strong statement about phonetics -- just that the letter
    ;; distributions are not the same, which they would be if the
    ;; inventories were not being drawn per language.
    (let [letters (fn [seed]
                    (set (mapcat (comp seq str/lower-case)
                                 (names-of (nm/language seed) 120))))]
      (is (not= (letters 1) (letters 2)))
      (is (not= (letters 3) (letters 7))))))

(deftest keying-test
  (let [lang (nm/language 9)]

    (testing "a name is a function of its key and nothing else"
      ;; So that adding a town, or naming them in another order, does not
      ;; rename everything else on the map.
      (is (= (nm/name-for lang 42) (nm/name-for lang 42)))
      (is (= (mapv #(nm/name-for lang %) [7 3 11])
             (mapv #(nm/name-for lang %) [7 3 11])))
      (is (= (nm/name-for lang 5)
             (last (mapv #(nm/name-for lang %) (range 6))))))

    (testing "different keys give different names, mostly"
      ;; Collisions are inevitable from a small inventory and are not a
      ;; bug; a map of five hundred places wants most of them distinct.
      (is (> (count (distinct (names-of lang 500))) 400)))))
