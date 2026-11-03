(ns allgo.language.parsing
  "Natural language: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 23.

    CYK parsing      CYK-PARSE (figure 23.5): the most probable derivation
                     of every span of the words from every nonterminal of
                     a probabilistic grammar in Chomsky normal form, spans
                     filled shortest first, a span a step; and the parse
                     tree read back from it
    sentence trees   bracketed trees as the Penn Treebank writes them
                     (figure 23.8), read and written
    generation       trees sampled from a probabilistic grammar, and
                     their sentences, each word with the part of speech it
                     was drawn as: tagged data whose tags are right by
                     construction, for training taggers
    learning         a grammar from a treebank, by counting each rule's
                     share of its left-hand side's nodes; and from unparsed
                     sentences, by the inside-outside algorithm -- EM, the
                     expected uses of each rule found from inside and
                     outside probabilities, an iteration a step

  A grammar is `{:rules [[X Y Z p] ...] :lexicon {word [[X p] ...]}}`:
  binary rules X -> Y Z with probability p, and lexical rules X -> word.
  A tree is `[label child ...]`, a leaf a string. `e-prob-chomsky` is
  aima-python's (MIT) small grammar of that shape."
  (:require [clojure.math :as math]
            [clojure.string :as str]))

(def e-prob-chomsky
  "aima-python's `E_Prob_Chomsky`."
  {:rules [[:S :NP :VP 1.0]
           [:NP :Article :Noun 0.6] [:NP :Adjective :Noun 0.4]
           [:VP :Verb :NP 0.5] [:VP :Verb :Adjective 0.5]]
   :lexicon {"the" [[:Article 0.5]] "a" [[:Article 0.25]] "an" [[:Article 0.25]]
             "robot" [[:Noun 0.4]] "sheep" [[:Noun 0.4]] "fence" [[:Noun 0.2]]
             "good" [[:Adjective 0.5]] "new" [[:Adjective 0.2]] "sad" [[:Adjective 0.3]]
             "is" [[:Verb 0.5]] "say" [[:Verb 0.3]] "are" [[:Verb 0.2]]}})

;; ---------------------------------------------------------------------------
;; CYK

(defn cyk-seq
  "Every state of CYK-PARSE on `words`: `{:P :back :span}`, `:P` mapping
  `[X start length]` to the probability of the most probable X over that
  span, `:back` how it was made; the first state has the lexical rules
  filled in, and each step fills one span, shortest spans first."
  [words {:keys [rules lexicon]}]
  (let [n (count words)
        lexical (reduce (fn [st [i w]]
                          (reduce (fn [st [x p]]
                                    (-> st (assoc-in [:P [x i 1]] p) (assoc-in [:back [x i 1]] w)))
                                  st (lexicon w)))
                        {:P {} :back {}}
                        (map-indexed vector words))
        spans (for [length (range 2 (inc n)) start (range 0 (inc (- n length)))] [start length])]
    (reductions (fn [{:keys [P] :as st} [start length]]
                  (reduce (fn [st [x y z p len1]]
                            (let [len2 (- length len1)
                                  q (* (get P [y start len1] 0.0) (get P [z (+ start len1) len2] 0.0) p)]
                              (if (> q (get-in st [:P [x start length]] 0.0))
                                (-> st (assoc-in [:P [x start length]] q)
                                    (assoc-in [:back [x start length]] [y z len1]))
                                st)))
                          (assoc st :span [start length])
                          (for [len1 (range 1 length) [x y z p] rules] [x y z p len1])))
                lexical spans)))

(defn cyk-parse
  "The CYK table of `words`: `{:P :back}`."
  [words grammar]
  (last (cyk-seq words grammar)))

(defn parse-tree
  "The most probable tree for the whole of `words` rooted at `root`
  (default :S), and its probability: `[tree p]`, or nil."
  ([words grammar] (parse-tree words grammar :S))
  ([words grammar root]
   (let [{:keys [P back]} (cyk-parse words grammar)
         n (count words)]
     (letfn [(build [x start length]
               (let [b (back [x start length])]
                 (if (string? b)
                   [x b]
                   (let [[y z len1] b]
                     [x (build y start len1) (build z (+ start len1) (- length len1))]))))]
       (when-let [p (P [root 0 n])]
         [(build root 0 n) p])))))

;; ---------------------------------------------------------------------------
;; Bracketed trees

(defn read-tree
  "A tree from Penn Treebank bracketing: `[S [NP-SBJ Her eyes] ...]`, each
  bracket a label and its children, words the leaves. Either bracket
  shape, square or round."
  [text]
  (let [toks (re-seq #"[\[\]()]|[^\s\[\]()]+" text)
        pos (volatile! 0)
        next-tok #(let [t (nth toks @pos nil)] (vswap! pos inc) t)
        peek-tok #(nth toks @pos nil)]
    (letfn [(node []
              (let [t (next-tok)]
                (if (#{"[" "("} t)
                  (let [label (when-not (#{"[" "("} (peek-tok)) (next-tok))]
                    (loop [kids []]
                      (if (#{"]" ")"} (peek-tok))
                        (do (next-tok) (into (if label [label] [nil]) kids))
                        (recur (conj kids (node))))))
                  t)))]
      (node))))

(defn write-tree
  "A tree in bracketed form, on one line."
  [tree]
  (if (string? tree)
    tree
    (let [[label & kids] tree]
      (str "[" (when label (str (if (keyword? label) (name label) label) " ")) (str/join " " (map write-tree kids)) "]"))))

(defn leaves
  "The words of a tree, in order."
  [tree]
  (if (string? tree) [tree] (vec (mapcat leaves (rest tree)))))

;; ---------------------------------------------------------------------------
;; Generation

(def ambiguous-grammar
  "A small grammar of this library's whose words are often more than one
  part of speech -- \"fish\", \"duck\", \"watch\" and \"saw\" nouns or
  verbs, \"can\" a noun or a modal -- so that tagging needs the context,
  not just the word. Nonterminals other than the start are free to
  appear only on the right; the lexicon gives P(word | tag)."
  {:rules [[:S :NP :VP 1.0]
           [:NP :Det :N 0.7] [:NP :Adj :N 0.3]
           [:VP :V :NP 0.5] [:VP :Modal :VB 0.3] [:VP :V :Adv 0.2]
           [:VB :V :NP 1.0]]
   :lexicon {"the" [[:Det 0.6]] "a" [[:Det 0.4]]
             "fish" [[:N 0.2] [:V 0.2]] "duck" [[:N 0.15] [:V 0.15]] "watch" [[:N 0.15] [:V 0.2]]
             "saw" [[:N 0.1] [:V 0.15]] "can" [[:N 0.1] [:Modal 0.7]]
             "dog" [[:N 0.15]] "cat" [[:N 0.15]] "sees" [[:V 0.15]] "likes" [[:V 0.15]]
             "young" [[:Adj 0.5]] "old" [[:Adj 0.5]] "will" [[:Modal 0.3]]
             "quickly" [[:Adv 0.5]] "well" [[:Adv 0.5]]}})

(defn- pick-weighted [rng choices]
  (let [total (reduce + (map second choices))
        r (* total (double (rng)))]
    (loop [[[x p] & more] choices acc 0.0]
      (if (or (empty? more) (< r (+ acc p))) x (recur more (+ acc p))))))

(defn generate-tree
  "A tree sampled from `grammar` from `start` (default :S): each
  nonterminal expanded by one of its rules chosen by probability, each
  part of speech by a word chosen by P(word | tag). `[X left right]`, a
  part of speech `[X word]`."
  ([grammar rng] (generate-tree grammar rng :S))
  ([{:keys [rules lexicon]} rng start]
   (let [by-head (group-by first rules)
         words-of (reduce (fn [m [w tags]] (reduce (fn [m [t p]] (update m t (fnil conj []) [w p])) m tags)) {} lexicon)]
     (letfn [(expand [x]
               (if-let [rs (by-head x)]
                 (let [[_ y z] (pick-weighted rng (map (fn [r] [r (peek r)]) rs))]
                   [x (expand y) (expand z)])
                 [x (pick-weighted rng (words-of x))]))]
       (expand start)))))

(defn tagged-leaves
  "A tree's words, each with the part of speech above it: `[[word tag]
  ...]`."
  [tree]
  (let [[label & kids] tree]
    (if (and (= 1 (count kids)) (string? (first kids)))
      [[(first kids) label]]
      (vec (mapcat tagged-leaves kids)))))

(defn generate
  "A sentence sampled from `grammar` from `start` (default :S), each word
  with the part of speech it was drawn as: `[[word tag] ...]`."
  ([grammar rng] (generate grammar rng :S))
  ([grammar rng start] (tagged-leaves (generate-tree grammar rng start))))

;; ---------------------------------------------------------------------------
;; Learning a grammar

(defn- normalize-grammar
  "The grammar whose rules, binary and lexical, have probabilities
  proportional to `counts` -- `{[X Y Z] c}` and `{[X word] c}` -- those of
  each left-hand side summing to 1."
  [binary lexical]
  (let [totals (reduce (fn [m [[x] c]] (update m x (fnil + 0.0) c)) {} (concat binary lexical))]
    {:rules (vec (for [[[x y z] c] (sort-by (comp str key) binary) :when (pos? (double c))] [x y z (/ (double c) (totals x))]))
     :lexicon (reduce (fn [m [[x w] c]] (if (pos? (double c)) (update m w (fnil conj []) [x (/ (double c) (totals x))]) m))
                      {} (sort-by (comp str key) lexical))}))

(defn learn-grammar
  "A probabilistic grammar from a treebank, `trees` in Chomsky normal
  form, by counting, as the chapter describes: each rule X -> Y Z or X ->
  word given the fraction of the X nodes it makes."
  [trees]
  (letfn [(nodes [t] (if (string? t) [] (cons t (mapcat nodes (rest t)))))]
    (let [all (mapcat nodes trees)
          lexical? (fn [[_ k]] (string? k))]
      (normalize-grammar (frequencies (map (fn [[x [y] [z]]] [x y z]) (remove lexical? all)))
                         (frequencies (map (fn [[x w]] [x w]) (filter lexical? all)))))))

(defn inside
  "The inside probabilities of `words`: `{[X start length] p}`, the
  probability that X derives exactly that span -- CYK with a sum where it
  has a max."
  [words {:keys [rules lexicon]}]
  (let [n (count words)
        lexical (into {} (for [[i w] (map-indexed vector words) [x p] (lexicon w)] [[x i 1] p]))]
    (reduce (fn [beta [start length]]
              (reduce (fn [beta [x y z p len1]]
                        (let [q (* (get beta [y start len1] 0.0) (get beta [z (+ start len1) (- length len1)] 0.0) (double p))]
                          (if (pos? q) (update beta [x start length] (fnil + 0.0) q) beta)))
                      beta
                      (for [len1 (range 1 length) [x y z p] rules] [x y z p len1])))
            lexical
            (for [length (range 2 (inc n)) start (range 0 (inc (- n length)))] [start length]))))

(defn outside
  "The outside probabilities of `words` given their inside ones `beta`:
  `{[X start length] p}`, the probability of everything else in the
  sentence with X spanning that, from :S over the whole -- the longest
  spans first, each passing its share down to the two parts it splits
  into."
  [words {:keys [rules]} beta root]
  (let [n (count words)]
    (reduce (fn [alpha [start length]]
              (reduce (fn [alpha [x y z p len1]]
                        (let [a (get alpha [x start length] 0.0)
                              ly [y start len1] lz [z (+ start len1) (- length len1)]
                              by (get beta ly 0.0) bz (get beta lz 0.0)]
                          (if (and (pos? a) (pos? (* by bz)))
                            (-> alpha
                                (update ly (fnil + 0.0) (* a (double p) bz))
                                (update lz (fnil + 0.0) (* a (double p) by)))
                            alpha)))
                      alpha
                      (for [len1 (range 1 length) [x y z p] rules] [x y z p len1])))
            {[root 0 n] 1.0}
            (for [length (range n 1 -1) start (range 0 (inc (- n length)))] [start length]))))

(defn- expected-counts
  "The E step over `sentences`: each rule's expected number of uses,
  `[binary lexical]`, and the sentences' total log likelihood. Sentences
  the grammar cannot derive are left out."
  [sentences {:keys [rules lexicon] :as grammar} root]
  (reduce (fn [[binary lexical ll] words]
            (let [n (count words)
                  beta (inside words grammar)
                  total (get beta [root 0 n] 0.0)]
              (if-not (pos? total)
                [binary lexical ll]
                (let [alpha (outside words grammar beta root)
                      binary (reduce (fn [m [x y z p]]
                                       (update m [x y z] (fnil + 0.0)
                                               (/ (reduce + (for [length (range 2 (inc n)) start (range 0 (inc (- n length))) len1 (range 1 length)]
                                                              (* (get alpha [x start length] 0.0) (double p)
                                                                 (get beta [y start len1] 0.0) (get beta [z (+ start len1) (- length len1)] 0.0))))
                                                  total)))
                                     binary rules)
                      lexical (reduce (fn [m [i w]]
                                        (reduce (fn [m [x p]] (update m [x w] (fnil + 0.0) (/ (* (get alpha [x i 1] 0.0) (double p)) total)))
                                                m (lexicon w)))
                                      lexical (map-indexed vector words))]
                  [binary lexical (+ ll (math/log total))]))))
          [{} {} 0.0] sentences))

(defn inside-outside-seq
  "The inside-outside algorithm (Baker, 1979): EM for the rule
  probabilities of `grammar` from unparsed `sentences`, as Baum-Welch is
  EM for an HMM's. The E step finds each rule's expected number of uses
  in each sentence's parses -- inside probabilities from the words up,
  outside ones from the root down -- and the M step makes each rule's
  probability its share of its left-hand side's. The rules themselves are
  fixed: one with probability 0 stays 0. An iteration a step, each state
  `{:grammar :log-likelihood}`, the likelihood the sentences have under
  that grammar, which never falls. `root` defaults to :S."
  ([grammar sentences] (inside-outside-seq grammar sentences :S))
  ([grammar sentences root]
   (letfn [(states [g]
             (lazy-seq
              (let [[binary lexical ll] (expected-counts sentences g root)]
                (cons {:grammar g :log-likelihood ll}
                      (states (normalize-grammar binary lexical))))))]
     (states grammar))))
