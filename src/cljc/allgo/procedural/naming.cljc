(ns allgo.procedural.naming
  "Invented names, by inventing the language first.

  The cheap way to name a thousand places is to draw syllables out of one
  bag, and it reads as exactly that: every name is plausible on its own
  and the set has no character, because there is nothing for a name to be
  characteristic *of*. Two towns a valley apart come out sounding like
  they were founded by different species.

  So this generates a language and then generates names out of it. A
  language here is a small phonology -- which consonants it uses, which
  vowels, what a syllable is allowed to look like, how many syllables a
  word runs to -- drawn once from a seed out of a much larger pool. Every
  name from one language shares that inventory, and the family
  resemblance is a consequence rather than an effect: if a language never
  drew `k`, none of its towns have a `k` in them.

  That is the whole trick, and it is Azgaar's: names belong to cultures,
  cultures belong to regions, so neighbors sound like neighbors and a
  border is audible.

  ## Seeds, not streams

  `name-for` takes the thing being named -- a cell id, a river mouth --
  and derives its own stream from that plus the language's seed. So a
  town's name does not depend on how many names were generated before it,
  which means adding a town, or naming them in a different order, does
  not rename everything else on the map.

  The stream is `allgo.random`'s Park-Miller, for the reason
  given there: it is arithmetic on doubles under 2^53, so the JVM and a
  browser agree about what a seed means. A name generator that produced
  different names on the server and the client would be a strange thing
  to debug."
  (:require [allgo.random :as random]
            [clojure.string :as str]))

;; ---------------------------------------------------------------------------
;; The pools a language is drawn from
;;
;; Deliberately wider than any one language uses. The pools are the space
;; of possible languages; a language is a handful of picks out of them,
;; and it is the narrowness of the picks that gives it a character.

(def ^:private simple-onsets
  ["b" "d" "f" "g" "h" "j" "k" "l" "m" "n" "p" "r" "s" "t" "v" "w" "z"
   "th" "sh" "ch" "kh" "gh" "ph" "zh"])

(def ^:private cluster-onsets
  ["br" "dr" "fr" "gr" "kr" "pr" "tr" "bl" "fl" "gl" "kl" "pl" "sl"
   "st" "sk" "sp" "sn" "sm" "sw" "tw" "thr" "shr" "str"])

(def ^:private simple-nuclei ["a" "e" "i" "o" "u" "y"])

(def ^:private diphthongs
  ["ae" "ai" "au" "ea" "ei" "eo" "ia" "ie" "io" "oa" "oi" "ou" "ua"])

(def ^:private codas
  ["n" "r" "l" "s" "m" "th" "k" "t" "d" "ng" "rn" "rk" "ld" "lm" "st"
   "sk" "nd" "nt" "rd" "sh" "ch" "ss" "ll" "rr"])

(defn- pick-n
  "`n` distinct items, in a stable order."
  [rng n coll]
  (loop [taken [] pool (vec coll) i 0]
    (if (or (= (count taken) (long n)) (empty? pool) (> i 200))
      taken
      (let [k (min (dec (count pool)) (long (* (rng) (count pool))))]
        (recur (conj taken (nth pool k))
               (into (subvec pool 0 k) (subvec pool (inc k)))
               (inc i))))))

(defn language
  "A phonology, drawn from `seed`.

  The inventories are small on purpose. A language with twenty consonants
  and every vowel is the same bag of syllables everyone else is drawing
  from, and its names will not sound like they belong together. Six to
  eleven consonants and three to five vowels is enough to build a
  thousand distinguishable names and few enough that they rhyme with each
  other."
  [seed]
  (let [rng (random/rng seed 1 1 977)
        cluster? (< (rng) 0.55)
        simple (pick-n rng (+ 6 (long (* 6 (rng)))) simple-onsets)
        onsets (into simple
                     (when cluster? (pick-n rng (+ 1 (long (* 4 (rng)))) cluster-onsets)))
        ;; Always mostly plain vowels, with diphthongs as seasoning. A
        ;; language whose whole inventory is `oi`, `au` and `ia` makes
        ;; every name a variation on the same noise, and drawing the
        ;; vowels from one pool lets that happen surprisingly often.
        vowels (into (pick-n rng (+ 3 (long (* 3 (rng)))) simple-nuclei)
                     (pick-n rng (long (* 3 (rng))) diphthongs))
        finals (pick-n rng (+ 2 (long (* 3 (rng)))) codas)
        ;; How often a syllable closes. A language that always closes is
        ;; clotted, one that never does is all open vowels; either extreme
        ;; is the same word over and over.
        coda-rate (+ 0.15 (* 0.5 (rng)))
        ;; Whether words may start on a vowel at all. Having this vary is
        ;; most of what makes two languages sound unrelated.
        onset-rate (+ 0.75 (* 0.25 (rng)))]
    {:seed seed
     :onsets (vec onsets)
     ;; Kept apart so a syllable following a closed one can reach for a
     ;; single consonant instead of a cluster. `plain` is single letters
     ;; only -- `sk` meeting `th` is still four consonants even though
     ;; both are "simple".
     :simple-onsets (vec simple)
     :plain-onsets (let [p (filterv #(= 1 (count %)) simple)]
                     (if (seq p) p (vec simple)))
     :nuclei (vec vowels)
     :codas (vec finals)
     :coda-rate coda-rate
     :onset-rate onset-rate
     :min-syllables 2
     :max-syllables (+ 2 (long (* 2 (rng))))}))

(defn- syllable
  "One syllable, and whether it closed.

  `after-coda?` is the whole reason this returns a pair. A closed
  syllable followed by a cluster onset gives `...rk` + `thr...`, which is
  five consonants in a row and not a word anybody would say. After a
  coda the next syllable takes a single consonant and is unlikely to
  close again."
  [{:keys [onsets plain-onsets nuclei codas coda-rate onset-rate]}
   rng first? after-coda? previous-onset]
  (let [pool (if after-coda? plain-onsets onsets)
        want-onset? (or (not first?) (< (rng) onset-rate))
        onset (when want-onset?
                ;; One retry when the draw repeats the last onset.
                ;; `Thiththi` is what happens without it, and a small
                ;; inventory makes the collision common rather than rare.
                (let [o (random/pick rng pool)]
                  (if (= o previous-onset) (random/pick rng pool) o)))
        nucleus (random/pick rng nuclei)
        close? (< (rng) (if after-coda? (* 0.3 coda-rate) coda-rate))]
    [(str onset nucleus (when close? (random/pick rng codas))) close? onset]))

(defn- tidy
  "Smooths the joins a syllable-at-a-time build leaves behind.

  Three of the same letter in a row is the common one, and a doubled
  vowel pair across a syllable boundary is the ugly one. Neither is a
  sound; both are an artifact of gluing."
  [w]
  (-> w
      (str/replace #"([a-z])\1{2,}" "$1$1")
      ;; A doubled vowel is always a join artifact, never a sound here:
      ;; `Phaarnpho` is `pha` meeting `arn`, and no language in the pools
      ;; has a long vowel to spell that way.
      (str/replace #"([aeiouy])\1" "$1")
      (str/replace #"([^aeiouy])\1([^aeiouy])" "$1$2")
      ;; Four consonants in a row is always a join and never a cluster
      ;; any of these languages has. Three is left alone, because `thr`
      ;; and `str` are real onsets in the pool.
      (str/replace #"([^aeiouy]{2})[^aeiouy]{2,}" "$1")))

(defn word
  "One word of `lang`, as its own repeatable function of `k`.

  Two different `k` give two different words; the same `k` always gives
  the same one, however many other names have been drawn."
  [{:keys [min-syllables max-syllables] :as lang} k]
  (let [rng (random/rng (:seed lang) (inc (long k)) 7 977)
        n (+ (long min-syllables)
             (long (* (rng) (inc (- (long max-syllables) (long min-syllables))))))
        raw (loop [i 0 acc "" after-coda? false prev nil]
              (if (= i (max 1 n))
                acc
                (let [[syl close? onset] (syllable lang rng (zero? i) after-coda? prev)]
                  (recur (inc i) (str acc syl) close? onset))))
        w (tidy raw)]
    (str (str/upper-case (subs w 0 1)) (subs w 1))))

(defn compound
  "Two short words of `lang` run together, for the places that want a
  longer name than a single word gives.

  Both halves are cut to a syllable or two first. Gluing two full-length
  words gives a six-syllable name that reads as a mistake rather than as
  a place."
  [lang k]
  (let [short (assoc lang :min-syllables 1 :max-syllables 2)
        a (word short k)
        b (str/lower-case (word short (+ 9973 (long k))))]
    (tidy (str a b))))

(defn name-for
  "A name for the thing keyed by `k`, occasionally a compound one.

  The occasional longer name is what stops a map reading as a list of
  interchangeable two-syllable tokens."
  [lang k]
  (let [rng (random/rng (:seed lang) (+ 31 (long k)) 3 977)]
    (if (< (rng) 0.22) (compound lang k) (word lang k))))
