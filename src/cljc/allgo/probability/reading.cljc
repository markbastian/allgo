(ns allgo.probability.reading
  "Reading degraded text with a generative program: Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, section
  15.4.

  The program of figure 15.11 makes an image of some letters:

    n             ~ Poisson(lambda)
    letters[i]    ~ uniform over a to z           (figure 15.11)
                    or each given the one before  (figure 15.15, a
                    bigram model estimated from English words)
    noise         ~ UniformReal(0.1, 1)
    pixel         ~ N(the rendered image's pixel, noise)

  Rendering here is a 5 by 7 bitmap font, each letter in a cell six
  pixels wide, ink 1 on 0, optionally blurred. How many letters there
  are is part of the unknown -- an open-universe model, the number of
  random choices depending on one of them.

  Inference is the book's: MCMC over execution traces. A sweep resamples
  each letter given the rest, proposes adding or removing the last letter
  (accepted by the Metropolis-Hastings rule), and resamples the noise on
  a grid. Without blur each letter's pixels are its own, so the
  posterior can also be computed exactly -- for independent letters
  position by position, for the bigram model by the forward algorithm
  -- which is what the MCMC is checked against."
  (:require [allgo.random :as random]
            [clojure.math :as math]
            [clojure.string :as str]))

(def letters (vec "abcdefghijklmnopqrstuvwxyz"))

(def ^:private glyph-rows
  {\a [".....", ".....", ".###.", "....#", ".####", "#...#", ".####"]
   \b ["#....", "#....", "####.", "#...#", "#...#", "#...#", "####."]
   \c [".....", ".....", ".####", "#....", "#....", "#....", ".####"]
   \d ["....#", "....#", ".####", "#...#", "#...#", "#...#", ".####"]
   \e [".....", ".....", ".###.", "#...#", "#####", "#....", ".####"]
   \f ["..##.", ".#...", "####.", ".#...", ".#...", ".#...", ".#..."]
   \g [".....", ".####", "#...#", ".####", "....#", "#...#", ".###."]
   \h ["#....", "#....", "####.", "#...#", "#...#", "#...#", "#...#"]
   \i ["..#..", ".....", ".##..", "..#..", "..#..", "..#..", ".###."]
   \j ["...#.", ".....", "..##.", "...#.", "...#.", "#..#.", ".##.."]
   \k ["#....", "#....", "#..#.", "#.#..", "##...", "#.#..", "#..#."]
   \l [".##..", "..#..", "..#..", "..#..", "..#..", "..#..", ".###."]
   \m [".....", ".....", "##.#.", "#.#.#", "#.#.#", "#.#.#", "#.#.#"]
   \n [".....", ".....", "####.", "#...#", "#...#", "#...#", "#...#"]
   \o [".....", ".....", ".###.", "#...#", "#...#", "#...#", ".###."]
   \p [".....", "####.", "#...#", "#...#", "####.", "#....", "#...."]
   \q [".....", ".####", "#...#", "#...#", ".####", "....#", "....#"]
   \r [".....", ".....", "#.##.", "##...", "#....", "#....", "#...."]
   \s [".....", ".....", ".####", "#....", ".###.", "....#", "####."]
   \t [".#...", ".#...", "####.", ".#...", ".#...", ".#..#", "..##."]
   \u [".....", ".....", "#...#", "#...#", "#...#", "#..##", ".##.#"]
   \v [".....", ".....", "#...#", "#...#", "#...#", ".#.#.", "..#.."]
   \w [".....", ".....", "#...#", "#...#", "#.#.#", "#.#.#", ".#.#."]
   \x [".....", ".....", "#...#", ".#.#.", "..#..", ".#.#.", "#...#"]
   \y [".....", "#...#", "#...#", ".####", "....#", "#...#", ".###."]
   \z [".....", ".....", "#####", "...#.", "..#..", ".#...", "#####"]})

(def glyphs
  "Each letter's 7 by 5 bitmap, rows of 0.0 and 1.0."
  (update-vals glyph-rows (fn [rows] (mapv (fn [r] (mapv #(if (= % \#) 1.0 0.0) r)) rows))))

(def rows 9)
(def ^:private cell 6)

(defn columns
  "The image width for up to `max-letters`: a margin, then a cell each."
  [max-letters]
  (+ 2 (* cell max-letters)))

(defn render
  "The clean image of the letters `text`, `rows` by (columns max): each
  glyph in its cell, a pixel of margin all round; blurred, if `blur?`, by
  averaging each pixel with its eight neighbors."
  [text max-letters blur?]
  (let [w (columns max-letters)
        sharp (reduce (fn [img [i ch]]
                        (reduce (fn [img [r c]] (assoc-in img [(inc r) (+ 1 (* cell i) c)] (get-in glyphs [ch r c])))
                                img
                                (for [r (range 7) c (range 5)] [r c])))
                      (vec (repeat rows (vec (repeat w 0.0))))
                      (map-indexed vector text))]
    (if-not blur?
      sharp
      (vec (for [r (range rows)]
             (vec (for [c (range w)]
                    (/ (reduce + (for [dr [-1 0 1] dc [-1 0 1]] (get-in sharp [(+ r dr) (+ c dc)] 0.0))) 9.0))))))))

(defn degrade
  "`image` with Gaussian noise of variance `noise` added to every pixel."
  [rng image noise]
  (let [sd (math/sqrt noise)]
    (mapv (fn [row] (mapv #(random/gaussian rng % sd) row)) image)))

;; ------------------------------------------------------------- the model

(defn- poisson-log [lambda n]
  (- (* n (math/log lambda)) lambda (reduce + 0.0 (map #(math/log %) (range 1 (inc n))))))

(def ^:private noise-grid (mapv #(+ 0.1 (* 0.05 %)) (range 19)))

(def common-words
  "Common English words, for estimating a bigram model of letters."
  ["the" "of" "and" "to" "in" "is" "it" "you" "that" "he" "was" "for" "on"
   "are" "with" "as" "his" "they" "be" "at" "one" "have" "this" "from" "or"
   "had" "by" "hot" "word" "but" "what" "some" "we" "can" "out" "other"
   "were" "all" "there" "when" "up" "use" "your" "how" "said" "an" "each"
   "she" "which" "do" "their" "time" "if" "will" "way" "about" "many" "then"
   "them" "write" "would" "like" "so" "these" "her" "long" "make" "thing"
   "see" "him" "two" "has" "look" "more" "day" "could" "go" "come" "did"
   "number" "sound" "no" "most" "people" "my" "over" "know" "water" "than"
   "call" "first" "who" "may" "down" "side" "been" "now" "find" "any" "new"
   "work" "part" "take" "get" "place" "made" "live" "where" "after" "back"
   "little" "only" "round" "man" "year" "came" "show" "every" "good" "me"
   "give" "our" "under" "name" "very" "through" "just" "form" "sentence"
   "great" "think" "say" "help" "low" "line" "differ" "turn" "cause" "much"
   "mean" "before" "move" "right" "boy" "old" "too" "same" "tell" "does"
   "set" "three" "want" "air" "well" "also" "play" "small" "end" "put"
   "home" "read" "hand" "port" "large" "spell" "add" "even" "land" "here"
   "must" "big" "high" "such" "follow" "act" "why" "ask" "men" "change"
   "went" "light" "kind" "off" "need" "house" "picture" "try" "us" "again"
   "animal" "point" "mother" "world" "near" "build" "self" "earth" "father"
   "head" "stand" "own" "page" "should" "country" "found" "answer" "school"
   "grow" "study" "still" "learn" "plant" "cover" "food" "sun" "four"
   "between" "state" "keep" "eye" "never" "last" "let" "thought" "city"
   "tree" "cross" "farm" "hard" "start" "might" "story" "saw" "far" "sea"
   "draw" "left" "late" "run" "while" "press" "close" "night" "real" "life"
   "few" "north" "open" "seem" "together" "next" "white" "children" "begin"
   "got" "walk" "example" "ease" "paper" "group" "always" "music" "those"
   "both" "mark" "often" "letter" "until" "mile" "river" "car" "feet" "care"
   "second" "book" "carry" "took" "science" "eat" "room" "friend" "began"
   "idea" "fish" "mountain" "stop" "once" "base" "hear" "horse" "cut" "sure"
   "watch" "color" "face" "wood" "main" "enough" "plain" "girl" "usual"
   "young" "ready" "above" "ever" "red" "list" "though" "feel" "talk" "bird"
   "soon" "body" "dog" "family" "direct" "pose" "leave" "song" "measure"
   "door" "product" "black" "short" "numeral" "class" "wind" "question"
   "happen" "complete" "ship" "area" "half" "rock" "order" "fire" "south"
   "problem" "piece" "told" "knew" "pass" "since" "top" "whole" "king"
   "space" "heard" "best" "hour" "better" "true" "during" "hundred" "five"
   "remember" "step" "early" "hold" "west" "ground" "interest" "reach"
   "fast" "verb" "sing" "listen" "six" "table" "travel" "less" "morning"
   "ten" "simple" "several" "vowel" "toward" "war" "lay" "against" "pattern"
   "slow" "center" "love" "person" "money" "serve" "appear" "road" "map"
   "rain" "rule" "govern" "pull" "cold" "notice" "voice" "unit" "power"
   "town" "fine" "certain" "fly" "fall" "lead" "cry" "dark" "machine" "note"
   "wait" "plan" "figure" "star" "box" "noun" "field" "rest" "correct"
   "able" "pound" "done" "beauty" "drive" "stood" "contain" "front" "teach"
   "week" "final" "gave" "green" "oh" "quick" "develop" "ocean" "warm"
   "free" "minute" "strong" "special" "mind" "behind" "clear" "tail"
   "produce" "fact" "street" "inch" "multiply" "nothing" "course" "stay"
   "wheel" "full" "force" "blue" "object" "decide" "surface" "deep" "moon"
   "island" "foot" "system" "busy" "test" "record" "boat" "common" "gold"
   "possible" "plane" "stead" "dry" "wonder" "laugh" "thousand" "ago" "ran"
   "check" "game" "shape" "equate" "miss" "brought" "heat" "snow" "tire"
   "bring" "yes" "distant" "fill" "east" "paint" "language" "among"])

(defn english-bigrams
  "A bigram model of letters estimated from `words`: `{:initial {letter
  p} :transition {letter {letter p}}}`, every count plus one."
  [words]
  (let [words (map str/lower-case words)
        initial (frequencies (map first words))
        pairs (frequencies (mapcat #(map vector % (rest %)) words))
        norm (fn [counts] (let [c (into {} (map (fn [l] [l (inc (get counts l 0))])) letters)
                                z (reduce + (vals c))]
                            (update-vals c #(/ (double %) z))))]
    {:initial (norm initial)
     :transition (into {} (for [a letters] [a (norm (into {} (for [b letters] [b (get pairs [a b] 0)])))]))}))

(defn model
  "A reading model: `:lambda` the expected number of letters (default
  10, as in the figure), `:max-letters` the image's capacity, `:blur?`,
  and `:bigrams` -- a model from `english-bigrams`, or nil for
  independent letters."
  [& {:keys [lambda max-letters blur? bigrams] :or {lambda 10.0 max-letters 14 blur? false}}]
  {:lambda lambda :max-letters max-letters :blur? blur? :bigrams bigrams})

(defn- text-log-prior
  "log p(text): its length's probability, then each letter's."
  [{:keys [lambda bigrams]} text]
  (+ (poisson-log lambda (count text))
     (if bigrams
       (reduce + (map-indexed (fn [i ch] (math/log (if (zero? i)
                                                     (get-in bigrams [:initial ch])
                                                     (get-in bigrams [:transition (nth text (dec i)) ch]))))
                              text))
       (* (count text) (math/log (/ 1.0 26))))))

(defn log-likelihood
  "log p(image | text, noise)."
  [{:keys [max-letters blur?]} image text noise]
  (let [clean (render text max-letters blur?)
        k (* -0.5 (math/log (* 2.0 math/PI noise)))]
    (reduce + (map (fn [orow crow]
                     (reduce + (map (fn [o c] (- k (/ (* (- o c) (- o c)) (* 2.0 noise)))) orow crow)))
                   image clean))))

(defn log-joint
  "log p(text, noise, image); the noise's prior is uniform, a constant."
  [m image text noise]
  (+ (text-log-prior m text) (log-likelihood m image text noise)))

;; ------------------------------------------------------------------- MCMC

(defn- log-sum-exp [xs]
  (let [top (reduce max ##-Inf xs)]
    (+ top (math/log (reduce + (map #(math/exp (- % top)) xs))))))

(defn- draw [rng log-weights]
  (let [z (log-sum-exp log-weights)]
    (random/weighted rng (map-indexed (fn [i w] [i (math/exp (- w z))]) log-weights))))

(defn- resample-letter
  "Letter `i` drawn from its conditional given the rest of the trace."
  [m image rng {:keys [text noise] :as trace} i]
  (let [options (mapv #(str (subs text 0 i) % (subs text (inc i))) letters)
        k (draw rng (mapv #(log-joint m image % noise) options))]
    (assoc trace :text (options k))))

(defn- birth-death
  "Propose adding a letter at the end or removing the last, each half
  the time (always adding to an empty text, always removing from a full
  one), accepted by the Metropolis-Hastings rule. The added letter is
  drawn from its conditional, the letters weighted by the joint
  probability with each, so the move is data-driven and its acceptance
  needs only their sum: for a birth, sum over l of p(text + l) over
  p(text), times the reverse move's chance over the forward's; a death
  the reverse."
  [{:keys [max-letters] :as m} image rng {:keys [text noise] :as trace}]
  (let [n (count text)
        p-birth (fn [k] (cond (zero? k) 1.0 (>= k max-letters) 0.0 :else 0.5))
        births (fn [t] (log-sum-exp (mapv #(log-joint m image (str t %) noise) letters)))
        grow? (< (rng) (p-birth n))]
    (if grow?
      (let [options (mapv #(str text %) letters)
            lws (mapv #(log-joint m image % noise) options)
            a (+ (- (log-sum-exp lws) (log-joint m image text noise))
                 (math/log (/ (- 1.0 (p-birth (inc n))) (p-birth n))))]
        (if (< (math/log (rng)) a) (assoc trace :text (options (draw rng lws))) trace))
      (let [shorter (subs text 0 (dec n))
            a (+ (- (log-joint m image shorter noise) (births shorter))
                 (math/log (/ (p-birth (dec n)) (- 1.0 (p-birth n)))))]
        (if (< (math/log (rng)) a) (assoc trace :text shorter) trace)))))

(defn- resample-noise
  "The noise drawn from its conditional on the grid."
  [m image rng {:keys [text] :as trace}]
  (assoc trace :noise (noise-grid (draw rng (mapv #(log-likelihood m image text %) noise-grid)))))

(defn- sweep
  "One MCMC sweep over the trace: each letter resampled, a few proposals
  to lengthen or shorten it, the noise resampled."
  [m image rng trace]
  (as-> trace t
    (reduce (partial resample-letter m image rng) t (range (count (:text t))))
    (reduce (fn [t _] (birth-death m image rng t)) t (range 4))
    (resample-noise m image rng t)))

(defn mcmc-seq
  "Traces `{:text :noise}` sampled by MCMC sweeps from the empty one."
  [m image rng]
  (iterate (partial sweep m image rng) {:text "" :noise 0.5}))

;; ------------------------------------------------------------ exactly

(defn- block-log-likelihood
  "log p(the pixels of cell `i` | the cell holds `ch`, or is blank if
  nil), unblurred: the cell's six columns."
  [image i ch noise]
  (let [k (* -0.5 (math/log (* 2.0 math/PI noise)))]
    (reduce + (for [r (range rows) c (range cell)
                    :let [x (+ 1 (* cell i) c)
                          clean (if (and ch (< 0 r 8) (< c 5)) (get-in glyphs [ch (dec r) c]) 0.0)
                          o (get-in image [r x])]]
                (- k (/ (* (- o clean) (- o clean)) (* 2.0 noise)))))))

(defn- margins-log-likelihood
  "log p(the pixels outside every cell, all blank)."
  [image noise]
  (let [k (* -0.5 (math/log (* 2.0 math/PI noise)))
        w (count (first image))]
    (reduce + (for [r (range rows) c [0 (dec w)]
                    :let [o (get-in image [r c])]]
                (- k (/ (* o o) (* 2.0 noise)))))))

(defn- forward
  "log p(cells 0..n-1 hold letters | noise) summed over the letters by
  the forward algorithm, and the forward messages, each a vector over
  `letters` of log probabilities."
  [{:keys [bigrams]} blocks n]
  (let [prior0 (if bigrams (mapv #(math/log (get-in bigrams [:initial %])) letters) (vec (repeat 26 (math/log (/ 1.0 26)))))
        trans (fn [a b] (if bigrams (math/log (get-in bigrams [:transition a b])) (math/log (/ 1.0 26))))
        alphas (reductions (fn [alpha i]
                             (mapv (fn [b] (+ (get-in blocks [i b])
                                              (log-sum-exp (map (fn [a] (+ (alpha a) (trans (letters a) (letters b)))) (range 26)))))
                                   (range 26)))
                           (mapv + prior0 (get blocks 0))
                           (range 1 n))]
    (vec (take n alphas))))

(defn- viterbi
  "The most probable letters for cells 0..n-1 given their block log
  likelihoods, by the Viterbi algorithm over the letter model."
  [{:keys [bigrams]} blocks n]
  (let [prior0 (if bigrams (mapv #(math/log (get-in bigrams [:initial %])) letters) (vec (repeat 26 (math/log (/ 1.0 26)))))
        trans (fn [a b] (if bigrams (math/log (get-in bigrams [:transition (letters a) (letters b)])) (math/log (/ 1.0 26))))
        step (fn [[delta back] i]
               (let [best (mapv (fn [b] (apply max-key #(+ (delta %) (trans % b)) (range 26))) (range 26))]
                 [(mapv (fn [b] (+ (get-in blocks [i b]) (delta (best b)) (trans (best b) b))) (range 26))
                  (conj back best)]))
        [delta back] (reduce step [(mapv + prior0 (get blocks 0)) []] (range 1 n))
        last-letter (apply max-key delta (range 26))]
    (apply str (map letters (reduce (fn [path best] (cons (best (first path)) path)) (list last-letter) (rseq back))))))

(defn exact-posterior
  "The exact posterior of an unblurred image: `{:length {n p}}`; for the
  most probable length and noise together, the most probable text `:text`
  (by the Viterbi algorithm); and, for independent letters, each
  position's letter probabilities at the most probable length,
  `:letters`."
  [{:keys [lambda max-letters bigrams] :as m} image]
  (let [per-noise
        (for [noise noise-grid
              :let [blocks (mapv (fn [i] (mapv #(block-log-likelihood image i % noise) letters)) (range max-letters))
                    blanks (mapv (fn [i] (block-log-likelihood image i nil noise)) (range max-letters))
                    margin (margins-log-likelihood image noise)
                    alphas (forward m blocks max-letters)]]
          {:noise noise :blocks blocks :blanks blanks
           :by-n (mapv (fn [n] (+ (poisson-log lambda n) margin
                                  (reduce + (subvec blanks n))
                                  (if (zero? n) 0.0 (log-sum-exp (alphas (dec n))))))
                       (range (inc max-letters)))})
        joint (for [{:keys [by-n] :as pn} per-noise n (range (inc max-letters))] [n (by-n n) pn])
        z (log-sum-exp (map second joint))
        length (reduce (fn [acc [n lw]] (update acc n (fnil + 0.0) (math/exp (- lw z)))) {} joint)
        best-n (key (apply max-key val length))
        ;; letter marginals at the best length, the noise summed out; with
        ;; independent letters each position is its own
        marginals (when (and (pos? best-n) (nil? bigrams))
                    (let [weights (map (fn [{:keys [by-n]}] (math/exp (- (by-n best-n) z))) per-noise)
                          wz (reduce + weights)]
                      (vec (for [i (range best-n)]
                             (let [per (map (fn [{:keys [blocks]} w]
                                              (let [lp (blocks i) s (log-sum-exp lp)]
                                                (mapv #(* (/ w wz) (math/exp (- % s))) lp)))
                                            per-noise weights)]
                               (zipmap letters (apply mapv + per)))))))
        [n _ {:keys [blocks]}] (apply max-key second joint)]
    {:length length
     :letters marginals
     :text (if (zero? n) "" (viterbi m blocks n))}))
