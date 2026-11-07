(ns allgo.vision.boundaries
  "Boundary detection by learning: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, section 25.3.4.

  Edges are where the brightness changes; a boundary is where one region
  meets another, and inside a striped region the brightness changes at
  every stripe. So the section's detector looks wider: a disk round each
  pixel, cut in two by a diameter at orientation theta, and the two
  halves' features compared -- here histograms of brightness and of
  gradient orientation (the texture), compared by their chi-squared
  distance. A boundary along theta makes the halves differ; a texture
  that is the same on both sides makes them alike however busy it is.
  The probability Pb(x, y, theta) that a boundary passes is learned from
  images whose boundaries are known: logistic regression
  (`allgo.learning.models`') on each channel's difference, the best
  orientation's.

  Images are `allgo.vision.image`'s, gray."
  (:require [allgo.learning.models :as models]
            [allgo.vision.edges :as edges]
            [allgo.vision.image :as img]
            [clojure.math :as math]))

(defn half-disks
  "The pixel offsets of a disk of radius `r` cut by the diameter at angle
  `theta`: `[one other]`, the offsets on each side of it (those on it
  left out)."
  [r theta]
  (let [c (math/cos theta) s (math/sin theta)
        disk (for [dy (range (- r) (inc r)) dx (range (- r) (inc r))
                   :when (<= (+ (* dx dx) (* dy dy)) (* r r))]
               [dx dy [(- (* c dy) (* s dx))]])
        side (fn [pred] (vec (keep (fn [[dx dy [d]]] (when (pred d) [dx dy])) disk)))]
    [(side #(> % 0.5)) (side #(< % -0.5))]))

(defn chi-squared
  "The chi-squared distance between two histograms: half the sum of (h -
  g)^2 / (h + g), 0 for the same and 1 for two with nothing in common."
  [h g]
  (* 0.5 (reduce + (map (fn [a b] (let [s (+ a b)] (if (pos? s) (/ (* (- a b) (- a b)) s) 0.0))) h g))))

(defn- histogram
  "A normalized histogram of `n` bins from the pairs `[bin weight]`."
  [n pairs]
  (let [h (reduce (fn [h [b w]] (update h b + w)) (vec (repeat n 0.0)) pairs)
        s (reduce + h)]
    (if (pos? s) (mapv #(/ % s) h) h)))

(defn channels
  "What the halves are compared on, per pixel: `{:brightness :texture}`,
  each pixel's `[bin weight]` -- its brightness in one of `bins` bins over
  [0, 1], weight 1, and its gradient orientation (after smoothing by
  `sigma`) in one of `bins` over [0, pi), weighted by the gradient's
  magnitude."
  [image bins sigma]
  (let [g (img/gradients (img/smooth image sigma))
        mag (img/magnitude g)
        orient (img/orientation g)
        bin (fn [v n] (max 0 (min (dec n) (long (* n v)))))]
    {:w (:w image) :h (:h image)
     :brightness (img/pixels image (fn [v] [(bin v bins) 1.0]))
     :texture (img/combine (fn [m a] [(bin (/ (mod a math/PI) math/PI) bins) m]) mag orient)
     :bins bins}))

(defn differences
  "The two halves' chi-squared distances at `(x, y)`, for the disk of
  radius `r` cut at each of `thetas`: `[[brightness texture] ...]`, one
  per orientation."
  [{:keys [brightness texture bins]} x y r thetas]
  (for [theta thetas]
    (let [[one other] (half-disks r theta)
          hist (fn [ch offsets] (histogram bins (map (fn [[dx dy]] (img/at ch (+ x dx) (+ y dy))) offsets)))]
      [(chi-squared (hist brightness one) (hist brightness other))
       (chi-squared (hist texture one) (hist texture other))])))

(defn features
  "The features the classifier reads at `(x, y)`: the two channels'
  differences at the orientation where they sum largest -- the way a
  boundary through the pixel, if there is one, most likely runs."
  [chs x y r thetas]
  (apply max-key #(reduce + %) (differences chs x y r thetas)))

(defn true-boundary
  "The boundary pixels of an image of region `labels`: those with a
  4-neighbor of another region. A set of `[x y]`."
  [{:keys [w h] :as labels}]
  (set (for [y (range h) x (range w)
             :let [l (img/at labels x y)]
             :when (some (fn [[dx dy]] (let [x2 (+ x dx) y2 (+ y dy)]
                                         (and (< -1 x2 w) (< -1 y2 h) (not= l (img/at labels x2 y2)))))
                         [[1 0] [-1 0] [0 1] [0 -1]])]
         [x y])))

(defn- near? [pixels d [x y]]
  (some (fn [[px py]] (and (<= (abs (- px x)) d) (<= (abs (- py y)) d))) pixels))

(defn examples
  "Training examples from images whose boundaries are known: `pairs` of
  `[image labels]`. From each, `n` pixels on a boundary, labeled 1, and
  `n` at least `margin` pixels from any, labeled 0: `[X ys]`."
  [pairs n margin r thetas bins sigma rng]
  (let [rows (for [[image labels] pairs
                   :let [chs (channels image bins sigma)
                         truth (true-boundary labels)
                         {:keys [w h]} image
                         inner (fn [[x y]] (and (<= r x (- w r 1)) (<= r y (- h r 1))))
                         on (vec (filter inner truth))
                         off (vec (for [y (range r (- h r)) x (range r (- w r))
                                        :when (not (near? truth margin [x y]))]
                                    [x y]))
                         pick (fn [ps] (repeatedly n #(ps (long (* (rng) (count ps))))))]
                   [label ps] [[1 (pick on)] [0 (pick off)]]
                   [x y] ps]
               [(features chs x y r thetas) label])]
    [(mapv first rows) (mapv second rows)]))

(defn train-seq
  "The boundary classifier learned from `[X ys]` by logistic regression,
  an epoch of gradient descent a step: the weights."
  [[X ys] alpha]
  (models/logistic-regression-seq X ys alpha))

(defn pb
  "The probability of a boundary at each pixel of `image` under the
  `weights`, Pb(x, y, theta) at its likeliest orientation: `{:p
  :theta}`, images of the probability and of that orientation, zero
  within `r` of the edge, where the disk would leave the picture."
  [weights image r thetas bins sigma]
  (let [chs (channels image bins sigma)
        {:keys [w h]} image
        inner? (fn [x y] (and (<= r x (- w r 1)) (<= r y (- h r 1))))
        best (img/image w h (fn [x y]
                              (if (inner? x y)
                                (apply max-key first (map (fn [d theta] [(models/predict-logistic weights d) theta])
                                                          (differences chs x y r thetas) thetas))
                                [0.0 0.0])))]
    {:p (img/pixels best first) :theta (img/pixels best second)}))

(defn thin
  "The boundary probabilities `pb` gives, kept only where they are a
  maximum across the boundary -- Canny's non-maximum suppression, the
  normal to the boundary standing for the gradient -- so a boundary is a
  pixel wide, not as wide as the disk."
  [{:keys [p theta]}]
  (edges/non-maximum-suppression p (img/pixels theta #(+ % (/ math/PI 2)))))
