(ns allgo.vision.edges
  "Edges and texture: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 25, which teaches them in the
  text.

  An edge is where the brightness changes sharply. Differentiating finds
  that, and finds the noise too, so the image is smoothed first; the
  edge is then where the gradient's magnitude is greatest across it --
  thinner than where it is merely large -- and edges are kept where they
  are strong, or weak but joined to strong ones. That is Canny's
  detector, a stage a step:

    smoothed    the image convolved with a Gaussian
    gradient    its derivatives, their magnitude
    thinned     non-maximum suppression: a pixel kept only if its
                magnitude is at least its two neighbors' along the
                gradient
    hysteresis  pixels over `high` start edges, and edges grow into
                touching pixels over `low`, a ring a step

  and then the edge pixels linked into curves: neighbors whose
  orientations agree belong to the same one.

  Texture is described by what edges there are, regardless of where:
  `orientation-histograms` counts, in each cell of a grid, the gradient
  orientations weighted by their magnitude -- the descriptor the chapter
  builds image features from."
  (:require [allgo.vision.image :as img]
            [clojure.math :as math]))

(defn- neighbors8 [[x y]]
  (for [dx [-1 0 1] dy [-1 0 1] :when (not= 0 dx dy)] [(+ x dx) (+ y dy)]))

(defn non-maximum-suppression
  "The magnitude where it is a maximum across the edge -- at least its two
  neighbors along the gradient direction, quantized to the nearest of
  four -- and zero elsewhere."
  [mag orient]
  (img/image (:w mag) (:h mag)
             (fn [x y]
               (let [m (img/at mag x y)
                     a (img/at orient x y)
                     ;; the gradient direction, folded to [0, pi) and
                     ;; rounded to 0, 45, 90 or 135 degrees
                     k (mod (long (math/round (/ (mod a math/PI) (/ math/PI 4)))) 4)
                     [dx dy] ([[1 0] [1 1] [0 1] [-1 1]] k)]
                 (if (and (>= m (img/at mag (+ x dx) (+ y dy))) (>= m (img/at mag (- x dx) (- y dy))))
                   m
                   0.0)))))

(defn hysteresis-seq
  "Edge growth from the pixels of `thinned` above `high` into touching
  pixels above `low`, a ring a step: each state `{:edges :frontier}`,
  ending when nothing more joins."
  [thinned low high]
  (let [{:keys [w h]} thinned
        strong (set (for [y (range h) x (range w) :when (> (img/at thinned x y) high)] [x y]))
        weak? (fn [[x y]] (and (< -1 x w) (< -1 y h) (> (img/at thinned x y) low)))]
    (take-while some?
                (iterate (fn [{:keys [edges frontier]}]
                           (when (seq frontier)
                             (let [new (set (for [p frontier q (neighbors8 p) :when (and (not (edges q)) (weak? q))] q))]
                               {:edges (into edges new) :frontier new})))
                         {:edges strong :frontier strong}))))

(defn canny-seq
  "Canny's edge detector on `image` a stage at a time: `{:stage :image}`
  for the smoothed image, the gradient magnitude and the thinned one,
  then `{:stage :hysteresis :edges}` as the edges grow. Options `:sigma`
  (default 1.0), and `:low` and `:high`, thresholds on the magnitude
  (default 0.04 and 0.1)."
  [image & {:keys [sigma low high] :or {sigma 1.0 low 0.04 high 0.1}}]
  (let [smoothed (img/smooth image sigma)
        g (img/gradients smoothed)
        mag (img/magnitude g)
        thinned (non-maximum-suppression mag (img/orientation g))]
    (concat [{:stage :smoothed :image smoothed}
             {:stage :gradient :image mag}
             {:stage :thinned :image thinned}]
            (map #(assoc % :stage :hysteresis) (hysteresis-seq thinned low high)))))

(defn canny
  "The edge pixels Canny's detector finds, a set of `[x y]`."
  [image & opts]
  (:edges (last (apply canny-seq image opts))))

(defn orientation-histograms
  "The image's texture: for each `cell` by `cell` square, a histogram of
  its gradient orientations in `bins` bins over [0, pi) -- an edge's
  orientation, regardless of which side is brighter -- each vote
  weighted by its gradient's magnitude, and the histogram normalized to
  sum to 1 (left at zero where there is no gradient). `{[cx cy] [h ...]}`."
  [image cell bins]
  (let [g (img/gradients image)
        mag (img/magnitude g)
        orient (img/orientation g)]
    (into {} (for [cy (range (quot (:h image) cell)) cx (range (quot (:w image) cell))]
               (let [hist (reduce (fn [h [x y]]
                                    (let [a (mod (img/at orient x y) math/PI)
                                          b (min (dec bins) (long (* bins (/ a math/PI))))]
                                      (update h b + (img/at mag x y))))
                                  (vec (repeat bins 0.0))
                                  (for [y (range (* cy cell) (* (inc cy) cell)) x (range (* cx cell) (* (inc cx) cell))] [x y]))
                     s (reduce + hist)]
                 [[cx cy] (if (pos? s) (mapv #(/ % s) hist) hist)])))))

(defn- turn
  "The angle between two edge orientations, which are lines, not
  arrows: folded into [0, pi / 2]."
  [a b]
  (let [d (mod (- a b) math/PI)] (min d (- math/PI d))))

(defn link-seq
  "Edge pixels linked into curves, as the chapter links them: two
  neighboring edge pixels belong to the same curve when their
  orientations (from `orient`, an image of gradient angles) differ by no
  more than `max-turn` radians. A curve grows a ring of neighbors a step
  from an unlinked pixel, and the next curve starts when it stops: each
  state `{:curves :curve :frontier :left}`."
  [edges orient max-turn]
  (let [angle (fn [[x y]] (img/at orient x y))
        start (fn [curves left]
                (if-let [p (first left)]
                  {:curves curves :curve #{p} :frontier #{p} :left (disj left p)}
                  {:curves curves :curve nil :frontier nil :left left}))]
    (take-while some?
                (iterate (fn [{:keys [curves curve frontier left]}]
                           (cond
                             (seq frontier)
                             (let [new (set (for [p frontier q (neighbors8 p)
                                                  :when (and (left q) (<= (turn (angle p) (angle q)) max-turn))]
                                              q))]
                               {:curves curves :curve (into curve new) :frontier new :left (reduce disj left new)})
                             curve (start (conj curves curve) left)
                             :else nil))
                         (start [] (into (sorted-set) edges))))))

(defn link
  "The curves `link-seq` finds, largest first: a vector of sets of `[x
  y]`."
  [edges orient max-turn]
  (vec (sort-by (comp - count) (:curves (last (link-seq edges orient max-turn))))))
