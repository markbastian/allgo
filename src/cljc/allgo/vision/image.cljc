(ns allgo.vision.image
  "Images for computer vision: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, chapter 25.

  An image is `{:w :h :px}`: `:px` a vector of `w` x `h` intensities,
  row by row from the top left -- grayscale, or for a color image a
  vector of `[r g b]` per pixel. Images are values; every operation
  returns a new one.

    convolution   with a kernel separable into a row and a column, so a k
                  by k filter costs 2k a pixel, not k^2; the edges
                  extended outward
    smoothing     by a Gaussian of standard deviation sigma, which the
                  chapter applies before differentiating, since a
                  derivative amplifies the noise a smoothing removes
    gradients     central differences, of a smoothed image the derivative
                  of a Gaussian; the magnitude and the orientation
    sampling      bilinearly between pixels, for the subpixel positions
                  optical flow asks about

  and scenes: images drawn from shapes, with noise, for the tests and
  demos -- whose right answers are known because they were drawn."
  (:require [allgo.numerics.splines :as splines]
            [allgo.random :as random]
            [clojure.math :as math]))

(defn image
  "A `w` by `h` image whose pixel at `(x, y)` is `(f x y)`."
  [w h f]
  {:w w :h h :px (vec (for [y (range h) x (range w)] (f x y)))})

(defn at
  "The pixel at `(x, y)`, the image extended beyond its edges by its edge
  pixels."
  [{:keys [w h px]} x y]
  (px (+ (* (max 0 (min (dec h) y)) w) (max 0 (min (dec w) x)))))

(defn pixels
  "The image with `f` applied to every pixel."
  [img f]
  (update img :px #(mapv f %)))

(defn combine
  "The image whose pixels are `(f a b)` of `a`'s and `b`'s."
  [f a b]
  (assoc a :px (mapv f (:px a) (:px b))))

(defn gray
  "A color image as gray: the mean of its channels."
  [img]
  (pixels img (fn [p] (if (vector? p) (/ (reduce + p) (count p)) p))))

(defn crop
  "The `w` by `h` part of `img` with its top left at `(x0, y0)`."
  [img x0 y0 w h]
  (image w h (fn [x y] (at img (+ x0 x) (+ y0 y)))))

(defn sample
  "The image at the real position `(x, y)`, bilinearly between the four
  pixels round it -- `allgo.numerics.splines`' interpolation."
  [img x y]
  (let [x0 (long (math/floor x)) y0 (long (math/floor y))]
    (splines/bilinear x0 (inc x0) y0 (inc y0)
                      (at img x0 y0) (at img (inc x0) y0) (at img x0 (inc y0)) (at img (inc x0) (inc y0))
                      x y)))

;; ---------------------------------------------------------------------------
;; Filtering

(defn channels
  "A color image as one gray image per channel."
  [{:keys [px] :as img}]
  (mapv (fn [c] (assoc img :px (mapv #(nth % c) px))) (range (count (first px)))))

(defn merge-channels
  "Gray images, one per channel, as a color image."
  [chans]
  (assoc (first chans) :px (apply mapv vector (map :px chans))))

(defn convolve-separable
  "`img` convolved with `row` along x and then `col` along y, each an odd
  number of taps centered on the pixel; a color image a channel at a
  time."
  [img row col]
  (if (vector? (first (:px img)))
    (merge-channels (mapv #(convolve-separable % row col) (channels img)))
    (let [pass (fn [im taps dx dy]
                 (let [r (quot (count taps) 2)]
                   (image (:w im) (:h im)
                          (fn [x y] (reduce + (map-indexed (fn [i t] (* t (at im (+ x (* dx (- i r))) (+ y (* dy (- i r)))))) taps))))))]
      (-> img (pass row 1 0) (pass col 0 1)))))

(defn gaussian-kernel
  "A Gaussian of standard deviation `sigma` sampled at the integers out to
  three sigma, normalized to sum to 1."
  [sigma]
  (let [r (max 1 (long (math/ceil (* 3 sigma))))
        ws (mapv (fn [i] (math/exp (- (/ (* i i) (* 2 sigma sigma))))) (range (- r) (inc r)))
        s (reduce + ws)]
    (mapv #(/ % s) ws)))

(defn smooth
  "`img` smoothed by a Gaussian of standard deviation `sigma`."
  [img sigma]
  (let [k (gaussian-kernel sigma)] (convolve-separable img k k)))

(defn gradients
  "The image's derivatives in x and y by central differences: `[gx gy]`,
  y increasing downward as the rows do."
  [img]
  (let [d [-0.5 0.0 0.5] one [1.0]]
    [(convolve-separable img d one) (convolve-separable img one d)]))

(defn magnitude [[gx gy]] (combine (fn [a b] (math/sqrt (+ (* a a) (* b b)))) gx gy))
(defn orientation [[gx gy]] (combine (fn [a b] (math/atan2 b a)) gx gy))

;; ---------------------------------------------------------------------------
;; Scenes

(defn- inside? [{:keys [kind] :as shape} x y]
  (case kind
    :rect (let [[x0 y0 x1 y1] (:box shape)] (and (<= x0 x) (< x x1) (<= y0 y) (< y y1)))
    :disk (let [[cx cy] (:center shape) r (:radius shape)]
            (<= (+ (* (- x cx) (- x cx)) (* (- y cy) (- y cy))) (* r r)))))

(defn render
  "A `w` by `h` image of `shapes` over `background`, later shapes in front:
  each `{:kind :rect :box [x0 y0 x1 y1]}` or `{:kind :disk :center [x y]
  :radius r}`, with `:value` its intensity (or color) -- or `:texture`,
  a function of the position, for a surface with a pattern on it. Pixels
  are sampled at their centers."
  [w h background shapes]
  (image w h (fn [x y]
               (let [cx (+ x 0.5) cy (+ y 0.5)]
                 (if-let [s (last (filter #(inside? % cx cy) shapes))]
                   (if-let [t (:texture s)] (t cx cy) (:value s))
                   (if (fn? background) (background cx cy) background))))))

(defn add-noise
  "`img` with Gaussian noise of standard deviation `sigma` added to each
  pixel (each channel)."
  [img sigma rng]
  (pixels img (fn [p] (if (vector? p)
                        (mapv #(random/gaussian rng % sigma) p)
                        (random/gaussian rng p sigma)))))

(defn pattern
  "A smooth random pattern with values in about [0, 1]: a few sinusoids
  of random direction and phase summed -- texture enough to match by,
  smooth enough to differentiate."
  [rng n wavelength]
  (let [waves (vec (repeatedly n (fn [] (let [a (* 2 math/PI (rng)) k (/ (* 2 math/PI) (* wavelength (+ 0.5 (rng))))]
                                          [(* k (math/cos a)) (* k (math/sin a)) (* 2 math/PI (rng))]))))]
    (fn [x y] (+ 0.5 (/ (reduce + (map (fn [[kx ky ph]] (math/sin (+ (* kx x) (* ky y) ph))) waves)) (* 2.0 (math/sqrt n)))))))
