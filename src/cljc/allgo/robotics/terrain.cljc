(ns allgo.robotics.terrain
  "Self-supervised learning of what is drivable: Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, section
  26.4.3.

  A laser range finder says reliably whether the ground just ahead is
  flat, but sees only a little way. A camera sees far, but a color is not
  a judgment. So the robot trains the camera on the laser: whatever
  patch the laser finds flat is an example of drivable ground, and a
  mixture of Gaussians, fitted by EM (`allgo.learning.unsupervised`'s),
  learns those patches' colors; every pixel in the image whose color the
  mixture finds likely enough is then drivable too, however far off. The
  robot collects its own labels -- self-supervised -- and as it drives
  onto new ground the laser's new examples teach the camera that too:
  drive from a road onto grass, and grass becomes drivable.

  An image is `{:w :h :px}`, `:px` row by row from the top, each pixel
  `[r g b]`; the bottom rows are nearest."
  (:require [allgo.learning.unsupervised :as unsupervised]
            [clojure.math :as math]))

(defn- log-density
  "The log of the mixture `{:weights :means :vars}`'s density at `x`."
  [{:keys [weights means vars]} x]
  (math/log (max 1e-300 (reduce + (map (fn [w mu var]
                                         (* w (reduce * (map (fn [xi m v] (/ (math/exp (- (/ (* (- xi m) (- xi m)) (* 2 v))))
                                                                             (math/sqrt (* 2 math/PI v))))
                                                             x mu var))))
                                       weights means vars)))))

(defn laser-examples
  "The colors of the pixels of `image` in the `patch`, `[x0 y0 x1 y1]`,
  the small area just ahead that the laser sweeps, that `flat?` (fn [x
  y]) says are flat: the drivable examples it labels."
  [{:keys [w px]} [x0 y0 x1 y1] flat?]
  (vec (for [y (range y0 y1) x (range x0 x1) :when (flat? x y)] (px (+ (* y w) x)))))

(defn learn
  "A model of drivable ground from `examples`, colors: a mixture of `k`
  Gaussians by `steps` iterations of EM, and a threshold on the log
  density, the `quantile` (default 0.02) of the examples' own -- so
  nearly all of what it was taught counts as drivable, and colors far
  from all of it do not."
  [examples k steps rng & {:keys [quantile] :or {quantile 0.02}}]
  (let [mixture (nth (unsupervised/em-gmm-seq examples k rng) steps)
        scores (vec (sort (map #(log-density mixture %) examples)))]
    {:mixture mixture :threshold (scores (long (* quantile (dec (count scores)))))}))

(defn drivable?
  "Whether `model` takes the color `c` for drivable ground."
  [{:keys [mixture threshold]} c]
  (>= (log-density mixture c) threshold))

(defn classify
  "The image as booleans: drivable or not, by color, all the way out."
  [model image]
  (update image :px #(mapv (partial drivable? model) %)))

(defn adapt-seq
  "The robot driving through `frames`, each `{:image :flat?}`: from each
  frame the examples the laser labels in `patch` are added to a window of
  the last `memory` examples, and the model learned again from them. Each
  state `{:model :examples}`."
  [frames patch k steps memory rng]
  (rest (reductions (fn [{:keys [examples]} {:keys [image flat?]}]
                      (let [examples (vec (take-last memory (into examples (laser-examples image patch flat?))))]
                        {:examples examples :model (learn examples k steps rng)}))
                    {:examples [] :model nil}
                    frames)))
