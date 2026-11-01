(ns allgo.robotics.localization
  "Monte Carlo localization: Russell and Norvig, *Artificial Intelligence:
  A Modern Approach*, 4th edition, chapter 26 (MONTE-CARLO-LOCALIZATION,
  figure 26.6), as steps.

  A robot in a known map of walls, with a range scanner and noisy
  wheels, keeps `N` guesses of where it is -- particles, each a pose `[x
  y heading]`. Each step every particle is moved by the motion model
  with its own noise, weighted by how well the scan it would see from
  there -- found by casting rays against the map -- matches the scan the
  robot did see, and the set resampled in proportion to the weights. The
  particles that survive cluster where the robot is.

  A map is a vector of wall segments `[[x0 y0] [x1 y1]]`; `scan-angles`
  are the beams, relative to the heading."
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.random :as random]
            [clojure.math :as math]))

(defn raycast
  "RAYCAST: the distance along the ray from `[x y]` at `angle` to the
  nearest wall of `walls`, or `max-range`."
  [walls [x y] angle max-range]
  (let [end [(+ x (* max-range (math/cos angle))) (+ y (* max-range (math/sin angle)))]]
    (reduce (fn [best [a b]]
              (if-let [[t u] (polygon/segment-intersection [x y] end a b)]
                (if (and (<= 0.0 t 1.0) (<= 0.0 u 1.0)) (min best (* t max-range)) best)
                best))
            max-range walls)))

(defn scan
  "The ranges the scanner sees from `pose`."
  [walls [x y heading] angles max-range]
  (mapv #(raycast walls [x y] (+ heading %) max-range) angles))

(defn move
  "The motion model: `pose` driven at speed `v` and turn rate `w` for
  `dt`, with Gaussian noise of `sigma-v` and `sigma-w` on each."
  [[x y heading] v w dt rng sigma-v sigma-w]
  (let [v (random/gaussian rng v sigma-v)
        w (random/gaussian rng w sigma-w)
        heading' (+ heading (* w dt))]
    [(+ x (* v dt (math/cos heading'))) (+ y (* v dt (math/sin heading'))) heading']))

(defn- inside? [[[x0 y0] [x1 y1]] [x y]] (and (< x0 x x1) (< y0 y y1)))

(defn- likelihood
  "P(z | z*): each beam Gaussian about the range expected, independently."
  [z z* sigma]
  (reduce * (map (fn [a b] (let [d (- a b)] (math/exp (- (/ (* d d) (* 2 sigma sigma)))))) z z*)))

(defn mcl-step
  "One MONTE-CARLO-LOCALIZATION update: the particles moved by `[v w]`,
  weighted by the scan `z`, and resampled."
  [{:keys [walls angles max-range dt sigma-v sigma-w sigma-z rng bounds random-fraction]
    :or {random-fraction 0.0}} particles [v w] z]
  (let [moved (mapv #(move % v w dt rng sigma-v sigma-w) particles)
        weights (mapv (fn [p] (if (inside? bounds p) (likelihood z (scan walls p angles max-range) sigma-z) 0.0)) moved)
        total (reduce + weights)]
    (if (pos? total)
      (let [cumulative (vec (reductions + (map #(/ % total) weights)))
            n (count moved)
            ;; Low-variance resampling: one random offset, n evenly spaced
            ;; pointers into the cumulative weights.
            r (/ (double (rng)) n)]
        (loop [i 0 j 0 out (transient [])]
          (if (= i n)
            ;; A few particles replaced by random ones: the remedy for
            ;; sample impoverishment, which otherwise leaves every particle a
            ;; copy of one early guess -- right or not.
            (let [[[x0 y0] [x1 y1]] bounds]
              (mapv (fn [p] (if (< (double (rng)) random-fraction)
                              [(+ x0 (* (- x1 x0) (rng))) (+ y0 (* (- y1 y0) (rng))) (* 2 math/PI (rng))]
                              p))
                    (persistent! out)))
            (let [u (+ r (/ (double i) n))
                  j (loop [j j] (if (and (< j (dec n)) (> u (cumulative j))) (recur (inc j)) j))]
              (recur (inc i) j (conj! out (moved j)))))))
      ;; Nothing fits: start again from anywhere.
      (let [[[x0 y0] [x1 y1]] bounds]
        (vec (repeatedly (count particles) #(vector (+ x0 (* (- x1 x0) (rng))) (+ y0 (* (- y1 y0) (rng))) (* 2 math/PI (rng)))))))))

(defn mcl-seq
  "Monte Carlo localization of a robot that starts at `start` and is driven
  by `controls`, each `[v w]`: each state `{:truth :particles :scan}`, the
  robot's true pose (moved by the same noisy motion model), the
  particles, and the scan it saw. `n` particles start spread over the
  map's `bounds` at random -- the robot does not know where it is."
  [{:keys [walls angles max-range dt sigma-v sigma-w sigma-z rng bounds] :as world} start controls n]
  (let [[[x0 y0] [x1 y1]] bounds
        particles (vec (repeatedly n #(vector (+ x0 (* (- x1 x0) (rng))) (+ y0 (* (- y1 y0) (rng))) (* 2 math/PI (rng)))))]
    (reductions (fn [{:keys [truth particles]} [v w]]
                  (let [truth' (move truth v w dt rng sigma-v sigma-w)
                        ;; The scanner's own noise.
                        z (mapv #(random/gaussian rng % sigma-z) (scan walls truth' angles max-range))]
                    {:truth truth' :scan z :particles (mcl-step world particles [v w] z)}))
                {:truth start :particles particles :scan nil}
                controls)))

(defn estimate
  "Where the particles say the robot is: the mean of those within a meter
  of the densest spot -- not of them all, since while the evidence is
  ambiguous they form more than one cluster, and the mean of two clusters
  is where neither is -- and the spread about it: `[[x y] spread]`."
  [particles]
  (let [cell (fn [[x y]] [(long (math/floor (* 2 x))) (long (math/floor (* 2 y)))])
        densest (key (apply max-key val (frequencies (map cell particles))))
        [cx cy] (map #(/ (+ % 0.5) 2.0) densest)
        near (filter (fn [[x y]] (< (+ (* (- x cx) (- x cx)) (* (- y cy) (- y cy))) 1.0)) particles)
        n (count near)
        mx (/ (reduce + (map first near)) n)
        my (/ (reduce + (map second near)) n)
        particles near]
    [[mx my] (math/sqrt (/ (reduce + (map (fn [[x y]] (+ (* (- x mx) (- x mx)) (* (- y my) (- y my)))) particles)) n))]))

(defn room
  "A test map: a 10 by 6 room with an L-shaped wall and a pillar inside,
  so that few places look alike."
  []
  (let [box (fn [[x0 y0] [x1 y1]] [[[x0 y0] [x1 y0]] [[x1 y0] [x1 y1]] [[x1 y1] [x0 y1]] [[x0 y1] [x0 y0]]])]
    {:walls (vec (concat (box [0.0 0.0] [10.0 6.0])
                         [[[3.0 0.0] [3.0 2.5]] [[3.0 2.5] [5.0 2.5]]]
                         (box [7.0 3.5] [8.0 4.5])))
     :bounds [[0.0 0.0] [10.0 6.0]]
     :angles (mapv #(* % (/ math/PI 4)) (range 8))
     :max-range 12.0 :dt 0.1 :sigma-v 0.3 :sigma-w 0.2 :sigma-z 0.6 :random-fraction 0.02}))
