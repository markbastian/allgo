(ns allgo.probability.association
  "Data association: tracking several objects when no observation says
  which object made it. Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, section 15.3.1.

  Each object moves by a linear-Gaussian model and is seen through one,
  as in a Kalman filter (`allgo.probability.temporal`), and each scan is
  a set of blips with no labels. Some blips may be false alarms
  (clutter), and an object may go unseen (a detection failure). For n
  objects and T scans there are (n!)^T ways to explain even the
  clutter-free blips, and no efficient exact filter. The section gives
  the approximations:

    nearest neighbor   each scan, repeatedly pair the closest predicted
                       position and blip; good when the objects are far
                       apart and the noise small
    most probable      each scan, the assignment of blips to objects that
                       maximizes their joint likelihood, by the
                       Hungarian method (`allgo.graph/min-cost-assignment`)
    particles          many association histories at once, each particle
                       drawing its scan's explanation in proportion to
                       how likely it is and weighted by how likely the
                       scan was, then resampled. Given its history, each
                       object's track is computed exactly by a Kalman
                       filter: the Rao-Blackwellization the section
                       recommends. A particle can recover from a wrong
                       guess the single-assignment methods commit to.

  `exact-posterior` enumerates every history, for checking the others on
  a few scans.

  A track is a Kalman belief `[mean cov]`. A model is `allgo.probability.
  temporal`'s `{:F :Q :H :R}` with `:p-detect`, the chance an object is
  seen, and `:clutter`, the density of false alarms over the plane (their
  expected number per unit area per scan)."
  (:require [allgo.graph :as graph]
            [allgo.numerics.linear :as lin]
            [allgo.probability.gaussian :as gaussian]
            [allgo.probability.temporal :as temporal]
            [allgo.random :as random]
            [clojure.math :as math]))

(defn constant-velocity
  "A model of objects moving at nearly constant velocity in the plane,
  state [x y vx vy], one time unit a scan, observed in position: `q` the
  acceleration noise's spectral density, `r` the observation noise's
  standard deviation."
  [q r]
  {:F [[1.0 0.0 1.0 0.0] [0.0 1.0 0.0 1.0] [0.0 0.0 1.0 0.0] [0.0 0.0 0.0 1.0]]
   :Q (lin/mat-scale [[(/ 1.0 3) 0.0 0.5 0.0] [0.0 (/ 1.0 3) 0.0 0.5] [0.5 0.0 1.0 0.0] [0.0 0.5 0.0 1.0]] q)
   :H [[1.0 0.0 0.0 0.0] [0.0 1.0 0.0 0.0]]
   :R [[(* r r) 0.0] [0.0 (* r r)]]
   :p-detect 1.0
   :clutter 0.0})

(defn- log-likelihood
  "log p(z | the forecast track)."
  [model track z]
  (let [[expected S] (temporal/kalman-observation model track)]
    (gaussian/log-density expected S z)))

(defn position
  "The mean position of `track`."
  [model [mean]]
  (lin/mat-vec (:H model) mean))

;; ------------------------------------------------------- one assignment

(defn- nearest-step
  "One scan by the nearest-neighbor filter: forecast every track, then
  repeatedly pair the closest remaining track and blip."
  [model tracks blips]
  (let [forecast (mapv #(temporal/kalman-forecast model %) tracks)
        pairs (sort-by first (for [i (range (count forecast)) k (range (count blips))]
                               [(lin/distance (position model (forecast i)) (blips k)) i k]))
        chosen (reduce (fn [{:keys [tracks blips] :as acc} [_ i k]]
                         (if (or (tracks i) (blips k)) acc (-> acc (assoc-in [:tracks i] k) (assoc-in [:blips k] i))))
                       {:tracks {} :blips {}} pairs)]
    (mapv (fn [i] (if-let [k (get-in chosen [:tracks i])]
                    (temporal/kalman-correct model (forecast i) (blips k))
                    (forecast i)))
          (range (count forecast)))))

(defn- probable-step
  "One scan by the most probable assignment: forecast every track, then
  the explanation of the scan of greatest probability, by the Hungarian
  method. A track's cost for a blip is minus the log of p-detect times
  the blip's likelihood, over the clutter density the blip would
  otherwise be put down to; each track also has a column of its own for
  going unseen, costing minus the log of 1 - p-detect."
  [{:keys [p-detect clutter] :as model} tracks blips]
  (let [forecast (mapv #(temporal/kalman-forecast model %) tracks)
        n (count forecast) m (count blips)
        never 1e12
        cost (vec (for [i (range n)]
                    (vec (concat (for [z blips]
                                   (- (if (pos? clutter) (math/log clutter) 0.0)
                                      (math/log p-detect) (log-likelihood model (forecast i) z)))
                                 (for [j (range n)]
                                   (if (and (= i j) (< p-detect 1.0)) (- (math/log (- 1.0 p-detect))) never))))))
        cols (graph/min-cost-assignment cost)]
    (mapv (fn [t k] (if (< k m) (temporal/kalman-correct model t (blips k)) t)) forecast cols)))

(defn nearest-neighbor-seq
  "The tracks after each scan by the nearest-neighbor filter, from the
  prior `tracks`, each scan a vector of blips. A track with no blip left
  for it keeps its forecast."
  [model tracks scans]
  (rest (reductions (partial nearest-step model) tracks scans)))

(defn most-probable-seq
  "The tracks after each scan by the most probable assignment."
  [model tracks scans]
  (rest (reductions (partial probable-step model) tracks scans)))

;; ------------------------------------------------- every explanation

(defn- explanations
  "Every way to explain a scan of `m` blips by `n` objects: for each
  object, the blip it made or nil if unseen, no blip twice; the blips
  left over are clutter."
  [n m]
  (reduce (fn [partials _]
            (for [p partials k (cons nil (remove (set p) (range m)))] (conj p k)))
          [[]]
          (range n)))

(defn- explanation-log-weight
  "log p(scan, explanation | forecast tracks): each seen object's blip
  likelihood times p-detect, each unseen object's 1 - p-detect, each
  false alarm the clutter density."
  [{:keys [p-detect clutter] :as model} forecast blips explanation]
  (let [seen (count (remove nil? explanation))]
    (+ (reduce + (map (fn [t k] (if k (+ (math/log p-detect) (log-likelihood model t (blips k))) (math/log (- 1.0 p-detect))))
                      forecast explanation))
       (let [false-alarms (- (count blips) seen)]
         (cond (zero? false-alarms) 0.0
               (pos? clutter) (* false-alarms (math/log clutter))
               :else ##-Inf)))))

(defn- log-sum-exp [xs]
  (let [m (reduce max ##-Inf xs)]
    (if (= m ##-Inf) m (+ m (math/log (reduce + (map #(math/exp (- % m)) xs)))))))

(defn- corrected [model forecast blips explanation]
  (mapv (fn [t k] (if k (temporal/kalman-correct model t (blips k)) t)) forecast explanation))

(defn- particle-step
  "One scan for one particle: its tracks forecast, its scan's explanation
  drawn in proportion to its weight, the tracks corrected by it, and the
  particle's log weight increased by the log of the scan's likelihood,
  the explanations' weights summed."
  [model rng {:keys [tracks history log-w]} blips]
  (let [forecast (mapv #(temporal/kalman-forecast model %) tracks)
        options (vec (explanations (count tracks) (count blips)))
        lws (mapv #(explanation-log-weight model forecast blips %) options)
        total (log-sum-exp lws)
        pick (random/weighted rng (zipmap (range (count options)) (map #(math/exp (- % total)) lws)))
        e (options pick)]
    {:tracks (corrected model forecast blips e) :history (conj history e) :log-w (+ log-w total)}))

(defn- resample
  "Particles drawn with replacement in proportion to their weights, the
  weights then made equal."
  [rng particles]
  (let [total (log-sum-exp (map :log-w particles))
        dist (frequencies (repeatedly (count particles)
                                      #(random/weighted rng (zipmap (range (count particles))
                                                                    (map (fn [p] (math/exp (- (:log-w p) total))) particles)))))]
    (vec (for [[i k] dist _ (range k)] (assoc (particles i) :log-w 0.0)))))

(defn- particle-scan [model rng particles blips]
  (resample rng (mapv #(particle-step model rng % blips) particles)))

(defn particle-seq
  "The particles after each scan, `n` of them from the prior `tracks`:
  each `{:tracks :history :log-w}`, the history its explanation of every
  scan so far -- for each object, the index of its blip or nil."
  [model tracks scans rng n]
  (rest (reductions (partial particle-scan model rng)
                    (vec (repeat n {:tracks tracks :history [] :log-w 0.0}))
                    scans)))

(defn history-posterior
  "The particles' estimate of the posterior over association histories:
  how often each history appears."
  [particles]
  (update-vals (frequencies (map :history particles)) #(/ (double %) (count particles))))

(defn- exact-step
  "Every history so far extended by every explanation of the next scan,
  each with its tracks and log weight."
  [model hypotheses blips]
  (for [{:keys [tracks history log-w]} hypotheses
        :let [forecast (mapv #(temporal/kalman-forecast model %) tracks)]
        e (explanations (count tracks) (count blips))]
    {:tracks (corrected model forecast blips e)
     :history (conj history e)
     :log-w (+ log-w (explanation-log-weight model forecast blips e))}))

(defn exact-posterior
  "The posterior probability of every association history, by
  enumerating them all -- exponentially many, so for a few scans."
  [model tracks scans]
  (let [final (reduce (partial exact-step model) [{:tracks tracks :history [] :log-w 0.0}] scans)
        total (log-sum-exp (map :log-w final))]
    (into {} (map (fn [{:keys [history log-w]}] [history (math/exp (- log-w total))])) final)))
