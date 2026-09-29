(ns allgo.astro.close-approach
  "When two orbiting objects pass closest, and how close: the method of
  Alfano and Negron, \"Determining satellite close approaches\", Journal
  of the Astronautical Sciences 41, 1993 (ANCAS; Vallado, algorithm 75).

  A search by sampling the distance would miss a pass between samples, and
  close passes at orbital speeds last seconds. ANCAS samples the two
  ephemerides coarsely -- minutes apart -- and over each interval fits
  their relative position with the cubic that matches the positions and
  velocities at both ends. The range-rate function f = rho . rho', whose
  zeros going from negative to positive are the moments of closest
  approach, is then a quintic in the interval's time, and its roots there
  can all be found. Each is reported, and refined on the ephemerides
  themselves.

  Ephemerides are functions of time (s) giving `[r v]`, km and km/s."
  (:require [allgo.geometry.vec3 :as v3]
            [allgo.numerics.roots :as roots]))

(defn- hermite
  "The cubic in tau in [0, 1] through `p0` `p1` with rates `d0` `d1` (per
  unit tau): its coefficients per component, constant term first."
  [p0 d0 p1 d1]
  (mapv (fn [a b c d]
          [a b (- (* 3.0 (- c a)) (* 2.0 b) d) (+ (* 2.0 (- a c)) b d)])
        p0 d0 p1 d1))

(defn- at [cs tau] (mapv (fn [[a b c d]] (+ a (* tau (+ b (* tau (+ c (* tau d))))))) cs))
(defn- rate [cs tau] (mapv (fn [[_ b c d]] (+ b (* tau (+ (* 2.0 c) (* 3.0 tau d))))) cs))

(defn- interval-minima
  "The taus in [0, 1) where the fitted range-rate turns from negative to
  positive: sampled finely enough to separate a quintic's roots, then
  bisected."
  [cs]
  (let [f #(v3/dot (at cs %) (rate cs %))
        taus (map #(/ % 64.0) (range 65))]
    (for [[a b] (partition 2 1 taus)
          :when (and (neg? (f a)) (not (neg? (f b))))]
      (roots/bisect f a b))))

(defn- refine
  "The time of least distance on the true ephemerides near `t`, by
  Newton's method on the range-rate, its derivative by central
  difference."
  [eph-a eph-b t h]
  (let [g (fn [t] (let [[ra va] (eph-a t) [rb vb] (eph-b t)]
                    (v3/dot (v3/sub ra rb) (v3/sub va vb))))
        dt (* 1e-3 h)]
    (loop [t t i 0]
      (let [gt (g t)
            dg (/ (- (g (+ t dt)) (g (- t dt))) (* 2.0 dt))
            step (/ gt dg)]
        (if (or (> i 30) (< (abs step) 1e-9) (not (pos? dg)))
          t
          (recur (- t step) (inc i)))))))

(defn close-approaches
  "Every close approach between the objects of ephemerides `eph-a` and
  `eph-b` from `t0` to `t1`, sampling each `step` seconds: `[{:t :miss
  :t-fit :miss-fit} ...]` -- the time and distance of each least
  separation found on the true ephemerides, and ANCAS's own estimate from
  the fitted cubics. With `:within` km, only those at least that close."
  ([eph-a eph-b t0 t1 step] (close-approaches eph-a eph-b t0 t1 step {}))
  ([eph-a eph-b t0 t1 step {:keys [within]}]
   (let [ts (conj (vec (range t0 t1 step)) t1)
         rel (fn [t] (let [[ra va] (eph-a t) [rb vb] (eph-b t)] [(v3/sub ra rb) (v3/sub va vb)]))
         samples (mapv (juxt identity rel) ts)]
     (vec (for [[[ta [pa va]] [tb [pb vb]]] (partition 2 1 samples)
                :let [h (- tb ta)
                      cs (hermite pa (v3/scale va h) pb (v3/scale vb h))]
                tau (interval-minima cs)
                :let [t-fit (+ ta (* tau h))
                      t (refine eph-a eph-b t-fit h)
                      [ra] (eph-a t) [rb] (eph-b t)
                      miss (v3/distance ra rb)]
                :when (or (nil? within) (<= miss within))]
            {:t t :miss miss :t-fit t-fit :miss-fit (v3/length (at cs tau))})))))
