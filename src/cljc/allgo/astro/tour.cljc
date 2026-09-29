(ns allgo.astro.tour
  "Tours through the solar system by gravity assist: a sequence of
  planets and the dates at each, flown as patched conics -- Lambert's
  transfer between each pair, and at each planet between a hyperbola
  that bends the excess velocity from the leg arriving to the leg
  leaving.

  An unpowered flyby keeps the excess speed and only turns it, by no more
  than a close pass allows. Real tours rarely match both legs exactly, so
  each flyby here is powered: a burn at periapsis joins an incoming
  hyperbola of the arriving excess speed to an outgoing one of the
  leaving speed. Each hyperbola turns the velocity by asin(1/e) on its
  side of periapsis, e = 1 + rp v^2/mu, and the periapsis radius is the
  one at which the two halves add to the turn wanted; the burn is the
  difference of the two periapsis speeds. (The model is the usual one of
  multiple-gravity-assist design; see, e.g., Izzo, \"Global optimization
  and space pruning for spacecraft trajectory design\", in Conway, ed.,
  Spacecraft Trajectory Optimization, 2010.) If that periapsis would be
  inside the planet or its atmosphere, no burn at periapsis can do the
  turn and the flyby is infeasible.

  Dates are TT MJD; speeds km/s; distances km."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.interplanetary :as ip]
            [allgo.astro.iod :as iod]
            [allgo.astro.planets :as planets]
            [allgo.astro.rotation :as rotation]
            [allgo.astro.universal :as universal]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.optimize :as opt]
            [allgo.numerics.roots :as roots]
            [clojure.math :as math]))

(defn- bend
  "The turn of a hyperbola of excess speed `v` and periapsis `rp` on one
  side of periapsis: asin(1/e), 1/e = a/(a + rp) with a = mu/v^2."
  [mu v rp]
  (let [a (/ mu (* v v))] (math/asin (/ a (+ a rp)))))

(defn powered-flyby
  "The flyby about a body of parameter `mu` that takes the incoming excess
  velocity `v-in` to the outgoing `v-out`, both relative to the body:
  `{:turn :rp :dv}` -- the angle between them, the periapsis radius at
  which the two half-hyperbolas bend through it, and the burn at
  periapsis. `:rp` is nil when the turn is too great for any passage
  (it must be under 180 degrees)."
  [mu v-in v-out]
  (let [vi (v3/length v-in) vo (v3/length v-out)
        turn (math/acos (max -1.0 (min 1.0 (/ (v3/dot v-in v-out) (* vi vo)))))
        f #(- (+ (bend mu vi %) (bend mu vo %)) turn)
        ;; f falls from pi - turn at rp = 0 to -turn far out; bracket it
        hi (loop [hi (/ mu (* vi vo))] (if (pos? (f hi)) (recur (* 4.0 hi)) hi))
        rp (when (< turn (- math/PI 1e-12))
             (roots/bisect f 0.0 hi))]
    {:turn turn :rp rp
     :dv (when rp
           (abs (- (math/sqrt (+ (* vo vo) (/ (* 2.0 mu) rp)))
                   (math/sqrt (+ (* vi vi) (/ (* 2.0 mu) rp))))))}))

(def ^:private default-altitude
  "Least flyby altitudes, km: clear of the atmospheres, and for Jupiter of
  the worst of its radiation belts."
  {:mercury 200.0 :venus 300.0 :earth 300.0 :mars 200.0
   :jupiter 600000.0 :saturn 10000.0 :uranus 5000.0 :neptune 5000.0})

(defn gm
  "The gravitational parameter a spacecraft feels close to `body`: the
  Earth's own for the Earth, whose Moon is 1.2% of the pair, and each
  other planet's system's, where the moons are a ten-thousandth or less."
  [body]
  (if (= body :earth) c/GM-earth (c/GM-planet body)))

(defn min-radius
  "The least periapsis radius allowed at `body`: its equatorial radius and
  `altitude` (default a safe margin above its atmosphere)."
  ([body] (min-radius body (default-altitude body 300.0)))
  ([body altitude] (+ (rotation/radius body) altitude)))

(defn- flyby-cost
  "A flyby's burn, with the shortfall of an infeasible one weighted as
  `cost` weighs it."
  [{:keys [dv rp min-rp feasible?]}]
  (+ (or dv 0.0) (if feasible? 0.0 (* 10.0 (if rp (/ (- min-rp rp) min-rp) 1.0)))))

(defn- resonant-leg
  "A return to `body` after `revs` revolutions of the spacecraft, from
  `t1` to `t2`: Lambert's problem there is singular -- the planet is back
  nearly where it was, and the plane of the transfer is fixed by almost
  nothing -- so the leg is flown as mission designers fly resonant orbits.
  The spacecraft's period is (t2 - t1)/revs, which sets its heliocentric
  speed as it leaves; its direction, two angles, is chosen to make the
  flybys at either end (`flyby-in` and `flyby-out`, functions of the
  incoming and outgoing excess velocities giving `powered-flyby` maps)
  cheapest, given `v-in`, the excess velocity arriving at the first, and
  `v-out`, the one leaving the second. It comes back to where it left; how
  far the planet has moved from there, `:miss`, is the patched conic's
  error, zero when the flight time is a whole number of the planet's
  years."
  [body t1 t2 revs v-in v-out flyby-in flyby-out ephemeris]
  (let [[r1 vb1] (ephemeris body t1)
        [rb2 vb2] (ephemeris body t2)
        dt (* 86400.0 (- t2 t1))
        period (/ dt revs)
        a (math/cbrt (* c/GM-sun (math/pow (/ period (* 2.0 math/PI)) 2.0)))
        speed (math/sqrt (* c/GM-sun (- (/ 2.0 (v3/length r1)) (/ 1.0 a))))
        e1 (v3/normalize vb1)
        h (v3/normalize (v3/cross r1 vb1))
        e3 (v3/cross h e1)
        velocity (fn [[lon lat]]
                   (v3/scale (v3/add (v3/scale (v3/add (v3/scale e1 (math/cos lon)) (v3/scale e3 (math/sin lon)))
                                               (math/cos lat))
                                     (v3/scale h (math/sin lat)))
                             speed))
        leg (fn [x]
              (let [v1 (velocity x)
                    [_ v2] (universal/propagate c/GM-sun r1 v1 dt)]
                {:v1 v1 :v2 v2 :v-inf-depart (v3/sub v1 vb1) :v-inf-arrive (v3/sub v2 vb2)}))
        f (fn [x] (let [{:keys [v-inf-depart v-inf-arrive]} (leg x)]
                    (+ (flyby-cost (flyby-in v-in v-inf-depart))
                       (flyby-cost (flyby-out v-inf-arrive v-out)))))
        starts (for [lon (range -1.2 1.25 0.3) lat [-0.3 0.0 0.3]] [lon lat])
        {:keys [x]} (apply min-key :f (map #(opt/nelder-mead f % {:step 0.05 :tol-x 1e-9 :tol-f 1e-12})
                                           (take 3 (sort-by f starts))))]
    (assoc (leg x) :resonant? true :revs revs :miss (v3/distance r1 rb2))))

(defn- dsm-leg
  "A leg from `from` at `t1` to `to` at `t2` broken by one deep-space
  maneuver: Lambert's arc to a point in space at an intermediate time, a
  burn there, and Lambert's arc on to the planet. The point and the time
  are chosen to make the whole cheapest -- the burn, and `start-cost` and
  `end-cost`, functions of the leg's excess velocities leaving and
  arriving (the flybys' costs, or the launch's and arrival's). Searched by
  the simplex method from points along the direct transfer, the leg's
  `:long?` and `:revs` applying to the arc that covers more of the
  way. The map is `interplanetary/transfer`'s with `:dsm` `{:t :r :dv}`
  added."
  [from to t1 t2 start-cost end-cost {:keys [ephemeris] :as opts}]
  (let [[r1 vb1] (ephemeris from t1)
        [r2 vb2] (ephemeris to t2)
        direct (ip/transfer from t1 to t2 opts)
        ;; a point in space as radius (AU) and ecliptic-ish angles in the
        ;; frame of the departure planet's orbit
        e1 (v3/normalize r1)
        h (v3/normalize (v3/cross r1 vb1))
        e2 (v3/cross h e1)
        ->r (fn [[rad lon lat]]
              (v3/scale (v3/add (v3/scale (v3/add (v3/scale e1 (math/cos lon)) (v3/scale e2 (math/sin lon)))
                                          (math/cos lat))
                                (v3/scale h (math/sin lat)))
                        (* c/AU rad)))
        <-r (fn [r] (let [x (v3/dot r e1) y (v3/dot r e2) z (v3/dot r h)]
                      [(/ (v3/length r) c/AU) (math/atan2 y x) (math/asin (/ z (v3/length r)))]))
        leg (fn [[eta & p]]
              (when (< 0.02 eta 0.98)
                (let [tm (+ t1 (* eta (- t2 t1)))
                      rm (->r p)
                      first? (> eta 0.5)
                      arc (fn [ra rb dt main?]
                            (iod/lambert c/GM-sun ra rb (* 86400.0 dt)
                                         (if main? opts (dissoc opts :long? :revs))))
                      a1 (arc r1 rm (- tm t1) first?)
                      a2 (arc rm r2 (- t2 tm) (not first?))]
                  (when (and a1 a2)
                    (let [[v1 vm-] a1 [vm+ v2] a2]
                      {:v1 v1 :v2 v2 :v-inf-depart (v3/sub v1 vb1) :v-inf-arrive (v3/sub v2 vb2)
                       :dsm {:t tm :r rm :dv (v3/distance vm+ vm-)}})))))
        f (fn [x] (when-let [{:keys [v-inf-depart v-inf-arrive dsm]} (leg x)]
                    (+ (start-cost v-inf-depart) (:dv dsm) (end-cost v-inf-arrive))))
        starts (for [eta [0.25 0.4 0.55 0.7]
                     :let [[r] (when direct (universal/propagate c/GM-sun r1 (:v1 direct) (* 86400.0 eta (- t2 t1))))]
                     :when r]
                 (into [eta] (<-r r)))
        best (when (seq starts)
               (apply min-key :f (map #(opt/nelder-mead f % {:step [0.05 0.05 0.05 0.02] :tol-x 1e-7 :tol-f 1e-10})
                                      starts)))]
    (when (and best (< (:f best) ##Inf))
      (let [{:keys [v-inf-depart] :as l} (leg (:x best))]
        (assoc l :c3 (v3/dot v-inf-depart v-inf-depart))))))

(defn fly
  "The tour through `bodies` (keywords of `allgo.astro.planets`) at the
  dates `dates`, one per body, each leg the short way with no complete
  revolutions unless its entries in `:long?` and `:revs` (vectors, one
  per leg) say otherwise. A leg from a planet back to itself between two
  flybys a whole number of its years later, give or take 5%, as
  Galileo's Earth-Earth leg, is flown as a resonant orbit (see
  `resonant-leg`), `:revs` then the spacecraft's revolutions (default
  1). A leg whose entry in `:dsm` is true is broken by a deep-space
  maneuver (see `dsm-leg`); resonant and maneuvering legs must have plain
  legs, or the ends, either side. `:ephemeris` is
  `interplanetary/transfer`'s.

  The result's `:dsm-dv` sums the deep-space maneuvers. Returns `{:legs :flybys
  :c3 :v-inf-arrive :flyby-dv :feasible?}` -- each leg the map of
  `interplanetary/transfer`, each flyby `{:body :date :v-in :v-out :turn
  :rp :min-rp :altitude :dv :feasible?}`, the launch energy, the excess speed at
  the last body, the sum of the flyby burns and whether every flyby
  clears its planet (`:min-altitude`, a map of body to km, overrides the
  defaults). nil when a leg has no transfer."
  ([bodies dates] (fly bodies dates {}))
  ([bodies dates {:keys [long? revs dsm min-altitude ephemeris]
                  :or {ephemeris planets/heliocentric-state} :as opts}]
   (let [n (dec (count bodies))
         ;; back to the same planet between two flybys after nearly a
         ;; whole number of its years: Lambert's plane is undefined there
         resonant? (fn [i] (and (= (bodies i) (bodies (inc i))) (< 0 i (dec n))
                                (let [years (/ (- (dates (inc i)) (dates i)) (planets/period (bodies i) (dates i)))]
                                  (and (>= years 0.5) (< (abs (- years (math/round years))) 0.05)))))
         dsm? (fn [i] (get dsm i false))
         lambert-legs (mapv (fn [i]
                              (when-not (or (resonant? i) (dsm? i))
                                (ip/transfer (bodies i) (dates i) (bodies (inc i)) (dates (inc i))
                                             (assoc opts :long? (get long? i false) :revs (get revs i 0)))))
                            (range n))
         flyby (fn [i v-in v-out]
                 (let [body (bodies i)
                       {:keys [rp] :as f} (powered-flyby (gm body) v-in v-out)
                       rmin (min-radius body (get min-altitude body (default-altitude body 300.0)))]
                   (assoc f :body body :date (dates i) :v-in v-in :v-out v-out :min-rp rmin
                          :altitude (some-> rp (- (rotation/radius body)))
                          :feasible? (boolean (and rp (>= rp rmin))))))
         plain? #(not (or (resonant? %) (dsm? %)))
         legs (when (every? #(some? (lambert-legs %)) (filter plain? (range n)))
                (mapv (fn [i]
                        (cond
                          (dsm? i)
                          (dsm-leg (bodies i) (bodies (inc i)) (dates i) (dates (inc i))
                                   (if (zero? i)
                                     v3/length
                                     #(flyby-cost (flyby i (:v-inf-arrive (lambert-legs (dec i))) %)))
                                   (if (= i (dec n))
                                     v3/length
                                     #(flyby-cost (flyby (inc i) % (:v-inf-depart (lambert-legs (inc i))))))
                                   (assoc opts :long? (get long? i false) :revs (get revs i 0)
                                          :ephemeris ephemeris))
                          (resonant? i)
                          (resonant-leg (bodies i) (dates i) (dates (inc i)) (get revs i 1)
                                        (:v-inf-arrive (lambert-legs (dec i)))
                                        (:v-inf-depart (lambert-legs (inc i)))
                                        #(flyby i %1 %2) #(flyby (inc i) %1 %2) ephemeris)
                          :else (lambert-legs i)))
                      (range n)))]
     (when (and legs (every? some? legs))
       (let [flybys (mapv #(flyby % (:v-inf-arrive (legs (dec %))) (:v-inf-depart (legs %)))
                          (range 1 n))]
         {:bodies bodies :dates dates :legs legs :flybys flybys
          :c3 (:c3 (first legs))
          :v-inf-arrive (v3/length (:v-inf-arrive (peek legs)))
          :flyby-dv (reduce + 0.0 (keep :dv flybys))
          :dsm-dv (reduce + 0.0 (keep (comp :dv :dsm) legs))
          :feasible? (every? :feasible? flybys)})))))

(defn cost
  "A tour's total cost in km/s: the launch's excess speed (or, with
  `:parking` a radius about the first body, the burn out of that circular
  orbit), the flyby burns and deep-space maneuvers, and the arrival's
  excess speed -- or with `:capture` `[rp ra]` the burn into that orbit
  about the last body, or with `:arrival :none` nothing, for a tour that
  ends in a flyby. Penalties, `:penalty` km/s (default 10) for each
  unit, let a search climb back to what can be flown: an infeasible
  flyby's shortfall as a fraction of its least radius, and a resonant
  return's miss beyond half the planet's sphere of influence, as a
  fraction of the sphere."
  [{:keys [bodies legs flybys c3 v-inf-arrive flyby-dv dsm-dv]}
   {:keys [parking capture arrival penalty] :or {penalty 10.0}}]
  (let [v-dep (math/sqrt c3)
        depart (if parking (:dv (ip/departure (gm (first bodies)) parking v-dep)) v-dep)
        arrive (cond (= arrival :none) 0.0
                     capture (let [[rp ra] capture] (ip/capture (gm (peek bodies)) rp v-inf-arrive ra))
                     :else v-inf-arrive)
        short (reduce + 0.0 (for [{:keys [rp min-rp feasible?]} flybys :when (not feasible?)]
                              (if rp (/ (- min-rp rp) min-rp) 1.0)))
        miss (reduce + 0.0 (for [[i {:keys [miss]}] (map-indexed vector legs) :when miss]
                             (let [soi (ip/sphere-of-influence (bodies i))]
                               (max 0.0 (/ (- miss (* 0.5 soi)) soi)))))]
    (+ depart flyby-dv dsm-dv arrive (* penalty (+ short miss)))))

(defn optimize
  "The cheapest tour through `bodies` near the dates `dates`, found by the
  simplex method over the launch date and each leg's flight time, the
  cost `cost` with `opts`: the `fly` result with `:cost` added. Options
  are `fly`'s and `cost`'s, and `:step`, the initial simplex's edge in
  days (default 10)."
  ([bodies dates] (optimize bodies dates {}))
  ([bodies dates {:keys [step] :or {step 10.0} :as opts}]
   (let [->dates (fn [[t0 & tofs]] (vec (reductions + t0 tofs)))
         x0 (into [(first dates)] (map - (rest dates) dates))
         f (fn [x] (when (every? pos? (rest x))
                     (some-> (fly bodies (->dates x) opts) (cost opts))))
         {:keys [x]} (opt/nelder-mead f x0 {:step step :tol-x 1e-5 :tol-f 1e-9
                                            :max-iter (* 400 (count x0))})
         t (fly bodies (->dates x) opts)]
     (assoc t :cost (cost t opts)))))
