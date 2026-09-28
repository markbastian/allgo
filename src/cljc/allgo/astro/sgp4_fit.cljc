(ns allgo.astro.sgp4-fit
  "Making a two-line element set: the SGP4 mean elements that best
  reproduce a set of states, by differential correction (Vallado and
  Crawford, \"SGP4 Orbit Determination\", AIAA 2008-6770 -- the fit by
  which element sets are made, run here on SGP4 itself).

  Mean elements are not osculating ones: SGP4's are its own averaged
  quantities, defined only by the theory, so the only way to them from
  states is to search -- propagate trial elements with `sgp4`, compare,
  correct, again. Each correction is a Gauss-Newton step, the partials
  taken by central differences and the linear problem solved by
  orthogonal reduction (`allgo.astro.estimation`); a step that would
  raise the residual is halved until it does not.

  The elements are solved for in equinoctial form -- mean motion,
  af = e cos(w + W), ag = e sin(w + W), chi = tan(i/2) sin W,
  psi = tan(i/2) cos W and the mean longitude M + w + W -- which has no
  singularity at the circular or equatorial orbits so many satellites
  fly (only at i = 180 degrees); B* is solved for as well on request.
  The start is the osculating elements of the observation nearest the
  epoch, close enough for the search to take hold -- its velocity found
  by Lambert's problem where the observations are positions alone.

  Observations are `{:t :r :v}`: minutes from the epoch, TEME position
  in km and, optionally, velocity in km/s."
  (:require [allgo.astro.estimation :as est]
            [allgo.astro.iod :as iod]
            [allgo.astro.kepler :as kepler]
            [allgo.astro.sgp4 :as sgp4]
            [allgo.math :as am]
            [allgo.numerics.differentiation :as diff]
            [clojure.math :as math]))

(defn elements->equinoctial
  "[n af ag chi psi L] of SGP4's `:no-kozai :ecco :inclo :nodeo :argpo
  :mo`."
  [{:keys [no-kozai ecco inclo nodeo argpo mo]}]
  (let [lp (+ nodeo argpo)
        t (math/tan (* 0.5 inclo))]
    [no-kozai (* ecco (math/cos lp)) (* ecco (math/sin lp))
     (* t (math/sin nodeo)) (* t (math/cos nodeo)) (+ lp mo)]))

(defn equinoctial->elements
  "SGP4's `:no-kozai :ecco :inclo :nodeo :argpo :mo` of [n af ag chi psi
  L]."
  [[n af ag chi psi L]]
  (let [lp (math/atan2 ag af)
        node (math/atan2 chi psi)]
    {:no-kozai n
     :ecco (math/hypot af ag)
     :inclo (* 2.0 (math/atan (math/hypot chi psi)))
     :nodeo (am/wrap-2pi node)
     :argpo (am/wrap-2pi (- lp node))
     :mo (am/wrap-2pi (- L lp))}))

(defn- velocity-of
  "The velocity of observation `o`: its own, or, for positions alone, by
  Lambert's problem to another 5 to 30 minutes away -- an arc long
  enough that noise matters little and short enough to be under half a
  revolution for any satellite SGP4 carries."
  [mus obs {:keys [t r v] :as o}]
  (or v
      (let [near (remove #(identical? o %) obs)
            within (filter #(<= 5.0 (abs (- (:t %) t)) 30.0) near)
            o2 (if (seq within)
                 (apply max-key #(abs (- (:t %) t)) within)
                 (apply min-key #(abs (- (:t %) t)) near))
            dt (* 60.0 (- (:t o2) t))
            [v1 v2] (if (pos? dt)
                      (iod/lambert mus r (:r o2) dt {})
                      (iod/lambert mus (:r o2) r (- dt) {}))]
        (if (pos? dt) v1 v2))))

(defn- initial-elements
  "Osculating elements of the observation nearest the epoch, taken as
  mean ones at the epoch."
  [obs mus]
  (let [{:keys [t r] :as o} (apply min-key #(abs (:t %)) obs)
        v (velocity-of mus obs o)
        {:keys [a e i raan argp M]} (kepler/state->elements mus r v)
        n (* 60.0 (math/sqrt (/ mus (* a a a))))]
    {:no-kozai n :ecco e :inclo i :nodeo raan :argpo argp
     :mo (am/wrap-2pi (- M (* n t)))}))

(defn fit
  "The element record, from `sgp4init`, that best fits the observations
  `obs` in the least-squares sense: `{:satrec :elements :rms :iterations
  :covariance}`, :elements the fitted mean elements, :rms the weighted
  residual's root mean square per component and :covariance that of [n af
  ag chi psi L] (and B*), in the units of the weights. `template` gives
  what is not fitted -- :epoch, and optionally :ndot, :nddot, :bstar,
  :whichconst, :opsmode and any keys to carry -- and may give a starting
  :no-kozai :ecco :inclo :nodeo :argpo :mo. Options:

    :bstar?    solve for B* too (default false)
    :sigma-r   the position weight's standard deviation, km (default 1e-3)
    :sigma-v   the velocity's, km/s (default 1e-6)
    :max-iter  default 25
    :tol       stop when an iteration improves the rms by less than this
               fraction of it (default 1e-10)"
  ([template obs] (fit template obs {}))
  ([template obs {:keys [bstar? sigma-r sigma-v max-iter tol]
                  :or {bstar? false sigma-r 1e-3 sigma-v 1e-6 max-iter 25 tol 1e-10}}]
   (let [{:keys [mus]} (sgp4/getgravconst (:whichconst template :wgs72))
         start (if (:no-kozai template)
                 (select-keys template [:no-kozai :ecco :inclo :nodeo :argpo :mo])
                 (initial-elements obs mus))
         template (merge {:bstar 0.0 :ndot 0.0 :nddot 0.0} template)
         x0 (cond-> (elements->equinoctial start) bstar? (conj (:bstar template)))
         np (count x0)
         satrec-of (fn [x]
                     (sgp4/sgp4init (merge template (equinoctial->elements (subvec x 0 6))
                                           (when bstar? {:bstar (x 6)}))))
         ;; the observations, and the model's prediction of them, weighted
         weigh (fn [{:keys [r v]}]
                 (into (mapv #(/ % sigma-r) r) (map #(/ % sigma-v)) v))
         observed (vec (mapcat weigh obs))
         predict (fn [x]
                   (let [satrec (satrec-of x)]
                     (vec (mapcat (fn [{:keys [t v]}]
                                    (let [s (sgp4/sgp4 satrec t)]
                                      (if (and (:r s) (zero? (:error s)))
                                        (weigh {:r (:r s) :v (when v (:v s))})
                                        (repeat (if v 6 3) ##NaN))))
                                  obs))))
         rms-of (fn [pred]
                  (let [s (reduce + (map (fn [o p] (let [d (- o p)] (* d d))) observed pred))]
                    (if (NaN? s) ##Inf (math/sqrt (/ s (count observed))))))
         steps (cond-> [(* 1e-7 (first x0)) 1e-7 1e-7 1e-7 1e-7 1e-7] bstar? (conj 1e-7))]
     (loop [x x0 pred (predict x0) k 0]
       (let [rms (rms-of pred)
             J (diff/jacobian predict x {:steps steps :richardson? false})
             rows (mapv (fn [H o p] {:H H :residual (- o p)}) J observed pred)
             {:keys [correction covariance]} (est/solve rows np)
             ;; halve a step that would make matters worse
             [x' pred' rms'] (loop [h 1.0 tries 0]
                               (let [x' (mapv + x (map #(* h %) correction))
                                     pred' (predict x')
                                     rms' (rms-of pred')]
                                 (if (or (<= rms' rms) (>= tries 20))
                                   [x' pred' rms']
                                   (recur (* 0.5 h) (inc tries)))))
             done? (or (> rms' rms) (<= (- rms rms') (* tol rms)) (>= (inc k) max-iter))]
         (if done?
           (let [[x rms] (if (<= rms' rms) [x' rms'] [x rms])
                 satrec (satrec-of x)]
             {:satrec satrec
              :elements (select-keys satrec [:no-kozai :ecco :inclo :nodeo :argpo :mo :bstar])
              :rms rms
              :iterations (inc k)
              :covariance covariance})
           (recur x' pred' (inc k))))))))
