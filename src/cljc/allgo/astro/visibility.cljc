(ns allgo.astro.visibility
  "Who can see what: lines of sight past the Earth, the Earth's shadow and
  a satellite's passages through it, and when a satellite can be seen from
  the ground (Vallado, *Fundamentals of Astrodynamics and Applications*,
  chapter 5: SIGHT, SHADOW and LIGHT).

  A satellite is seen with the naked eye only when three things hold at
  once: it is above the observer's horizon, it is in sunlight, and the
  observer is not -- the sky must be dark enough, which in practice means
  the Sun at least a few degrees below the horizon. The shadow itself is
  `allgo.astro.srp/shadow`'s cone, umbra and penumbra both; the cylinder
  is here too, as the simpler model many analyses use.

  Positions are in any one frame centered on the Earth, km; times in
  whatever unit the caller's position functions take."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.srp :as srp]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.roots :as roots]
            [clojure.math :as math]))

;; ------------------------------------------------------------ sight

(def ^:private polar-flattening
  "The WGS-84 flattening, for the oblate line of sight."
  (/ 1.0 298.257223563))

(defn sight?
  "Whether `r1` and `r2` can see each other past the Earth: the segment
  between them clears a sphere of the Earth's radius -- or, with
  `:oblate? true`, the ellipsoid, whose poles are 21 km lower. The
  ellipsoid is turned into a sphere by stretching z by a/b, which keeps
  straight lines straight."
  ([r1 r2] (sight? r1 r2 {}))
  ([r1 r2 {:keys [oblate?]}]
   (let [k (if oblate? (/ 1.0 (- 1.0 polar-flattening)) 1.0)
         stretch (fn [[x y z]] [x y (* k z)])
         a (stretch r1) b (stretch r2)
         d (v3/sub b a)
         ;; the point of the segment nearest the center
         t (let [dd (v3/dot d d)] (if (zero? dd) 0.0 (- (/ (v3/dot a d) dd))))
         nearest (if (<= 0.0 t 1.0) (v3/add-scaled a d t) (if (< t 0.0) a b))
         R c/R-earth]
     (and (>= (v3/length a) R) (>= (v3/length b) R)
          (>= (v3/length nearest) R)))))

;; ----------------------------------------------------------- shadow

(defn shadow
  "Where `r` is in the Earth's shadow with the Sun at `r-sun`: `{:state
  :fraction}`, the state :sunlit, :penumbra or :umbra and the fraction of
  the Sun's disk visible -- the conical shadow of the Sun's finite disk.
  With `:model :cylinder` the shadow is instead a cylinder of the Earth's
  radius along the anti-Sun line, all or nothing."
  ([r r-sun] (shadow r r-sun {}))
  ([r r-sun {:keys [model] :or {model :cone}}]
   (if (= model :cylinder)
     (let [s (v3/normalize r-sun)
           along (v3/dot r s)
           across (v3/length (v3/add-scaled r s (- along)))
           dark? (and (neg? along) (< across c/R-earth))]
       {:state (if dark? :umbra :sunlit) :fraction (if dark? 0.0 1.0)})
     (let [f (srp/shadow r r-sun)]
       {:state (cond (>= f 1.0) :sunlit (<= f 0.0) :umbra :else :penumbra)
        :fraction f}))))

(defn sunlit?
  "Whether any of the Sun is visible from `r`."
  [r r-sun]
  (pos? (:fraction (shadow r r-sun))))

(defn eclipses
  "The satellite's passages through the Earth's shadow between `t0` and
  `t1`: `(position t)` gives its position and `(sun t)` the Sun's, in one
  frame. Returns `[{:penumbra-in :umbra-in :umbra-out :penumbra-out}
  ...]`, the times it enters and leaves each, nil where a passage starts
  or ends outside the interval (or never reaches the umbra). `step` must
  be shorter than the shortest penumbra crossing -- tens of seconds for a
  low orbit -- and `tol` is how finely each instant is found."
  [position sun t0 t1 step tol]
  (let [dark? (fn [t] (< (:fraction (shadow (position t) (sun t))) 1.0))
        umbra? (fn [t] (<= (:fraction (shadow (position t) (sun t))) 0.0))
        edges (sort-by first
                       (concat (for [[t was] (roots/transitions dark? t0 t1 step tol)]
                                 [t (if was :penumbra-out :penumbra-in)])
                               (for [[t was] (roots/transitions umbra? t0 t1 step tol)]
                                 [t (if was :umbra-out :umbra-in)])))
        start (when (dark? t0) {})]
    (loop [edges edges current start out []]
      (if-let [[[t kind] & more] (seq edges)]
        (case kind
          :penumbra-in (recur more {:penumbra-in t} (if current (conj out current) out))
          :penumbra-out (recur more nil (conj out (assoc (or current {}) :penumbra-out t)))
          (recur more (assoc (or current {}) kind t) out))
        (if current (conj out current) out)))))

;; ------------------------------------------------------------ light

(defn elevation
  "The elevation of `r` above the horizon of a site at `site`, radians, the
  horizon taken normal to the site's geocentric radius."
  [site r]
  (let [d (v3/sub r site)
        up (v3/normalize site)
        vertical (v3/dot d up)]
    ;; atan2 of the vertical and horizontal parts: asin loses its digits
    ;; near the zenith
    (math/atan2 vertical (v3/length (v3/add-scaled d up (- vertical))))))

(defn visible?
  "Whether a satellite at `r` can be seen with the naked eye from a site
  at `site`, the Sun at `r-sun`: above the horizon by `:mask` (default 0),
  in sunlight, and the Sun at least `:twilight` below the site's horizon
  (default 6 degrees, the end of civil twilight)."
  ([site r r-sun] (visible? site r r-sun {}))
  ([site r r-sun {:keys [mask twilight] :or {mask 0.0 twilight (* 6.0 c/degrees)}}]
   (and (> (elevation site r) mask)
        (sunlit? r r-sun)
        (< (elevation site r-sun) (- twilight)))))
