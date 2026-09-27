(ns allgo.astro.stars
  "The naked-eye sky: the Bright Star Catalogue, in the frame everything
  else here is in.

  The catalogue is Hoffleit and Warren's fifth revised edition, VizieR
  V/50 -- nine thousand stars, every one brighter than about magnitude
  6.5, which is as faint as an eye sees on a dark night. It gives each a
  place at J2000 on the FK5 system, and FK5 J2000 is the mean equator and
  equinox of J2000 -- EME2000, the frame `allgo.astro` already refers the
  planets, the Moon and every satellite to. So a star needs no turning to
  sit behind them; its direction is simply the unit vector at its right
  ascension and declination. (FK5 and the ICRS differ by a few
  hundredths of an arcsecond, which no picture here can show.)

  The places are for J2000 and the stars move: `direction` carries each
  along its proper motion to any date, linearly, which is good for
  centuries. The fastest here, Arcturus, crosses a Moon's width in about
  eight hundred years.

  `color` is what a star of a B-V index would look like: B-V gives a
  temperature by Ballesteros's fit to blackbodies (2012), and the
  temperature a color by the usual fit to the blackbody's chromaticity.
  `flux` is brightness on the linear scale, from the magnitude's
  logarithmic one: five magnitudes is a factor of a hundred."
  (:require [allgo.math :as am]
            [clojure.math :as math]
            [clojure.string :as str]))

;; ---------------------------------------------------------------------------
;; The catalogue

(defn- number [s]
  (when-not (str/blank? s) (parse-double (str/trim s))))

(defn parse
  "The stars in the text of the catalogue, as `resources/public/data/bsc5.tsv`
  lays it out: a header of `#` comments and column names, then one star a
  line, tab-separated. Each is

    {:hr 2491 :name \"9Alp CMa\" :ra radians :dec radians :vmag -1.46
     :bv 0.0 :sp \"A1Vm\" :pm-ra radians/yr :pm-dec radians/yr}

  `:name` is the catalogue's own compact Flamsteed-and-Bayer form, nil
  where a star has none; `:bv` and the proper motions are nil where the
  catalogue gives none. `:pm-ra` is the projected motion, `cos(dec)
  dRA/dt`, as the catalogue's is."
  [text]
  (let [arcsec (/ math/PI 180.0 3600.0)
        deg (/ math/PI 180.0)]
    (into []
          (keep (fn [line]
                  (let [[hr nm ra dec v bv sp pmra pmde] (str/split line #"\t" -1)]
                    (when (and (seq hr) (re-matches #"\d+" hr) (number ra))
                      {:hr (parse-long hr)
                       :name (when-not (str/blank? nm) nm)
                       :ra (* deg (number ra))
                       :dec (* deg (number dec))
                       :vmag (number v)
                       :bv (number bv)
                       :sp (when-not (str/blank? sp) sp)
                       :pm-ra (some-> (number pmra) (* arcsec))
                       :pm-dec (some-> (number pmde) (* arcsec))}))))
          (str/split-lines text))))

(defn parse-names
  "The proper names in `resources/public/data/star-names.tsv` -- the IAU's
  approved names for Bright Star Catalogue stars -- as a map from HR
  number to name."
  [text]
  (into {}
        (keep (fn [line]
                (let [[hr nm] (str/split line #"\t")]
                  (when (and hr (re-matches #"\d+" hr) (not (str/blank? nm)))
                    [(parse-long hr) nm]))))
        (str/split-lines text)))

(defn parse-constellations
  "The constellations in `resources/public/data/constellations.tsv`, each

    {:abbr \"Ori\" :name \"Orion\" :label [ra dec] :lines [[[ra dec] ...] ...]}

  in radians: where to write its name, and its figure as polylines of
  places on the sky, to be drawn star to star."
  [text]
  (let [deg (/ math/PI 180.0)
        point (fn [s] (let [[a d] (str/split s #",")] [(* deg (number a)) (* deg (number d))]))]
    (into []
          (keep (fn [line]
                  (let [[abbr nm lra ldec lines] (str/split line #"\t" -1)]
                    (when (and (seq abbr) (not (str/starts-with? abbr "#"))
                               (not= "Abbr" abbr) (number lra))
                      {:abbr abbr
                       :name nm
                       :label [(* deg (number lra)) (* deg (number ldec))]
                       :lines (if (str/blank? lines)
                                []
                                (mapv (fn [pl] (mapv point (str/split pl #";")))
                                      (str/split lines #"\|")))}))))
          (str/split-lines text))))

;; ---------------------------------------------------------------------------
;; Where

(defn unit
  "The unit vector at right ascension `ra` and declination `dec`, radians."
  [ra dec]
  (let [cd (math/cos (double dec))]
    [(* cd (math/cos (double ra))) (* cd (math/sin (double ra))) (math/sin (double dec))]))

(defn direction
  "The unit vector toward `star` in EME2000, at `mjd` -- or at J2000, the
  catalogue's own epoch, without one.

  Proper motion is applied as a straight line on the sky, which is what
  the catalogue's two rates are: over the centuries this is used for the
  curvature of a star's true path is far under anything drawn."
  ([star] (direction star nil))
  ([{:keys [ra dec pm-ra pm-dec]} mjd]
   (let [years (if mjd (/ (- (double mjd) 51544.5) 365.25) 0.0)
         dec (+ (double dec) (* (double (or pm-dec 0.0)) years))
         ;; The catalogue's RA rate is along the great circle, so the
         ;; change in RA itself is that over cos(dec).
         cd (math/cos dec)
         ra (+ (double ra) (if (> (abs cd) 1e-9)
                             (/ (* (double (or pm-ra 0.0)) years) cd)
                             0.0))]
     (unit ra dec))))

;; ---------------------------------------------------------------------------
;; How it looks

(defn temperature
  "A star's effective temperature in kelvin from its B-V color index, by
  Ballesteros's fit to blackbodies (2012), good from the hottest stars to
  about B-V 2. Nil B-V -- the catalogue lacks it for a few -- is taken as
  a white star."
  ^double [bv]
  (let [bv (double (if (nil? bv) 0.4 bv))]
    (* 4600.0 (+ (/ 1.0 (+ (* 0.92 bv) 1.7))
                 (/ 1.0 (+ (* 0.92 bv) 0.62))))))

(defn color
  "The color a star of B-V index `bv` shows, as `[r g b]` in 0..1: its
  temperature's blackbody color, by the usual fit to the Planckian locus
  in sRGB (Tanner Helland's), then scaled so its brightest channel is 1.
  Red for Betelgeuse, blue-white for Rigel, white for Vega."
  [bv]
  (let [t (/ (temperature bv) 100.0)
        r (if (<= t 66.0) 255.0 (am/clamp (* 329.698727446 (math/pow (- t 60.0) -0.1332047592)) 0.0 255.0))
        g (am/clamp (if (<= t 66.0)
                      (- (* 99.4708025861 (math/log t)) 161.1195681661)
                      (* 288.1221695283 (math/pow (- t 60.0) -0.0755148492))) 0.0 255.0)
        b (cond (>= t 66.0) 255.0
                (<= t 19.0) 0.0
                :else (am/clamp (- (* 138.5177312231 (math/log (- t 10.0))) 305.0447927307) 0.0 255.0))
        m (max r g b)]
    [(/ r m) (/ g m) (/ b m)]))

(defn flux
  "Brightness relative to a star of magnitude zero: five magnitudes is a
  factor of a hundred, and brighter is smaller."
  ^double [vmag]
  (math/pow 10.0 (* -0.4 (double vmag))))
