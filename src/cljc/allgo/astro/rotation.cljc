(ns allgo.astro.rotation
  "Which way the Sun, the planets and the Moon face: the IAU rotation
  models.

  Each body's spin is given as three angles in the ICRF, which for
  everything drawn here is the same frame as EME2000 to far better than a
  pixel: the right ascension and declination of its north pole, and `W`,
  how far its prime meridian has turned past the node where its equator
  crosses the ICRF's. The pole drifts slowly, as a polynomial in Julian
  centuries `T`; `W` runs at the spin rate, as a polynomial in days `d`;
  and some bodies add periodic terms, sines and cosines of angles that are
  themselves linear in `T` -- thirteen for the Moon, whose librations they
  are, a handful for the giant planets.

  The numbers are the IAU Working Group on Cartographic Coordinates and
  Rotational Elements' (Archinal et al. 2018, the 2015 report) as NASA's
  NAIF distributes them in the planetary constants kernel `pck00011.tpc`,
  transcribed by a script from that file's data sections rather than by
  hand. The Earth's entry is the IAU's simple model, good to a fraction
  of a degree and so for drawing, not for pointing an antenna; that is
  what `allgo.astro.frames` is for.

  North is the IAU's: the pole on the north side of the solar system's
  invariable plane, so Venus and Uranus, which turn backward, have `W`
  decreasing."
  (:require [clojure.math :as math]))

(def ^:private angles
  "The barycenters' nutation-precession angles, degrees, each as
  coefficients of powers of T."
  {1 [[174.7910857 149472.53587500003]
      [349.5821714 298945.07175000006]
      [164.3732571 448417.60762500006]
      [339.1643429 597890.1435000001]
      [153.9554286 747362.679375]]
   3 [[125.045 -1935.5364525]
      [250.089 -3871.072905]
      [260.008 475263.3328725]
      [176.625 487269.629985]
      [357.529 35999.0509575]
      [311.589 964468.49931]
      [134.963 477198.869325]
      [276.617 12006.300765]
      [34.226 63863.5132425]
      [15.134 -5806.6093575]
      [119.743 131.84064]
      [239.961 6003.1503825]
      [25.053 473327.79642]]
   4 [[190.72646643 15917.10818695 0.0]
      [21.4689247 31834.27934054 0.0]
      [332.86082793 19139.89694742 0.0]
      [394.93256437 38280.79631835 0.0]
      [189.6327156 41215158.1842005 12.711923222]
      [121.46893664 660.22803474 0.0]
      [231.05028581 660.9912354 0.0]
      [251.37314025 1320.50145245 0.0]
      [217.98635955 38279.9612555 0.0]
      [196.19729402 19139.83628608 0.0]
      [198.991226 19139.4819985 0.0]
      [226.292679 38280.8511281 0.0]
      [249.663391 57420.7251593 0.0]
      [266.18351 76560.636795 0.0]
      [79.398797 0.5042615 0.0]
      [122.433576 19139.9407476 0.0]
      [43.058401 38280.8753272 0.0]
      [57.663379 57420.7517205 0.0]
      [79.476401 76560.6495004 0.0]
      [166.325722 0.5042615 0.0]
      [129.071773 19140.0328244 0.0]
      [36.352167 38281.0473591 0.0]
      [56.668646 57420.929536 0.0]
      [67.364003 76560.2552215 0.0]
      [104.79268 95700.4387578 0.0]
      [95.391654 0.5042615 0.0]]
   5 [[73.32 91472.9]
      [24.62 45137.2]
      [283.9 4850.7]
      [355.8 1191.3]
      [119.9 262.1]
      [229.8 64.3]
      [352.25 2382.6]
      [113.35 6070.0]
      [146.64 182945.8]
      [49.24 90274.4]
      [99.360714 4850.4046]
      [175.895369 1191.9605]
      [300.323162 262.5475]
      [114.012305 6070.2476]
      [49.511251 64.3]]
   8 [[357.85 52.316]
      [323.92 62606.6]
      [220.51 55064.2]
      [354.27 46564.5]
      [75.31 26109.4]
      [35.36 14325.4]
      [142.61 2824.6]
      [177.85 52.316]
      [647.84 125213.2]
      [355.7 104.632]
      [533.55 156.948]
      [711.4 209.264]
      [889.25 261.58]
      [1067.1 313.896]
      [1244.95 366.212]
      [1422.8 418.528]
      [1600.65 470.844]]})

(def models
  "Each body's IAU rotation model, from the PCK: pole right ascension and
  declination as polynomials in T, prime meridian as one in d, the
  periodic terms against the angles of the barycenter named by
  `:system`, and the triaxial radii in km."
  {:sun
   {:radii [695700.0 695700.0 695700.0]
    :ra [286.13 0.0 0.0]
    :dec [63.87 0.0 0.0]
    :pm [84.176 14.1844 0.0]}
   :mercury
   {:radii [2440.53 2440.53 2438.26]
    :ra [281.0103 -0.0328 0.0]
    :dec [61.4155 -0.0049 0.0]
    :pm [329.5988 6.1385108 0.0]
    :system 1
    :ra-terms [0.0 0.0 0.0 0.0 0.0]
    :dec-terms [0.0 0.0 0.0 0.0 0.0]
    :pm-terms [0.01067257 -0.00112309 -0.0001104 -2.539e-05 -5.71e-06]}
   :venus
   {:radii [6051.8 6051.8 6051.8]
    :ra [272.76 0.0 0.0]
    :dec [67.16 0.0 0.0]
    :pm [160.2 -1.4813688 0.0]}
   :earth
   {:radii [6378.1366 6378.1366 6356.7519]
    :ra [0.0 -0.641 0.0]
    :dec [90.0 -0.557 0.0]
    :pm [190.147 360.9856235 0.0]}
   :moon
   {:radii [1737.4 1737.4 1737.4]
    :ra [269.9949 0.0031 0.0]
    :dec [66.5392 0.013 0.0]
    :pm [38.3213 13.17635815 -1.4e-12]
    :system 3
    :ra-terms [-3.8787 -0.1204 0.07 -0.0172 0.0 0.0072 0.0 0.0 0.0 -0.0052 0.0 0.0 0.0043]
    :dec-terms [1.5419 0.0239 -0.0278 0.0068 0.0 -0.0029 0.0009 0.0 0.0 0.0008 0.0 0.0 -0.0009]
    :pm-terms [3.561 0.1208 -0.0642 0.0158 0.0252 -0.0066 -0.0047 -0.0046 0.0028 0.0052 0.004 0.0019 -0.0044]}
   :mars
   {:radii [3396.19 3396.19 3376.2]
    :ra [317.269202 -0.10927547 0.0]
    :dec [54.432516 -0.05827105 0.0]
    :pm [176.049863 350.891982443297 0.0]
    :system 4
    :ra-terms [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 6.8e-05 0.000238 5.2e-05 9e-06 0.419057]
    :dec-terms [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.1e-05 0.000141 3.1e-05 5e-06 1.591274]
    :pm-terms [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.000145 0.000157 4e-05 1e-06 1e-06 0.584542]}
   :jupiter
   {:radii [71492.0 71492.0 66854.0]
    :ra [268.056595 -0.006499 0.0]
    :dec [64.495303 0.002413 0.0]
    :pm [284.95 870.536 0.0]
    :system 5
    :ra-terms [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.000117 0.000938 0.001432 3e-05 0.00215]
    :dec-terms [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 5e-05 0.000404 0.000617 -1.3e-05 0.000926]
    :pm-terms [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]}
   :saturn
   {:radii [60268.0 60268.0 54364.0]
    :ra [40.589 -0.036 0.0]
    :dec [83.537 -0.004 0.0]
    :pm [38.9 810.7939024 0.0]}
   :uranus
   {:radii [25559.0 25559.0 24973.0]
    :ra [257.311 0.0 0.0]
    :dec [-15.175 0.0 0.0]
    :pm [203.81 -501.1600928 0.0]}
   :neptune
   {:radii [24764.0 24764.0 24341.0]
    :ra [299.36 0.0 0.0]
    :dec [43.46 0.0 0.0]
    :pm [249.978 541.1397757 0.0]
    :system 8
    :ra-terms [0.7 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    :dec-terms [-0.51 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    :pm-terms [-0.48 0.0 0.0 0.0 0.0 0.0 0.0 0.0]}})

(defn- poly [cs x] (reduce (fn [acc c] (+ (* acc x) c)) 0.0 (rseq (vec cs))))

(defn angles-at
  "The pole's right ascension and declination and the prime meridian's
  angle `W` for `body` at `mjd` (TDB, for which TT will do), in radians."
  [body mjd]
  (let [{:keys [ra dec pm system ra-terms dec-terms pm-terms]} (models body)
        d  (- (double mjd) 51544.5)
        T  (/ d 36525.0)
        th (mapv #(math/to-radians (poly % T)) (get angles system))
        periodic (fn [terms f]
                   (reduce + 0.0 (map (fn [c t] (* c (f t))) terms th)))]
    {:ra  (math/to-radians (+ (poly ra T) (periodic ra-terms math/sin)))
     :dec (math/to-radians (+ (poly dec T) (periodic dec-terms math/cos)))
     :w   (math/to-radians (mod (+ (poly pm d) (periodic pm-terms math/sin)) 360.0))}))

(defn- rz [a] (let [c (math/cos a) s (math/sin a)] [[c (- s) 0.0] [s c 0.0] [0.0 0.0 1.0]]))
(defn- rx [a] (let [c (math/cos a) s (math/sin a)] [[1.0 0.0 0.0] [0.0 c (- s)] [0.0 s c]]))

(defn- mat-mul [a b]
  (mapv (fn [row] (mapv (fn [j] (reduce + (map * row (map #(nth % j) b)))) (range 3))) a))

(defn body->icrf
  "The rotation taking body-fixed coordinates -- x through the prime
  meridian, z through the north pole -- into the ICRF, as rows of a 3x3
  matrix: Rz(ra + 90 deg) Rx(90 deg - dec) Rz(W)."
  [body mjd]
  (let [{:keys [ra dec w]} (angles-at body mjd)]
    (mat-mul (mat-mul (rz (+ ra (/ math/PI 2.0))) (rx (- (/ math/PI 2.0) dec))) (rz w))))

(defn pole
  "The north pole's direction in the ICRF, a unit vector."
  [body mjd]
  (let [{:keys [ra dec]} (angles-at body mjd)]
    [(* (math/cos dec) (math/cos ra)) (* (math/cos dec) (math/sin ra)) (math/sin dec)]))

(defn sidereal-day
  "How long `body` takes to turn once against the stars, in hours;
  negative when it turns backward."
  [body]
  (/ (* 24.0 360.0) (second (:pm (models body)))))

(defn radius
  "The equatorial radius, km."
  [body]
  (first (:radii (models body))))
