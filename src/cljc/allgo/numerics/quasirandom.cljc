(ns allgo.numerics.quasirandom
  "Sampling that covers a space more evenly than chance: LaValle,
  *Planning Algorithms*, Cambridge University Press, 2006, sections 5.2.3
  and 5.2.4 -- the samples a roadmap planner or an integral is better
  served by.

    van der Corput  the integers' digits in a base, reflected about the
                    point: 1/2, 1/4, 3/4, 1/8, ... -- each new point in
                    the largest gap
    Halton          a van der Corput sequence per dimension, in coprime
                    bases; open-ended, low discrepancy in any dimension
    Hammersley      for a known count n: i / n in the first dimension,
                    Halton in the rest -- lower discrepancy still
    Sukharev grid   the grid of cell centers, the least dispersion of any
                    n points in the unit cube

  and the measures: dispersion, the radius of the largest empty ball --
  how far any point can be from the nearest sample -- and the star
  discrepancy, the worst difference between the share of points in a box
  from the origin and the box's volume. All in the unit cube [0, 1]^d.")

(def primes [2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53])

(defn van-der-corput
  "The `i`th term of the van der Corput sequence in `base`, from i = 1."
  ^double [^long i ^long base]
  (loop [i i f (/ 1.0 base) r 0.0]
    (if (zero? i) r (recur (quot i base) (/ f base) (+ r (* f (rem i base)))))))

(defn halton
  "The first `n` points of the Halton sequence in `d` dimensions, bases
  the first d primes."
  [n d]
  (vec (for [i (range 1 (inc n))] (mapv #(van-der-corput i (primes %)) (range d)))))

(defn hammersley
  "The `n` points of the Hammersley set in `d` dimensions: (i + 1/2) / n,
  then Halton's first d - 1 coordinates."
  [n d]
  (vec (for [i (range n)] (into [(/ (+ i 0.5) n)] (map #(van-der-corput (inc i) (primes %)) (range (dec d)))))))

(defn sukharev
  "The Sukharev grid of `k` points a side in `d` dimensions: the centers
  of the k^d equal cells."
  [k d]
  (letfn [(grid [d] (if (zero? d) [[]] (for [p (grid (dec d)) i (range k)] (conj p (/ (+ i 0.5) k)))))]
    (vec (grid d))))

(defn uniform
  "`n` points drawn uniformly at random, for comparison."
  [n d rng]
  (vec (repeatedly n (fn [] (vec (repeatedly d #(double (rng))))))))

(defn dispersion
  "The dispersion of `points`, estimated over `probes` (points of the cube
  to test, a fine grid say): the largest distance from a probe to its
  nearest point, under the max norm, as LaValle's definition allows."
  [points probes]
  (reduce max (for [q probes]
                (reduce min (map (fn [p] (reduce max (map #(abs (- %1 %2)) p q))) points)))))

(defn star-discrepancy
  "The star discrepancy of `points`, estimated over boxes [0, corner) for
  every corner in `corners`: the largest difference between the share of
  points in the box and its volume."
  [points corners]
  (let [n (count points)]
    (reduce max (for [c corners]
                  (abs (- (/ (count (filter (fn [p] (every? true? (map < p c))) points)) (double n))
                          (reduce * c)))))))
