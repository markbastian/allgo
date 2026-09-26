(ns allgo.numerics.tridiagonal
  "Solving a tridiagonal system in linear time.

  The Thomas algorithm: Gaussian elimination with everything that is
  structurally zero left out. A general `n` by `n` solve is `O(n^3)`;
  a matrix with only three occupied diagonals needs one sweep down to
  eliminate the sub-diagonal and one sweep back to substitute, so it is
  `O(n)` and allocates two vectors.

  These turn up wherever a problem is one-dimensional and implicit: a
  second difference along a line, an implicit diffusion step, a natural
  cubic spline. `allgo.physics.sphere-fluid` gets one per Fourier
  wavenumber, which is what makes its Poisson solve direct rather than
  iterative.

  Arrays rather than the vectors of `allgo.numerics.linear`, because the
  systems here are long and solved in bulk -- a sphere solve runs one per
  wavenumber per step -- and because the diagonals are naturally separate
  arrays rather than a matrix anyone would write down.

  ## When it is safe

  Elimination without pivoting, which is not always allowed. It is
  guaranteed for diagonally dominant systems and for symmetric positive
  definite ones, which between them cover every system that comes from a
  second difference. `solve` returns nil rather than infinities if a
  pivot vanishes, so a caller that is unsure can find out."
  (:refer-clojure :exclude [solve])
  (:require [allgo.array :as a]))

(defn solve
  "Solves the tridiagonal system with the given diagonals, into `out`.

  `sub` is the sub-diagonal with `sub[0]` unused, `diag` the diagonal,
  `super` the super-diagonal with `super[n-1]` unused, and `rhs` the right
  hand side. `out` receives the solution and is returned; `c-scratch` is
  working space of the same length, supplied by the caller so that a loop
  over many systems allocates nothing.

  Returns nil if a pivot vanishes, leaving `out` in an unspecified state."
  [^doubles sub ^doubles diag ^doubles super ^doubles rhs ^doubles out ^doubles c-scratch]
  (let [n (alength diag)]
    (if (zero? n)
      out
      (let [d0 (aget diag 0)]
        (if (zero? d0)
          nil
          (do
            ;; Forward sweep: eliminate the sub-diagonal, carrying the
            ;; modified super-diagonal in scratch and the modified right
            ;; hand side in the output.
            (aset c-scratch 0 (/ (aget super 0) d0))
            (aset out 0 (/ (aget rhs 0) d0))
            (loop [i 1]
              (if (= i n)
                ;; Back substitution.
                (do
                  (loop [i (- n 2)]
                    (when (>= i 0)
                      (aset out i (- (aget out i) (* (aget c-scratch i) (aget out (inc i)))))
                      (recur (dec i))))
                  out)
                (let [m (- (aget diag i) (* (aget sub i) (aget c-scratch (dec i))))]
                  (if (zero? m)
                    nil
                    (do
                      (when (< i (dec n))
                        (aset c-scratch i (/ (aget super i) m)))
                      (aset out i (/ (- (aget rhs i) (* (aget sub i) (aget out (dec i)))) m))
                      (recur (inc i)))))))))))))

(defn solve-vectors
  "`solve` over ordinary sequences, returning a vector.

  For tests and for callers with one small system, where the allocation is
  beneath notice and the arrays are a nuisance."
  [sub diag super rhs]
  (let [n (count diag)
        arr (fn [xs] (a/f64 (map double xs)))
        out (a/f64 n)]
    (when (solve (arr sub) (arr diag) (arr super) (arr rhs) out (a/f64 n))
      (vec out))))

(defn solve-complex
  "Thomas for a system with real off-diagonals and a complex diagonal.

  The shape a semi-implicit wave term takes: discretizing a first
  derivative along a periodic axis and transforming it puts a pure
  imaginary number on the diagonal and leaves the off-diagonals, which
  came from a second derivative across the other axis, real. The
  elimination is unchanged -- only the arithmetic is complex, and only on
  the diagonal and the right hand side.

  Real and imaginary parts travel in separate arrays, as in
  `allgo.numerics.fft`. Returns `out-re`, or nil if a pivot vanishes."
  [^doubles sub ^doubles diag-re ^doubles diag-im ^doubles super
   ^doubles rhs-re ^doubles rhs-im ^doubles out-re ^doubles out-im
   ^doubles c-re ^doubles c-im]
  (let [n (alength diag-re)]
    (if (zero? n)
      out-re
      ;; One complex division, written out: (a+bi)/(c+di).
      ;; No primitive hints: seven arguments, and Clojure's primitive
      ;; function interfaces stop at four.
      (letfn [(divide! [ar ai br bi k ^doubles dst-re ^doubles dst-im]
                (let [ar (double ar) ai (double ai)
                      br (double br) bi (double bi) k (long k)
                      den (+ (* br br) (* bi bi))]
                  (if (zero? den)
                    false
                    (do (aset dst-re k (/ (+ (* ar br) (* ai bi)) den))
                        (aset dst-im k (/ (- (* ai br) (* ar bi)) den))
                        true))))]
        (if-not (and (divide! (aget super 0) 0.0 (aget diag-re 0) (aget diag-im 0) 0 c-re c-im)
                     (divide! (aget rhs-re 0) (aget rhs-im 0)
                              (aget diag-re 0) (aget diag-im 0) 0 out-re out-im))
          nil
          (loop [i 1]
            (if (= i n)
              (do (loop [i (- n 2)]
                    (when (>= i 0)
                      ;; out_i -= c_i * out_{i+1}, complex.
                      (let [cr (aget c-re i) ci (aget c-im i)
                            xr (aget out-re (inc i)) xi (aget out-im (inc i))]
                        (aset out-re i (- (aget out-re i) (- (* cr xr) (* ci xi))))
                        (aset out-im i (- (aget out-im i) (+ (* cr xi) (* ci xr)))))
                      (recur (dec i))))
                  out-re)
              (let [a (aget sub i)
                    mr (- (aget diag-re i) (* a (aget c-re (dec i))))
                    mi (- (aget diag-im i) (* a (aget c-im (dec i))))
                    rr (- (aget rhs-re i) (* a (aget out-re (dec i))))
                    ri (- (aget rhs-im i) (* a (aget out-im (dec i))))]
                (if-not (and (or (= i (dec n))
                                 (divide! (aget super i) 0.0 mr mi i c-re c-im))
                             (divide! rr ri mr mi i out-re out-im))
                  nil
                  (recur (inc i)))))))))))
