(ns allgo.numerics.fft
  "The fast Fourier transform, radix 2.

  Here because a separable elliptic problem on a periodic axis stops being
  a linear system and becomes arithmetic once you change basis. A second
  difference around a ring is a *convolution*, so the discrete Fourier
  basis diagonalizes it: each wavenumber is scaled and no wavenumber talks
  to another. `allgo.physics.sphere-fluid` uses exactly that to turn a
  Poisson solve over a whole sphere into one small tridiagonal solve per
  wavenumber, which is direct, exact to round-off, and fast enough to run
  every frame.

  ## What makes it fast

  A transform of length `n` is a matrix of `n^2` products written down
  honestly. Split the sequence into its even and odd halves and each half
  is a transform of length `n/2` over the same roots of unity, so the work
  halves and the depth is `log n`. That is the whole of Cooley and Tukey.

  This is the iterative form rather than the recursive one: permute the
  input into bit-reversed order first, and the butterflies then run
  bottom-up over contiguous blocks with no call stack and no allocation.

  ## Conventions

  Complex data is two `double` arrays, real and imaginary, transformed in
  place. Not an array of pairs, and not interleaved: the arrays are the
  shape every caller here already has, and keeping them apart means a
  real-valued caller can hand over a zeroed imaginary array and never
  touch it again.

  `forward!` carries no scaling and `inverse!` carries all of it, so
  `(inverse! (forward! x))` is `x`. Lengths must be powers of two."
  (:require [allgo.array :as a]
            [clojure.math :as math]))

(defn power-of-two?
  "Whether `n` is a positive power of two -- one bit set, so `n` and `n-1`
  share none."
  [n]
  (let [n (long n)]
    (and (pos? n) (zero? (bit-and n (dec n))))))

(defn- twiddle-table
  "`exp(sign * 2*pi*i*j/n)` for `j` in `[0, n/2)`, as `[cos sin]`.

  One table serves every stage: the factor a stage of half-length `h`
  wants for its `k`th butterfly is this table at `k * n / (2h)`. Built
  from `cos` and `sin` directly rather than by repeatedly multiplying a
  root of unity, which is faster and drifts -- the error of the recurrence
  grows with the transform where this stays at one ulp."
  [n sign]
  (let [n (long n)
        h (quot n 2)
        c (a/f64 h)
        s (a/f64 h)
        step (* (double sign) 2.0 math/PI (/ 1.0 n))]
    (dotimes [j h]
      (aset c j (math/cos (* step j)))
      (aset s j (math/sin (* step j))))
    [c s]))

(def ^:private cached-table
  (memoize twiddle-table))

(defn- bit-reverse!
  "Permutes both arrays into bit-reversed index order.

  The butterflies below combine elements that are adjacent only after this
  shuffle; doing it up front is what lets the rest run over contiguous
  blocks."
  [^doubles re ^doubles im]
  (let [n (alength re)]
    (loop [i 0 j 0]
      (when (< i (dec n))
        (when (< i j)
          (let [tr (aget re i) ti (aget im i)]
            (aset re i (aget re j)) (aset im i (aget im j))
            (aset re j tr) (aset im j ti)))
        (let [j (long
                 (loop [j j k (quot n 2)]
                   (if (and (pos? k) (<= k j))
                     (recur (- j k) (quot k 2))
                     (+ j k))))]
          (recur (inc i) j))))))

(defn transform!
  "In-place unscaled transform. `sign` is -1 forward, +1 inverse.

  Returns the real array, so the pair can be threaded."
  ^doubles [^doubles re ^doubles im sign]
  (let [n (alength re)]
    (when-not (power-of-two? n)
      (throw (ex-info "FFT length must be a power of two" {:n n})))
    (when-not (= n (alength im))
      (throw (ex-info "real and imaginary parts must be the same length"
                      {:re n :im (alength im)})))
    (bit-reverse! re im)
    (let [[^doubles tc ^doubles ts] (cached-table n sign)]
      (loop [len 2]
        (when (<= len n)
          (let [h (quot len 2)
                stride (quot n len)]
            (loop [base 0]
              (when (< base n)
                (loop [k 0]
                  (when (< k h)
                    (let [t (* k stride)
                          wr (aget tc t) wi (aget ts t)
                          a (+ base k)
                          b (+ a h)
                          br (aget re b) bi (aget im b)
                          vr (- (* br wr) (* bi wi))
                          vi (+ (* br wi) (* bi wr))
                          ur (aget re a) ui (aget im a)]
                      (aset re a (+ ur vr))
                      (aset im a (+ ui vi))
                      (aset re b (- ur vr))
                      (aset im b (- ui vi)))
                    (recur (inc k))))
                (recur (+ base len)))))
          (recur (* 2 len)))))
    re))

(defn forward!
  "Forward transform, in place and unscaled."
  [^doubles re ^doubles im]
  (transform! re im -1.0))

(defn inverse!
  "Inverse transform, in place, carrying the whole `1/n`."
  [^doubles re ^doubles im]
  (transform! re im 1.0)
  (let [n (alength re)
        k (/ 1.0 n)]
    (dotimes [i n]
      (aset re i (* k (aget re i)))
      (aset im i (* k (aget im i)))))
  re)

(defn spectrum
  "The transform of a real sequence, as `[[re im] ...]`.

  A convenience for looking at a signal rather than solving with one;
  everything hot goes through `forward!` on arrays it already owns."
  [xs]
  (let [n (count xs)
        re (a/f64 (map double xs))
        im (a/f64 n)]
    (forward! re im)
    (mapv (fn [i] [(aget re i) (aget im i)]) (range n))))
