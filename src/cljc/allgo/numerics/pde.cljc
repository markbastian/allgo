(ns allgo.numerics.pde
  "Partial differential equations by finite differences (Chapra and
  Canale, *Numerical Methods for Engineers*, chapters 29-30).

  Elliptic: Laplace's equation on a rectangular plate by Liebmann's method
  -- Gauss-Seidel on the five-point difference equations, over-relaxed --
  with fixed edges, edges of given derivative (insulated ones among
  them, handled by the imaginary node past the edge), and irregular
  boundaries that cut a grid line short of the next node (the Laplacian
  with unequal arms, eq. 29.24); and the heat flux, the secondary
  variable, from the solution.

  Parabolic: the heat-conduction equation dT/dt = k d2T/dx2 in one
  dimension by the explicit method (stable only for k dt/dx^2 <= 1/2),
  the simple implicit method and Crank-Nicolson (both tridiagonal and
  unconditionally stable, the latter second order in time too); in two by
  the explicit method and the alternating-direction implicit (ADI) scheme,
  each half step implicit along one axis only, so tridiagonal."
  (:require [allgo.numerics.tridiagonal :as tri]
            [clojure.math :as math]))

;; ---------------------------------------------------------------- elliptic

(defn- side-value [v k] (if (fn? v) (v k) (double v)))

(defn liebmann
  "The plate's temperatures by Liebmann's method. The plate has `:nx` by
  `:ny` intervals of `:dx` by `:dy` (default 1 each); the grid is indexed
  `[j i]`, j = 0 the bottom edge and i = 0 the left. Options:

    :boundary    {:top :bottom :left :right} the fixed edge temperatures,
                 each a number or a function of the node index along it
    :derivative  {side g} edges held at dT/dn = g instead (outward normal
                 for :right and :top, inward for :left and :bottom, as
                 dT/dx and dT/dy) -- 0 for an insulated edge; their nodes
                 become unknowns, balanced against an imaginary node
                 beyond the edge
    :arms        (fn [i j dir]) -> nil, or [fraction value]: the boundary
                 met at that fraction of the spacing in direction dir
                 (:left :right :down :up) with that temperature, for an
                 irregular boundary (Chapra and Canale eq. 29.24)
    :inside?     (fn [i j]) -> whether a node is part of the plate
                 (default every interior node), for irregular shapes
    :lambda      over-relaxation weight (default 1.5)
    :tol         stop when no node changes by more than this relative
                 (default 1e-10)

  Returns `{:T :iterations}`."
  [{:keys [nx ny dx dy boundary derivative arms inside? lambda tol max-iter]
    :or {dx 1.0 dy 1.0 lambda 1.5 tol 1e-10 max-iter 100000}}]
  (let [{:keys [top bottom left right]} boundary
        dv (or derivative {})
        free-left? (contains? dv :left) free-right? (contains? dv :right)
        free-bottom? (contains? dv :bottom) free-top? (contains? dv :top)
        T0 (vec (for [j (range (inc ny))]
                  (vec (for [i (range (inc nx))]
                         (cond
                           (and (= j 0) bottom (not free-bottom?)) (side-value bottom i)
                           (and (= j ny) top (not free-top?)) (side-value top i)
                           (and (= i 0) left (not free-left?)) (side-value left j)
                           (and (= i nx) right (not free-right?)) (side-value right j)
                           :else 0.0)))))
        unknown? (fn [i j]
                   (and (or (< 0 i nx) (and (= i 0) free-left?) (and (= i nx) free-right?))
                        (or (< 0 j ny) (and (= j 0) free-bottom?) (and (= j ny) free-top?))
                        (or (nil? inside?) (inside? i j))))
        nodes (vec (for [j (range (inc ny)) i (range (inc nx)) :when (unknown? i j)] [i j]))
        rx (/ 1.0 (* dx dx)) ry (/ 1.0 (* dy dy))
        ;; one node's new value from its neighbors
        update-node
        (fn [T i j]
          (let [arm (fn [dir] (when arms (arms i j dir)))
                ;; each direction: [arm length in spacings, neighbor value]
                nb (fn [dir di dj g-side]
                     (if-let [[a v] (arm dir)]
                       [a v]
                       (let [ii (+ i di) jj (+ j dj)]
                         (if (and (<= 0 ii nx) (<= 0 jj ny))
                           [1.0 (get-in T [jj ii])]
                           ;; past a derivative edge: the imaginary node's mirror
                           (let [g (dv g-side)
                                 h (if (#{:left :right} g-side) dx dy)
                                 mirror (get-in T [(- j dj) (- i di)])]
                             [1.0 (+ mirror (* 2.0 h g (if (#{:left :bottom} g-side) -1.0 1.0)))])))))
                [al vl] (nb :left -1 0 :left) [ar vr] (nb :right 1 0 :right)
                [ad vd] (nb :down 0 -1 :bottom) [au vu] (nb :up 0 1 :top)
                ;; eq. 29.24: the Laplacian with arms a1 a2 in x and b1 b2 in y
                cl (/ (* 2.0 rx) (* al (+ al ar))) cr (/ (* 2.0 rx) (* ar (+ al ar)))
                cd (/ (* 2.0 ry) (* ad (+ ad au))) cu (/ (* 2.0 ry) (* au (+ ad au)))]
            (/ (+ (* cl vl) (* cr vr) (* cd vd) (* cu vu)) (+ cl cr cd cu))))]
    (loop [T T0 k 1]
      (let [[T' change]
            (reduce (fn [[T change] [i j]]
                      (let [old (get-in T [j i])
                            new (+ (* lambda (update-node T i j)) (* (- 1.0 lambda) old))]
                        [(assoc-in T [j i] new)
                         (max change (if (zero? new) (abs (- new old)) (abs (/ (- new old) new))))]))
                    [T 0.0] nodes)]
        (if (or (< change tol) (>= k max-iter))
          {:T T' :iterations k}
          (recur T' (inc k)))))))

(defn flux
  "The heat flux at the interior node (i, j) of the plate `T`, `k'` the
  conductivity: q = -k' grad T by centered differences, `{:qx :qy :qn
  :theta}` -- the components, the magnitude and the direction (radians
  from the x axis) (Chapra and Canale 29.2.3)."
  [T i j k' dx dy]
  (let [qx (- (* k' (/ (- (get-in T [j (inc i)]) (get-in T [j (dec i)])) (* 2.0 dx))))
        qy (- (* k' (/ (- (get-in T [(inc j) i]) (get-in T [(dec j) i])) (* 2.0 dy))))]
    {:qx qx :qy qy :qn (math/hypot qx qy) :theta (math/atan2 qy qx)}))

;; ------------------------------------------------------ one dimension, in time

(defn- tridiagonal [sub diag super rhs] (tri/solve-vectors sub diag super rhs))

(defn heat-1d
  "The temperature of a rod, dT/dt = k d2T/dx2, from `:initial` (a number
  or a function of x) on `:n` intervals of `:dx`, stepping `:dt` to
  `:t-end`, by `:scheme` -- `:explicit` (forward in time, centered in
  space; unstable if k dt/dx^2 > 1/2), `:implicit` (backward in time) or
  `:crank-nicolson` (centered in time). The ends are held at `:left` and
  `:right` (numbers, or functions of t), or, with `:derivative {side g}`,
  at dT/dx = g there. Returns `[[t [T0 ... Tn]] ...]` every step."
  [{:keys [k n dx dt t-end initial left right derivative scheme]
    :or {scheme :explicit derivative {}}}]
  (let [lam (/ (* k dt) (* dx dx))
        at (fn [v t] (if (fn? v) (v t) (double v)))
        steps (long (math/round (/ t-end dt)))
        T0 (mapv #(if (fn? initial) (initial (* % dx)) (double initial)) (range (inc n)))
        T0 (cond-> T0
             (not (contains? derivative :left)) (assoc 0 (at left 0.0))
             (not (contains? derivative :right)) (assoc n (at right 0.0)))
        dl (contains? derivative :left) dr (contains? derivative :right)
        ;; the unknowns: interior nodes, and the derivative ends
        lo (if dl 0 1) hi (if dr n (dec n))
        ;; the value an end node's missing neighbor takes (imaginary node)
        ghost (fn [T side] (if (= side :left)
                             (- (T 1) (* 2.0 dx (derivative :left)))
                             (+ (T (dec n)) (* 2.0 dx (derivative :right)))))
        nbrs (fn [T i] [(if (= i 0) (ghost T :left) (T (dec i))) (if (= i n) (ghost T :right) (T (inc i)))])
        step
        (fn [T t]
          (let [t1 (+ t dt)
                T (cond-> T
                    (not dl) (assoc 0 (at left t))
                    (not dr) (assoc n (at right t)))
                fixed (fn [T'] (cond-> T' (not dl) (assoc 0 (at left t1)) (not dr) (assoc n (at right t1))))]
            (case scheme
              :explicit
              (fixed (reduce (fn [T' i] (let [[a b] (nbrs T i)] (assoc T' i (+ (T i) (* lam (+ a (* -2.0 (T i)) b))))))
                             T (range lo (inc hi))))
              (:implicit :crank-nicolson)
              (let [cn? (= scheme :crank-nicolson)
                    ;; implicit: -l T_i-1 + (1 + 2l) T_i - l T_i+1 = T_i
                    ;; Crank-Nicolson: -l T_i-1 + 2(1 + l) T_i - l T_i+1 = l T_i-1 + 2(1 - l) T_i + l T_i+1
                    d (if cn? (* 2.0 (+ 1.0 lam)) (+ 1.0 (* 2.0 lam)))
                    idx (vec (range lo (inc hi)))
                    m (count idx)
                    rhs0 (mapv (fn [i] (if cn?
                                         (let [[a b] (nbrs T i)] (+ (* lam a) (* 2.0 (- 1.0 lam) (T i)) (* lam b)))
                                         (T i)))
                               idx)
                    Tn (fixed T)
                    ;; the known ends at t1 move to the right side; a derivative
                    ;; end's imaginary node doubles its inner neighbor's weight
                    sub (mapv (fn [k] (let [i (idx k)] (cond (and (= i n) dr) (* -2.0 lam) :else (- lam)))) (range m))
                    super (mapv (fn [k] (let [i (idx k)] (cond (and (= i 0) dl) (* -2.0 lam) :else (- lam)))) (range m))
                    rhs (cond-> rhs0
                          (not dl) (update 0 + (* lam (Tn 0)))
                          (not dr) (update (dec m) + (* lam (Tn n)))
                          dl (update 0 - (* 2.0 lam dx (derivative :left)))
                          dr (update (dec m) + (* 2.0 lam dx (derivative :right))))
                    sol (tridiagonal sub (vec (repeat m d)) super rhs)]
                (reduce (fn [T' k] (assoc T' (idx k) (sol k))) Tn (range m))))))]
    (vec (reductions (fn [[t T] _] [(+ t dt) (step T t)]) [0.0 T0] (range steps)))))

;; ------------------------------------------------------ two dimensions, in time

(defn heat-2d
  "The temperature of a plate, dT/dt = k (d2T/dx2 + d2T/dy2), on `:nx` by
  `:ny` intervals of `:dx` (equal in both directions), from `:initial`
  (default 0) with the edges held at `:boundary` {:top :bottom :left
  :right} (numbers), stepping `:dt` to `:t-end` by `:scheme` -- `:explicit`
  (stable for k dt/dx^2 <= 1/4) or `:adi` (Peaceman and Rachford: each
  step two halves, the first implicit along y and explicit along x, the
  second the other way, each a set of tridiagonal systems; Chapra and
  Canale 30.5.2). Returns `[[t T] ...]`, T indexed `[j i]`.

  The explicit half of each ADI step takes the neighboring lines at the
  old time level, as the scheme's equations (30.20, 30.22) have it. The
  book's worked Example 30.5 instead takes each neighbor as just updated,
  sweep by sweep; `:updated-neighbors? true` does that, to reproduce it.
  Both converge to the same steady state."
  [{:keys [k nx ny dx dt t-end initial boundary scheme updated-neighbors?] :or {initial 0.0 scheme :adi}}]
  (let [{:keys [top bottom left right]} boundary
        lam (/ (* k dt) (* dx dx))
        steps (long (math/round (/ t-end dt)))
        T0 (vec (for [j (range (inc ny))]
                  (vec (for [i (range (inc nx))]
                         (cond (= j 0) (double bottom) (= j ny) (double top)
                               (= i 0) (double left) (= i nx) (double right)
                               :else (double initial))))))
        explicit (fn [T]
                   (reduce (fn [T' [i j]]
                             (assoc-in T' [j i] (+ (get-in T [j i])
                                                   (* lam (+ (get-in T [j (inc i)]) (get-in T [j (dec i)])
                                                             (get-in T [(inc j) i]) (get-in T [(dec j) i])
                                                             (* -4.0 (get-in T [j i])))))))
                           T (for [j (range 1 ny) i (range 1 nx)] [i j])))
        ;; implicit along one axis for each line of the other:
        ;; -l T_-1 + 2(1 + l) T - l T_+1 = l T'_-1 + 2(1 - l) T + l T'_+1
        half (fn [T along-y?]
               (let [lines (if along-y? (range 1 nx) (range 1 ny))
                     m (if along-y? (dec ny) (dec nx))]
                 (reduce (fn [T' line]
                           (let [src (if updated-neighbors? T' T)
                                 get (fn [p q] (if along-y? (get-in src [q p]) (get-in src [p q])))
                                 ;; p the fixed line index, q runs along it
                                 own (fn [p q] (if along-y? (get-in T [q p]) (get-in T [p q])))
                                 rhs (mapv (fn [q] (+ (* lam (get (dec line) q)) (* 2.0 (- 1.0 lam) (own line q)) (* lam (get (inc line) q))))
                                           (range 1 (inc m)))
                                 rhs (-> rhs (update 0 + (* lam (own line 0))) (update (dec m) + (* lam (own line (inc m)))))
                                 sol (tridiagonal (vec (repeat m (- lam))) (vec (repeat m (* 2.0 (+ 1.0 lam))))
                                                  (vec (repeat m (- lam))) rhs)]
                             (reduce (fn [T' q] (if along-y? (assoc-in T' [(inc q) line] (sol q)) (assoc-in T' [line (inc q)] (sol q))))
                                     T' (range m))))
                         T lines)))
        step (if (= scheme :explicit) explicit (fn [T] (half (half T true) false)))]
    (vec (reductions (fn [[t T] _] [(+ t dt) (step T)]) [0.0 T0] (range steps)))))
