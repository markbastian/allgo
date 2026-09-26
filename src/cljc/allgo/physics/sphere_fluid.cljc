(ns allgo.physics.sphere-fluid
  "Incompressible flow on the surface of a sphere.

  A companion to `allgo.physics.fluid`, which solves the same physics on a
  flat grid, and a deliberately different answer to it. That one tracks
  velocity and enforces incompressibility by solving for a pressure that
  cancels the divergence. This one never represents divergence at all.

  ## Why vorticity, and not velocity and pressure

  A sphere is a closed surface of genus zero. By the Hodge decomposition,
  a tangent field on it splits into a divergence-free part, a curl-free
  part, and a harmonic part -- and on a sphere the harmonic part is empty,
  because there is no loop you cannot shrink to a point. So *every*
  divergence-free field on a sphere is the rotated gradient of some scalar
  streamfunction, with nothing left over:

      u = grad(psi) x n        exactly divergence-free, by construction
      omega = laplacian(psi)   the vorticity that goes with it

  That turns the whole simulation into one scalar field. Incompressibility
  stops being a constraint to enforce and becomes a property of the
  representation: there is no pressure, no projection, and no residual --
  the divergence is zero to round-off rather than to solver tolerance.
  The checcurboard that forces a flat solver onto a staggered grid cannot
  form either, and it is worth being exact about why: that mode lives in
  the pressure, and there is no pressure. Centered differences do have an
  odd-even null space of their own -- an alternating streamfunction
  differences to nothing -- but nothing here drives it, because the
  streamfunction is not iterated toward; it is solved for directly, from
  a Laplacian whose checcurboard eigenvalue is the largest it has rather
  than zero.

  What it costs: this is two-dimensional incompressible flow and nothing
  else. No free surface, no obstacles, no compressibility. On a sphere
  that is the interesting case anyway -- it is the barotropic vorticity
  equation, which is the simplest model of a planetary atmosphere.

  ## Where this sits

  Yang, Corse, Lu, Wolper and Jiang, *Real-Time Fluid Simulation on the
  Surface of a Sphere* (2019), solve the same problem the other way:
  velocity and pressure on a staggered spherical grid, with an FFT
  pressure solve and a careful treatment of the pole, where the velocity
  there is recovered by averaging the tangent vectors around it. Their
  contribution is largely that treatment -- it removes the spectral
  filtering earlier work needed to keep the pole quiet.

  This takes the other road, the one the more recent surface-flow work
  has converged on: carry vorticity, not velocity. Elcott and others'
  circulation-preserving simplicial fluids, and the streamfunction and
  vorticity formulations for general surfaces that followed, all rest on
  the same observation -- that incompressibility is cheaper to *represent*
  than to enforce. On a sphere the payoff is unusually complete, because
  the harmonic part of the decomposition is empty and the streamfunction
  is therefore the whole story.

  The trade is real and worth stating. Velocity and pressure generalize:
  to three dimensions, to free surfaces, to obstacles, none of which this
  can do. Vorticity does not generalize, and in exchange it is exact.

  ## The equation

      D(omega + f)/Dt = nu * laplacian(omega)

  `f = 2*rotation*sin(latitude)` is the planetary vorticity: the vorticity
  the fluid already has by sitting on a spinning ball. Advecting the
  *absolute* vorticity `omega + f` rather than `omega` alone is the whole
  of the Coriolis effect, and it is what organizes the flow into zonal
  bands -- a fluid parcel pushed north must lose relative vorticity to
  keep the sum, which turns it back. Set `:rotation` to zero and the
  bands disappear into ordinary two-dimensional turbulence.

  ## The grid, and the poles

  Latitude and longitude, with latitudes at cell *centers*:

      phi_j = -pi/2 + (j + 1/2) * dphi     j in [0, nlat)
      lam_i = i * dlam                     i in [0, nlon), periodic

  The half-cell offset matters: no sample sits on a pole, so `cos(phi)`
  never vanishes and nothing is ever divided by zero. The poles end up on
  cell *faces*, where `cos` is zero, and a zero at a face is exactly
  right -- it is the statement that no flux crosses a point.

  The classical pole problem is then almost entirely absent, and it is
  worth being clear about why, because it is a consequence of the
  formulation rather than a trick. The singularity in these coordinates is
  a singularity of the *frame*, not of the physics: a vector has no
  well-defined east at the pole, but a scalar has a perfectly good value
  there. This solver advects only scalars, and does it by rotating points
  along great circles in 3D, where no frame appears at all. What is left
  is interpolation near the pole, handled by the one rule that a point
  just past the pole at longitude `lam` is a point just before it at
  `lam + pi`.

  `allgo.procedural.planet` argues at length against latitude and
  longitude, and this is the case it does not cover: there the grid was a
  lossy projection of a function that had no seams, while here the grid
  *is* the state. Nothing is projected, and the axis buys a fast direct
  solver that no other arrangement offers.

  ## Solving

  Each step:

    1. advect the vorticity, semi-Lagrangian, along great circles
    2. diffuse, implicitly, if there is any viscosity
    3. solve for the streamfunction, with the Coriolis term folded in
    4. read the velocity off `psi`

  Steps 2 and 3 are the *same solve* with a different diagonal, and all of
  them are direct. The operator is separable and the longitude axis is
  periodic, so a Fourier transform along it diagonalizes that half of the
  Laplacian and leaves one tridiagonal system per zonal wavenumber.
  Transform, solve, transform back: `O(n log n)`, exact to round-off, no
  iteration and no tolerance. It is the classical fast spherical Poisson
  solver, and it is the reason the axis is worth keeping.

  Three operators go through that one piece of machinery:

      laplacian                       Poisson, for the streamfunction
      laplacian - alpha               Helmholtz, for implicit diffusion
      laplacian + c d/dlam            the semi-implicit Coriolis term

  The third is the one that lets the planet spin. Rotation makes the
  vorticity-streamfunction loop oscillate -- those are Rossby waves -- and
  the fastest of them are the largest-scale ones, with frequency near the
  rotation rate itself. Any scheme that steps them from the start of the
  interval is doing forward Euler on an oscillation, which grows at every
  step size there is; measured here, energy conserved to one percent at a
  rotation of 20 and ten thousand times too large at 80, at every step
  size tried. Solving the term implicitly instead removes the restriction,
  and it costs nothing extra because `d/dlam` is diagonal in the same
  Fourier basis -- the tridiagonal system merely becomes complex. See
  `solve-rossby!`."
  (:require [allgo.array :as a]
            [allgo.numerics.fft :as fft]
            [allgo.numerics.tridiagonal :as tri]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; The grid

(defn grid
  "The geometry, with everything that depends only on it worked out once.

  `:nlon` must be a power of two -- the longitude transform is radix two.
  Cell area is `R^2 cos(phi) dphi dlam`, which is what every integral over
  the sphere below weights by."
  ([] (grid {}))
  ([{:keys [nlat nlon radius] :or {nlat 64 nlon 128 radius 1.0}}]
   (let [nlat (long nlat) nlon (long nlon) radius (double radius)]
     (when-not (fft/power-of-two? nlon)
       (throw (ex-info "nlon must be a power of two" {:nlon nlon})))
     (let [dphi (/ math/PI nlat)
           dlam (/ (* 2.0 math/PI) nlon)
           phi (a/f64 nlat)
           cos-phi (a/f64 nlat)
           sin-phi (a/f64 nlat)
           ;; cos(latitude) at the cell faces, so face 0 and face nlat are
           ;; the poles and are exactly zero.
           ^doubles cos-face (a/f64 (inc nlat))
           area (a/f64 nlat)
           ;; Longitude trig depends only on the column, and advection
           ;; wants it for every cell of every row on every step.
           cos-lam (a/f64 nlon)
           sin-lam (a/f64 nlon)]
       (dotimes [i nlon]
         (aset cos-lam i (math/cos (* i dlam)))
         (aset sin-lam i (math/sin (* i dlam))))
       (dotimes [j nlat]
         (let [p (+ (- (* 0.5 math/PI)) (* (+ j 0.5) dphi))]
           (aset phi j p)
           (aset cos-phi j (math/cos p))
           (aset sin-phi j (math/sin p))
           (aset area j (* radius radius (math/cos p) dphi dlam))))
       (dotimes [j (inc nlat)]
         (aset cos-face j (math/cos (+ (- (* 0.5 math/PI)) (* j dphi)))))
       ;; The two ends are the poles; make them exactly zero rather than
       ;; whatever cos returns for an argument near pi/2.
       (aset cos-face 0 0.0)
       (aset cos-face nlat 0.0)
       {:nlat nlat :nlon nlon :radius radius
        :dphi dphi :dlam dlam
        :phi phi :cos-phi cos-phi :sin-phi sin-phi
        :cos-face cos-face :area area
        :cos-lam cos-lam :sin-lam sin-lam
        :cells (* nlat nlon)}))))

(defn idx
  "Flat index of cell `[i j]` -- longitude `i`, latitude `j`."
  ^long [{:keys [nlon]} i j]
  (+ (* (long j) (long nlon)) (long i)))

(defn direction
  "The unit vector of cell `[i j]`, in the world frame the rest of the
  library uses: `y` toward the north pole."
  [{:keys [^doubles cos-phi ^doubles sin-phi ^double dlam]} i j]
  (let [j (long j)
        lam (* (long i) dlam)
        c (aget cos-phi j)]
    [(* c (math/cos lam)) (aget sin-phi j) (* c (math/sin lam))]))

;; ---------------------------------------------------------------------------
;; The separable solve
;;
;; Fourier along longitude, tridiagonal along latitude. Both the Poisson
;; solve for the streamfunction and the implicit diffusion step are the
;; same operator with a different diagonal, so they are the same code.

(defn- copy-into!
  "`dst <- src`. Not `System/arraycopy`, which is JVM only, and this
  namespace has to run in a browser."
  [^doubles dst ^doubles src]
  (dotimes [i (alength dst)] (aset dst i (aget src i)))
  dst)

(defn- workspace
  "Scratch for one solve, sized to the grid and reused across steps."
  [{:keys [nlat nlon]}]
  (let [nlat (long nlat) nlon (long nlon)
        cells (* nlat nlon)]
    {:re (a/f64 cells)
     :im (a/f64 cells)
     :row-re (a/f64 nlon)
     :row-im (a/f64 nlon)
     :sub (a/f64 nlat)
     :diag (a/f64 nlat)
     :super (a/f64 nlat)
     :rhs (a/f64 nlat)
     :sol (a/f64 nlat)
     :scratch (a/f64 nlat)
     ;; The imaginary halves, for the semi-implicit solve.
     :diag-im (a/f64 nlat)
     :rhs-im (a/f64 nlat)
     :sol-im (a/f64 nlat)
     :scratch-im (a/f64 nlat)}))

(defn solve-helmholtz!
  "Solves `(laplacian - alpha) x = b` on the sphere, into `out`.

  `alpha` zero gives Poisson's equation, and any positive `alpha` gives
  the implicit diffusion step -- which is the same system with a heavier
  diagonal, and is nonsingular where Poisson's is not.

  At `alpha` zero the operator has the constant functions in its null
  space, so the solution is fixed only up to a constant and `b` must
  integrate to zero over the sphere. Neither matters here: the caller
  removes the mean of `b` first, and a constant added to a streamfunction
  changes no velocity."
  [grid ws ^doubles b alpha ^doubles out]
  (let [alpha (double alpha)
        {:keys [nlat nlon ^double radius ^double dphi ^double dlam
                ^doubles cos-phi ^doubles cos-face]} grid
        {:keys [^doubles re ^doubles im ^doubles row-re ^doubles row-im
                ^doubles sub ^doubles diag ^doubles super
                ^doubles rhs ^doubles sol ^doubles scratch]} ws
        nlat (long nlat) nlon (long nlon)
        r2 (* radius radius)
        dphi2 (* dphi dphi)]
    ;; 1. Transform every latitude row along longitude. A row is
    ;;    contiguous, which is why the flat index runs longitude fastest.
    (dotimes [j nlat]
      (let [base (* j nlon)]
        (dotimes [i nlon]
          (aset row-re i (aget b (+ base i)))
          (aset row-im i 0.0))
        (fft/forward! row-re row-im)
        (dotimes [i nlon]
          (aset re (+ base i) (aget row-re i))
          (aset im (+ base i) (aget row-im i)))))
    ;; 2. One tridiagonal system per zonal wavenumber, solved twice --
    ;;    once for the real part and once for the imaginary. The matrix is
    ;;    real, so the two halves never mix and complex arithmetic is not
    ;;    needed anywhere.
    (dotimes [m nlon]
      (let [;; The second difference around a ring of nlon points has
            ;; eigenvalue -(2 - 2cos(2 pi m / nlon)) / dlam^2.
            k-m (/ (- 2.0 (* 2.0 (math/cos (/ (* 2.0 math/PI m) nlon))))
                   (* dlam dlam))
            singular? (and (zero? m) (zero? alpha))]
        (dotimes [j nlat]
          (let [cj (aget cos-phi j)
                cm (aget cos-face j)
                cp (aget cos-face (inc j))]
            (aset sub j cm)
            (aset super j cp)
            (aset diag j (- (- (+ cm cp))
                            (* k-m (/ dphi2 cj))
                            (* alpha r2 cj dphi2)))))
        (when singular?
          ;; Pin the last unknown. The equations are linearly dependent
          ;; here -- they sum to zero, which is why the mean of b has to
          ;; vanish -- so dropping one and fixing the constant is not an
          ;; approximation, it is choosing which solution to return.
          (aset sub (dec nlat) 0.0)
          (aset diag (dec nlat) 1.0))
        (dotimes [part 2]
          (let [^doubles src (if (zero? part) re im)]
            (dotimes [j nlat]
              (aset rhs j (* r2 (aget cos-phi j) dphi2 (aget src (+ (* j nlon) m)))))
            (when singular?
              (aset rhs (dec nlat) 0.0))
            (if (tri/solve sub diag super rhs sol scratch)
              (dotimes [j nlat]
                (aset src (+ (* j nlon) m) (aget sol j)))
              (throw (ex-info "tridiagonal solve failed" {:wavenumber m :alpha alpha})))))))
    ;; 3. Back to longitude.
    (dotimes [j nlat]
      (let [base (* j nlon)]
        (dotimes [i nlon]
          (aset row-re i (aget re (+ base i)))
          (aset row-im i (aget im (+ base i))))
        (fft/inverse! row-re row-im)
        (dotimes [i nlon]
          (aset out (+ base i) (aget row-re i)))))
    out))

;; ---------------------------------------------------------------------------
;; Integrals over the sphere

(defn d-lambda!
  "Centered longitude derivative of a scalar, into `out`. Periodic, so no
  boundary and no pole ever enters it."
  [{:keys [nlat nlon ^double dlam]} ^doubles f ^doubles out]
  (let [nlat (long nlat) nlon (long nlon)]
    (dotimes [j nlat]
      (let [base (* j nlon)]
        (dotimes [i nlon]
          (aset out (+ base i)
                (/ (- (aget f (+ base (mod (inc i) nlon)))
                      (aget f (+ base (mod (+ i nlon -1) nlon))))
                   (* 2.0 dlam))))))
    out))

(defn solve-rossby!
  "Solves `(laplacian + c * d/dlam) psi = b`, the semi-implicit operator.

  This is what lets the planet spin fast. The Coriolis term makes the
  vorticity, streamfunction and velocity loop oscillate -- those are
  Rossby waves -- and the fastest of them are the *largest* scale ones,
  with frequency near the rotation rate itself. Any scheme that steps them
  with values from the start of the interval is doing forward Euler on an
  oscillation, which grows without bound however small the step. Centring
  the trajectory reduces the growth from `O((sigma dt)^2)` per step to
  `O((sigma dt)^4)`, and that is still hopeless once `sigma dt` reaches
  one, which at a rotation of 80 it does.

  Treating the term implicitly removes the restriction outright, and the
  trapezoidal rule is the implicit rule to use. The instinct that implicit
  means stable is right, but the usual implicit rule is not: backward
  Euler is stable *because* it damps, and what it would damp here are
  precisely the waves that organize the flow into bands. The trapezoidal
  rule has an amplification factor of exactly one at every step size --
  neutral, not damped. What it gets wrong at a large step is the phase
  speed of the fast waves, not their existence.

  The cost is almost nothing, which is a happy accident of this solver's
  structure. `d/dlam` is diagonal in the same Fourier basis that already
  diagonalizes the longitude half of the Laplacian, so the operator is
  still one tridiagonal system per zonal wavenumber -- with an imaginary
  number added to its diagonal, and therefore a complex solve rather than
  a real one. Same algorithm, same order of work.

  (The diffusion in this namespace has been implicit from the start, for
  the opposite reason: diffusion really is dissipative, so backward Euler
  is exactly right for it. See `solve-helmholtz!`.)"
  [grid ws ^doubles b c ^doubles out]
  (let [c (double c)
        {:keys [nlat nlon ^double radius ^double dphi ^double dlam
                ^doubles cos-phi ^doubles cos-face]} grid
        {:keys [^doubles re ^doubles im ^doubles row-re ^doubles row-im
                ^doubles sub ^doubles diag ^doubles super
                ^doubles rhs ^doubles sol ^doubles scratch
                ^doubles diag-im ^doubles rhs-im ^doubles sol-im
                ^doubles scratch-im]} ws
        nlat (long nlat) nlon (long nlon)
        r2 (* radius radius)
        dphi2 (* dphi dphi)]
    (dotimes [j nlat]
      (let [base (* j nlon)]
        (dotimes [i nlon]
          (aset row-re i (aget b (+ base i)))
          (aset row-im i 0.0))
        (fft/forward! row-re row-im)
        (dotimes [i nlon]
          (aset re (+ base i) (aget row-re i))
          (aset im (+ base i) (aget row-im i)))))
    (dotimes [m nlon]
      (let [k-m (/ (- 2.0 (* 2.0 (math/cos (/ (* 2.0 math/PI m) nlon))))
                   (* dlam dlam))
            ;; A centered first difference around the ring transforms to
            ;; `i sin(2 pi m / n) / dlam` -- pure imaginary, and diagonal.
            dl (/ (math/sin (/ (* 2.0 math/PI m) nlon)) dlam)
            singular? (zero? m)]
        (dotimes [j nlat]
          (let [cj (aget cos-phi j)
                cm (aget cos-face j)
                cp (aget cos-face (inc j))]
            (aset sub j cm)
            (aset super j cp)
            (aset diag j (- (- (+ cm cp)) (* k-m (/ dphi2 cj))))
            (aset diag-im j (* c dl r2 cj dphi2))))
        (when singular?
          (aset sub (dec nlat) 0.0)
          (aset diag (dec nlat) 1.0)
          (aset diag-im (dec nlat) 0.0))
        (dotimes [j nlat]
          (let [w (* r2 (aget cos-phi j) dphi2)]
            (aset rhs j (* w (aget re (+ (* j nlon) m))))
            (aset rhs-im j (* w (aget im (+ (* j nlon) m))))))
        (when singular?
          (aset rhs (dec nlat) 0.0)
          (aset rhs-im (dec nlat) 0.0))
        (if (tri/solve-complex sub diag diag-im super rhs rhs-im sol sol-im
                               scratch scratch-im)
          (dotimes [j nlat]
            (aset re (+ (* j nlon) m) (aget sol j))
            (aset im (+ (* j nlon) m) (aget sol-im j)))
          (throw (ex-info "semi-implicit solve failed" {:wavenumber m :c c})))))
    (dotimes [j nlat]
      (let [base (* j nlon)]
        (dotimes [i nlon]
          (aset row-re i (aget re (+ base i)))
          (aset row-im i (aget im (+ base i))))
        (fft/inverse! row-re row-im)
        (dotimes [i nlon]
          (aset out (+ base i) (aget row-re i)))))
    out))

(defn integral
  "Area-weighted integral of a field over the sphere."
  ^double [{:keys [nlat nlon ^doubles area]} ^doubles f]
  (let [nlat (long nlat) nlon (long nlon)]
    (loop [j 0 acc 0.0]
      (if (= j nlat)
        acc
        (let [a (aget area j)
              base (* j nlon)]
          (recur (inc j)
                 (+ acc (* a (loop [i 0 s 0.0]
                               (if (= i nlon) s (recur (inc i) (+ s (aget f (+ base i))))))))))))))

(defn- remove-mean!
  "Subtracts the area-weighted mean, so that the field integrates to zero.

  Vorticity is a curl, so its integral over a closed surface is zero as a
  matter of fact and not of taste; what accumulates here is the drift of
  repeated interpolation. Removing it is also exactly the condition the
  Poisson solve needs to be solvable."
  [{:keys [^double radius] :as grid} ^doubles f]
  (let [total (integral grid f)
        mean (/ total (* 4.0 math/PI radius radius))]
    (dotimes [i (alength f)] (aset f i (- (aget f i) mean)))
    f))

;; ---------------------------------------------------------------------------
;; Velocity from the streamfunction

(defn- pole-wrapped
  "The value of a scalar field at latitude row `j`, allowing one row past
  either pole.

  Going over a pole and continuing brings you back on the far side, so the
  row just beyond the last one is the last one again, half a world along
  in longitude. Exact for a scalar, and the only thing either pole needs."
  ;; Five arguments, so no primitive hints: Clojure's primitive function
  ;; interfaces stop at four.
  [^doubles f nlat nlon i j]
  (let [nlat (long nlat) nlon (long nlon) j (long j)
        i (mod (long i) nlon)]
    (cond
      (neg? j) (aget f (mod (+ i (quot nlon 2)) nlon))
      (>= j nlat) (aget f (+ (* (dec nlat) nlon) (mod (+ i (quot nlon 2)) nlon)))
      :else (aget f (+ (* j nlon) i)))))

(defn velocity!
  "Reads the velocity off the streamfunction, into `east` and `north`.

      u_east  = -(1/R) d(psi)/d(phi)
      u_north =  (1/(R cos phi)) d(psi)/d(lam)

  which is `grad(psi) x n` written out. Centered differences, with the
  pole rule standing in for the rows that do not exist."
  [{:keys [nlat nlon ^double radius ^double dphi ^double dlam ^doubles cos-phi]}
   ^doubles psi ^doubles east ^doubles north]
  (let [nlat (long nlat) nlon (long nlon)]
    (dotimes [j nlat]
      (let [cj (aget cos-phi j)
            base (* j nlon)]
        (dotimes [i nlon]
          (let [dpsi-dphi (/ (- (pole-wrapped psi nlat nlon i (inc j))
                                (pole-wrapped psi nlat nlon i (dec j)))
                             (* 2.0 dphi))
                dpsi-dlam (/ (- (pole-wrapped psi nlat nlon (inc i) j)
                                (pole-wrapped psi nlat nlon (dec i) j))
                             (* 2.0 dlam))]
            (aset east (+ base i) (/ (- dpsi-dphi) radius))
            (aset north (+ base i) (/ dpsi-dlam (* radius cj)))))))
    [east north]))

;; ---------------------------------------------------------------------------
;; Advection

(defn- bilinear
  "Bilinear value at a latitude and longitude, and -- when `lo` and `hi`
  are given -- the smallest and largest of the four stencil values,
  written at index `k`.

  Those bounds are what a corrected advection has to be clamped to. An
  interpolation can never leave the range of what it interpolated, which
  is why plain semi-Lagrangian advection cannot go unstable; a correction
  added on top can, unless it is held inside the same range."
  [{:keys [nlat nlon dphi dlam]} ^doubles f phi lam ^doubles lo ^doubles hi k]
  (let [nlat (long nlat) nlon (long nlon)
        dphi (double dphi) dlam (double dlam)
        phi (double phi) lam (double lam)
        two-pi (* 2.0 math/PI)
        lam (let [l (rem lam two-pi)] (if (neg? l) (+ l two-pi) l))
        fx (/ lam dlam)
        fy (- (/ (+ phi (* 0.5 math/PI)) dphi) 0.5)
        i0 (long (math/floor fx))
        j0 (long (math/floor fy))
        tx (- fx i0)
        ty (- fy j0)
        ;; `pole-wrapped` has to give up its primitive return hint to fit
        ;; inside the four-argument limit, so pin the values here instead.
        v00 (double (pole-wrapped f nlat nlon i0 j0))
        v10 (double (pole-wrapped f nlat nlon (inc i0) j0))
        v01 (double (pole-wrapped f nlat nlon i0 (inc j0)))
        v11 (double (pole-wrapped f nlat nlon (inc i0) (inc j0)))]
    (when lo
      (let [k (long k)]
        (aset lo k (min v00 (min v10 (min v01 v11))))
        (aset hi k (max v00 (max v10 (max v01 v11))))))
    (+ (* (- 1.0 ty) (+ (* (- 1.0 tx) v00) (* tx v10)))
       (* ty (+ (* (- 1.0 tx) v01) (* tx v11))))))

(defn sample
  "Bilinear value of a scalar field at latitude `phi`, longitude `lam`.

  Periodic in longitude and pole-wrapped in latitude, so any direction on
  the sphere is in range and no case needs guarding by the caller."
  ^double [grid ^doubles f ^double phi ^double lam]
  (double (bilinear grid f phi lam nil nil 0)))

(defn departures!
  "For each cell, where the material arriving there came from.

  The backtrace, and the only part of advection that touches geometry.
  It depends on the velocity and not on what is being carried, so a step
  that advects several fields computes this once and samples many times --
  which matters, because this is where all the trigonometry is.

  On a sphere the departure point is exact rather than approximate. The
  velocity is tangent, so the path is a great circle, and stepping back
  along one is a rotation: a point `p` moved back by arc angle `theta` in
  direction `u` lands at `p cos(theta) - u sin(theta)`. No coordinate
  frame appears anywhere in that, which is precisely why the poles give no
  trouble -- there is no frame there to be singular.

  A negative `dt` traces forward instead, which is what the corrected
  scheme needs."
  [{:keys [nlat nlon ^double radius ^doubles cos-phi ^doubles sin-phi
           ^doubles cos-lam ^doubles sin-lam]}
   ^doubles east ^doubles north dt ^doubles dep-phi ^doubles dep-lam]
  (let [dt (double dt)
        nlat (long nlat) nlon (long nlon)]
    (dotimes [j nlat]
      (let [base (* j nlon)
            sp (aget sin-phi j)
            cp (aget cos-phi j)]
        (dotimes [i nlon]
          (let [k (+ base i)
                clam (aget cos-lam i)
                slam (aget sin-lam i)
                ;; The east and north unit vectors here, written out
                ;; rather than called for: this is the innermost loop of
                ;; the whole solver.
                ue (aget east k)
                un (aget north k)
                vx (+ (* ue (- slam)) (* un (- sp) clam))
                vy (* un cp)
                vz (+ (* ue clam) (* un (- sp) slam))
                speed (math/sqrt (+ (* vx vx) (* vy vy) (* vz vz)))]
            (if (< speed 1e-14)
              (do (aset dep-phi k (math/asin (max -1.0 (min 1.0 sp))))
                  (aset dep-lam k (math/atan2 (* cp slam) (* cp clam))))
              (let [theta (/ (* speed dt) radius)
                    c (math/cos theta)
                    sc (/ (math/sin theta) speed)
                    qx (- (* cp clam c) (* vx sc))
                    qy (- (* sp c) (* vy sc))
                    qz (- (* cp slam c) (* vz sc))
                    len (math/sqrt (+ (* qx qx) (* qy qy) (* qz qz)))]
                (aset dep-phi k (math/asin (max -1.0 (min 1.0 (/ qy len)))))
                (aset dep-lam k (math/atan2 qz qx))))))))
    dep-phi))

(defn sample-at!
  "Interpolates a field at the departure points, into `out`.

  With `lo` and `hi`, also records the range of each stencil, which a
  corrected scheme needs to clamp itself to."
  ([grid field dep-phi dep-lam out] (sample-at! grid field dep-phi dep-lam out nil nil))
  ([grid ^doubles field ^doubles dep-phi ^doubles dep-lam
    ^doubles out ^doubles lo ^doubles hi]
   (dotimes [k (alength out)]
     (aset out k (double (bilinear grid field (aget dep-phi k) (aget dep-lam k) lo hi k))))
   out))

(defn advect!
  "Semi-Lagrangian advection of a scalar, into `out`.

  Stam's stable fluids: for each cell ask where the material arriving
  there came from, and interpolate the field at that point. It cannot go
  unstable however large the step, because the answer is always an
  interpolation between values that already exist -- and for the same
  reason it cannot create a new maximum, which is what keeps a tracer
  inside the range it was painted with.

  It allocates its own departure arrays, so it is for occasional use;
  `step!` computes them once and shares them across fields."
  [{:keys [cells] :as grid} ^doubles field ^doubles east ^doubles north dt ^doubles out]
  (let [dep-phi (a/f64 (long cells))
        dep-lam (a/f64 (long cells))]
    (departures! grid east north dt dep-phi dep-lam)
    (sample-at! grid field dep-phi dep-lam out)))

(defn correct!
  "The MacCormack correction, given a forward and a round-trip pass.

  Plain semi-Lagrangian advection is stable and first order, and it pays
  for the stability in dissipation -- every step interpolates, every
  interpolation averages, and a vortex quietly turns into a smudge. The
  fix is to measure that error rather than accept it: advect forward, then
  advect the result *backward*, and whatever fails to come back where it
  started is an estimate of what the forward pass lost. Subtract half of
  it and the scheme is second order.

  `forward` is the plain advected field, `round-trip` is that field traced
  back again, and `field` is what was started with. The correction is an
  extrapolation rather than an interpolation and so can overshoot, which
  is why it is clamped to the range of the stencil the forward pass
  actually read: inside those bounds the result is still bounded by values
  that existed, which is the property the whole method rests on. Where the
  clamp bites, it falls back to the plain answer.

  ## One honest caveat, and when it bites

  This is no longer strictly dissipative. Energy can come out above where
  it went in, which a first-order scheme can never do, and the clamp does
  not prevent it: bounding every value by its own stencil bounds the
  *field*, not its gradients, and the energy lives in the gradients.

  In smooth flow it does not matter -- a banded rotating simulation comes
  back within a few percent. It matters in fully developed turbulence,
  where structures reach the grid scale and the correction pushes against
  the clamp everywhere at once. Measured on jets well below their
  stability threshold, over fourteen time units: energy ratio 0.45 with
  the plain scheme and 2.87 with this one, and *worse* with resolution --
  0.72, 2.87, 4.55 at 32, 64 and 128 latitudes. Growing with resolution
  is the tell that it is the scheme and not the physics.

  The fix is the one every two-dimensional turbulence code uses, and it is
  not a workaround: a flow with no dissipation at all is not a well posed
  numerical problem, because the energy that cascades to the grid scale
  has nowhere to go. `:viscosity` gives it somewhere. On the same turbulent
  case, a viscosity of 5e-4 takes the ratio from 2.87 to 0.56 at 64
  latitudes and from 4.55 to 0.67 at 128 -- and, more to the point, stops
  it growing with resolution.

  There is no one value that suits both regimes, so this is a knob and not
  a default. That same 5e-4 takes the *banded* case from 0.91 down to
  0.39, which is no longer a simulation of anything. Banded flow wants
  zero or 1e-4; turbulent flow wants five times that. Which regime you are
  in is decided by the rotation rate -- see `zonal-jets!`."
  [^doubles field ^doubles forward ^doubles round-trip ^doubles lo ^doubles hi ^doubles out]
  (dotimes [k (alength field)]
    (let [corrected (+ (aget forward k) (* 0.5 (- (aget field k) (aget round-trip k))))]
      (aset out k (min (aget hi k) (max (aget lo k) corrected)))))
  out)

;; ---------------------------------------------------------------------------
;; The simulation

(defn simulation
  "A flow on a sphere, ready to step.

  * `:grid` -- from `grid`, or its options.
  * `:rotation` -- the sphere's angular rate, which sets the Coriolis
    term. Zero for plain two-dimensional turbulence.
  * `:viscosity` -- diffuses vorticity. Zero leaves it inviscid, in which
    case the only dissipation is whatever the interpolation loses.
  * `:tracer?` -- carry a passive scalar along with the flow, for seeing
    what the velocity is doing.
  * `:advection` -- `:corrected` (the default) for MacCormack, or
    `:semi-lagrangian` for the plain first-order trace. The plain one is
    half the cost and loses a great deal more; it is kept because the
    difference between the two is the clearest thing in this namespace to
    look at."
  ([] (simulation {}))
  ([{:keys [rotation viscosity tracer? advection]
     :or {rotation 0.0 viscosity 0.0 tracer? true advection :corrected}
     :as opts}]
   ;; `:grid` is read from opts rather than destructured, so that the
   ;; binding does not shadow the var of the same name.
   (let [spec (:grid opts)
         ;; A built grid carries :cells and an options map does not, which
         ;; is the only reliable way to tell them apart -- both have :nlat.
         g (if (:cells spec) spec (grid (or spec {})))
         cells (long (:cells g))]
     {:grid g
      :rotation (double rotation)
      :viscosity (double viscosity)
      :vorticity (a/f64 cells)
      :streamfunction (a/f64 cells)
      :east (a/f64 cells)
      :north (a/f64 cells)
      :tracer (when tracer? (a/f64 cells))
      :advection advection
      :coriolis (:coriolis opts :semi-implicit)
      :trajectory (:trajectory opts :centered)
      :scratch-a (a/f64 cells)
      :scratch-b (a/f64 cells)
      :scratch-c (a/f64 cells)
      :scratch-lo (a/f64 cells)
      :scratch-hi (a/f64 cells)
      :dep-phi (a/f64 cells)
      :dep-lam (a/f64 cells)
      :fwd-phi (a/f64 cells)
      :fwd-lam (a/f64 cells)
      :saved-vorticity (a/f64 cells)
      :saved-tracer (when tracer? (a/f64 cells))
      :saved-east (a/f64 cells)
      :saved-north (a/f64 cells)
      :mid-east (a/f64 cells)
      :mid-north (a/f64 cells)
      :psi-tilde (a/f64 cells)
      :dlam-a (a/f64 cells)
      :dlam-b (a/f64 cells)
      :workspace (workspace g)
      :time 0.0
      :opts opts})))

(defn planetary-vorticity
  "`f = 2 * rotation * sin(latitude)`: the vorticity the fluid has merely
  by sitting on a spinning sphere."
  ^double [^double rotation ^double sin-phi]
  (* 2.0 rotation sin-phi))

(defn- beta-source!
  "Applies the planetary vorticity gradient as a source term, in place.

      D(omega)/Dt = -v * beta,   beta = (1/R) df/dphi = 2 * rotation * cos(phi) / R

  The alternative to carrying `omega + f` through the advection, and much
  the better conditioned of the two. `f` is the vorticity the fluid has by
  sitting on a spinning ball, and on a fast planet it dwarfs the flow's
  own -- at a rotation of 80 it runs to +/-160 where the vorticity of
  interest is +/-8. Advecting the sum then means interpolating a field
  twenty times larger than the answer and subtracting almost all of it
  back off, so the interpolation error on the part that does not matter
  swamps the part that does. The error grows with rotation, which is the
  tell.

  Here nothing large is ever interpolated. `f` never enters the advection;
  only its gradient enters, as a source, and the source is the size of the
  effect rather than the size of the background."
  [{:keys [nlat nlon ^double radius ^doubles cos-phi]}
   ^doubles vorticity ^doubles north rotation dt]
  (let [nlat (long nlat) nlon (long nlon)
        rotation (double rotation)
        dt (double dt)]
    (dotimes [j nlat]
      (let [beta (/ (* 2.0 rotation (aget cos-phi j)) radius)
            base (* j nlon)]
        (dotimes [i nlon]
          (let [k (+ base i)]
            (aset vorticity k (- (aget vorticity k) (* dt beta (aget north k))))))))
    vorticity))

(defn- add-planetary!
  "Adds or removes the planetary vorticity, in place."
  [{:keys [nlat nlon ^doubles sin-phi]} ^doubles f ^double rotation ^double sign]
  (let [nlat (long nlat) nlon (long nlon)]
    (dotimes [j nlat]
      (let [fj (* sign (planetary-vorticity rotation (aget sin-phi j)))
            base (* j nlon)]
        (dotimes [i nlon]
          (aset f (+ base i) (+ (aget f (+ base i)) fj)))))
    f))

(defn update-velocity!
  "Solves for the streamfunction and reads the velocity off it."
  [{:keys [grid workspace ^doubles vorticity ^doubles streamfunction
           ^doubles east ^doubles north] :as sim}]
  (remove-mean! grid vorticity)
  (solve-helmholtz! grid workspace vorticity 0.0 streamfunction)
  (velocity! grid streamfunction east north)
  sim)

(defn- carrier
  "Computes the departure points for one step and returns a function that
  advects a field along them, in place.

  The backtrace depends on the velocity and not on what is being carried,
  so a step that moves several fields pays for it once."
  [{:keys [grid advection ^doubles scratch-a ^doubles scratch-b ^doubles scratch-c
           ^doubles scratch-lo ^doubles scratch-hi] :as sim}
   ^doubles east ^doubles north dt]
  (let [corrected? (not= advection :semi-lagrangian)]
    (departures! grid east north dt (:dep-phi sim) (:dep-lam sim))
    (when corrected?
      (departures! grid east north (- dt) (:fwd-phi sim) (:fwd-lam sim)))
    (fn [^doubles field]
      (if corrected?
        (do (sample-at! grid field (:dep-phi sim) (:dep-lam sim)
                        scratch-a scratch-lo scratch-hi)
            (sample-at! grid scratch-a (:fwd-phi sim) (:fwd-lam sim) scratch-c)
            (correct! field scratch-a scratch-c scratch-lo scratch-hi scratch-b)
            (copy-into! field scratch-b))
        (do (sample-at! grid field (:dep-phi sim) (:dep-lam sim) scratch-a)
            (copy-into! field scratch-a))))))

(defn- diffuse!
  "Implicit diffusion of the vorticity, if there is any viscosity."
  [{:keys [grid workspace ^double viscosity ^doubles vorticity ^doubles scratch-c]} dt]
  (when (pos? viscosity)
    ;; (I - nu dt L) w' = w, which is (L - 1/(nu dt)) w' = -w/(nu dt):
    ;; the Poisson solver with a heavier diagonal.
    (let [alpha (/ 1.0 (* viscosity (double dt)))]
      (dotimes [i (alength vorticity)]
        (aset scratch-c i (* (- alpha) (aget vorticity i))))
      (solve-helmholtz! grid workspace scratch-c alpha vorticity))))

(defn- advance!
  "One advance of the carried fields by a given velocity, in place.

  Split out of `step!` because the corrected trajectory below has to run
  it twice with two different velocities and the same starting state."
  [{:keys [grid coriolis ^double rotation ^doubles vorticity tracer] :as sim}
   ^doubles east ^doubles north dt]
  (let [dt (double dt)
        carry! (carrier sim east north dt)]
    ;; The Coriolis effect: a parcel keeps omega + f, so moving north
    ;; costs it relative vorticity.
    (if (or (zero? rotation) (= coriolis :beta))
      (do (carry! vorticity)
          (when-not (zero? rotation)
            (beta-source! grid vorticity north rotation dt)))
      (do (add-planetary! grid vorticity rotation 1.0)
          (carry! vorticity)
          (add-planetary! grid vorticity rotation -1.0)))
    (diffuse! sim dt)
    (when tracer (carry! tracer))
    sim))

(defn- semi-implicit-step!
  "One step with the Coriolis term solved implicitly.

  Advect, then solve the whole linear wave part in one go rather than
  stepping it. See `solve-rossby!` for why this is the difference between
  a planet that can spin and one that cannot."
  [{:keys [grid workspace ^double rotation ^doubles vorticity
           ^doubles streamfunction ^doubles east ^doubles north
           ^doubles scratch-b ^doubles psi-tilde ^doubles dlam-a ^doubles dlam-b
           tracer] :as sim}
   dt]
  (let [dt (double dt)
        radius (double (:radius grid))
        ;; The trapezoidal weight: half of dt * 2 * rotation / R^2.
        c (/ (* dt rotation) (* radius radius))
        carry! (carrier sim east north dt)
        n (alength vorticity)]
    (copy-into! psi-tilde streamfunction)
    (carry! vorticity)                  ; omega-tilde
    (carry! psi-tilde)                  ; psi-tilde, along the same paths
    (diffuse! sim dt)
    (remove-mean! grid vorticity)
    (d-lambda! grid psi-tilde dlam-a)
    (dotimes [k n]
      (aset scratch-b k (- (aget vorticity k) (* c (aget dlam-a k)))))
    (solve-rossby! grid workspace scratch-b c streamfunction)
    ;; The same relation the solve enforces, used to read the new
    ;; vorticity off without applying a Laplacian to get it.
    (d-lambda! grid streamfunction dlam-b)
    (dotimes [k n]
      (aset vorticity k (- (aget vorticity k)
                           (* c (+ (aget dlam-b k) (aget dlam-a k))))))
    (velocity! grid streamfunction east north)
    (when tracer (carry! tracer))
    sim))

(defn step!
  "One step. Advect, diffuse, solve, and read the velocity back off.

  Mutates the simulation's arrays, like `allgo.physics.fluid` and for the
  same reason: the fields are large and a step touches all of them, so a
  version that returned fresh ones would spend its time in the allocator.

  ## Why the trajectory is corrected

  The obvious scheme traces back along the velocity as it stands at the
  start of the step. For pure advection that is fine. Here it is not,
  and the reason is worth stating because it is invisible until the
  planet spins fast.

  Vorticity, the streamfunction and the velocity form a closed loop, and
  on a rotating sphere that loop *oscillates*: it carries Rossby waves,
  whose frequency grows with the rotation rate. Stepping an oscillation
  with the value at the start of the interval is forward Euler, and
  forward Euler on an oscillation has an amplification factor of
  `sqrt(1 + (sigma dt)^2)`, which is greater than one for every step size
  there is. It does not blow up quickly -- it grows like
  `exp(sigma^2 dt t / 2)` -- so at a slow rotation nothing looks wrong,
  and at a fast one the whole field explodes after a while regardless of
  how small the step is. Measured here: energy conserved to one percent
  at a rotation of 20, and ten thousand times too large at 80, at every
  step size tried.

  So the trajectory is centered instead. Advect with the velocity at hand,
  solve for where that leaves the velocity, and then throw the provisional
  answer away and advect the *original* state again with the average of
  the two. That is the trapezoidal rule on the loop, whose amplification
  error is `O((sigma dt)^4)` rather than `O((sigma dt)^2)`, and it costs
  one extra solve and one extra advection.

  None of this applies under the default `:coriolis :semi-implicit`, which
  takes the oscillation out of the explicit part altogether and so runs
  the cheap trajectory regardless; `:trajectory` is read only by the
  explicit `:coriolis :absolute` path, where it is the difference between
  usable to a rotation of about 30 and usable to about 10. Both are kept
  because the comparison between them is the clearest demonstration in
  this namespace of why the implicit treatment is worth having."
  [{:keys [trajectory ^double rotation ^doubles vorticity ^doubles east ^doubles north
           ^doubles saved-vorticity ^doubles saved-tracer ^doubles saved-east
           ^doubles saved-north ^doubles mid-east ^doubles mid-north tracer] :as sim}
   dt]
  (let [dt (double dt)]
    (if (= (:coriolis sim) :semi-implicit)
      (semi-implicit-step! sim dt)
      (if (or (= trajectory :lagged) (zero? rotation))
        (do (advance! sim east north dt)
            (update-velocity! sim))
        (do
        ;; Remember where the step started; the corrector redoes it.
          (copy-into! saved-vorticity vorticity)
          (copy-into! saved-east east)
          (copy-into! saved-north north)
          (when tracer (copy-into! saved-tracer tracer))
        ;; Predictor, purely to find out where the velocity is going.
          (advance! sim east north dt)
          (update-velocity! sim)
          (dotimes [k (alength east)]
            (aset mid-east k (* 0.5 (+ (aget saved-east k) (aget east k))))
            (aset mid-north k (* 0.5 (+ (aget saved-north k) (aget north k)))))
        ;; Corrector, from the original state along the averaged velocity.
          (copy-into! vorticity saved-vorticity)
          (when tracer (copy-into! tracer saved-tracer))
          (advance! sim mid-east mid-north dt)
          (update-velocity! sim))))
    (assoc sim :time (+ (double (:time sim)) dt))))

;; ---------------------------------------------------------------------------
;; Setting one going

(defn add-vortex!
  "Adds a Gaussian patch of vorticity centered on a direction.

  `strength` is signed -- positive spins one way, negative the other --
  and `width` is the patch's angular radius in radians. A pair of opposite
  sign is a dipole, which propagates; a single one just spins."
  [{:keys [nlat nlon ^doubles cos-phi ^doubles sin-phi ^double dlam]}
   ^doubles field [cx cy cz] strength width]
  (let [nlat (long nlat) nlon (long nlon)
        strength (double strength)
        width (double width)
        cx (double cx) cy (double cy) cz (double cz)
        len (math/sqrt (+ (* cx cx) (* cy cy) (* cz cz)))
        cx (/ cx len) cy (/ cy len) cz (/ cz len)
        inv (/ 1.0 (* 2.0 width width))]
    (dotimes [j nlat]
      (let [sp (aget sin-phi j) cp (aget cos-phi j) base (* j nlon)]
        (dotimes [i nlon]
          (let [lam (* i dlam)
                dot (+ (* cx cp (math/cos lam)) (* cy sp) (* cz cp (math/sin lam)))
                ;; Great-circle angle to the center.
                ang (math/acos (max -1.0 (min 1.0 dot)))
                k (+ base i)]
            (aset field k (+ (aget field k)
                             (* strength (math/exp (- (* ang ang inv))))))))))
    field))

(defn zonal-jets!
  "Alternating east-west jets, as a function of latitude, plus a small
  wave to break the symmetry.

  The companion to `add-vortex!`, and the initial condition for anything
  meant to look like a planetary atmosphere. Banded flow does not arise
  from scattered vortices in any reasonable time -- on a real planet the
  bands are maintained by forcing from below over geological time, and a
  decaying simulation started from turbulence just decays. So the bands
  are put in, and what the simulation then shows is what the rotation
  does *to* them: holds them together, and rolls their shear layers into
  the chains of vortices that run along every boundary.

  `n-jets` is how many alternations from pole to pole, `amplitude` the
  vorticity of the jets, and `perturbation` the size of the wavy nudge as
  a fraction of that -- a few percent is enough, and much more swamps the
  instability it is meant to trigger. `wave` is how many wavelengths of
  the nudge go round, which sets how many vortices form along each shear
  line before they start merging."
  ([grid vorticity] (zonal-jets! grid vorticity {}))
  ([{:keys [nlat nlon ^doubles phi ^double dlam]} ^doubles vorticity
    {:keys [n-jets amplitude perturbation wave]
     :or {n-jets 9 amplitude 8.0 perturbation 0.03 wave 6}}]
   (let [nlat (long nlat) nlon (long nlon)
         n-jets (double n-jets)
         amplitude (double amplitude)
         perturbation (double perturbation)
         wave (double wave)]
     (dotimes [j nlat]
       (let [p (aget phi j)
             band (* amplitude (math/sin (* n-jets p)))
             ;; The nudge rides on the jets and vanishes with them, so it
             ;; perturbs the shear rather than adding a separate flow.
             swell (* perturbation amplitude (math/cos (* n-jets p)))
             base (* j nlon)]
         (dotimes [i nlon]
           (aset vorticity (+ base i)
                 (+ band (* swell (math/sin (* wave i dlam))))))))
     vorticity)))

(defn seed-tracer!
  "Paints the tracer in latitude bands, so that the flow's stretching and
  folding is visible in it."
  [{:keys [nlat nlon]} ^doubles tracer bands]
  (let [nlat (long nlat) nlon (long nlon) bands (double bands)]
    (dotimes [j nlat]
      (let [v (if (even? (long (math/floor (* bands (/ (double j) nlat))))) 1.0 0.0)
            base (* j nlon)]
        (dotimes [i nlon] (aset tracer (+ base i) v))))
    tracer))

;; ---------------------------------------------------------------------------
;; Diagnostics

(defn total-vorticity
  "Should be zero: vorticity is a curl, and a closed surface has no edge
  for it to leak out of."
  ^double [{:keys [grid ^doubles vorticity]}]
  (integral grid vorticity))

(defn kinetic-energy
  "`1/2 integral |u|^2 dA`. Conserved by the inviscid equation, so its
  drift measures what the interpolation is costing."
  ^double [{:keys [grid ^doubles east ^doubles north]}]
  (let [{:keys [nlat nlon ^doubles area]} grid
        nlat (long nlat) nlon (long nlon)]
    (loop [j 0 acc 0.0]
      (if (= j nlat)
        (* 0.5 acc)
        (let [a (aget area j) base (* j nlon)]
          (recur (inc j)
                 (+ acc (* a (loop [i 0 s 0.0]
                               (if (= i nlon)
                                 s
                                 (let [e (aget east (+ base i))
                                       n (aget north (+ base i))]
                                   (recur (inc i) (+ s (* e e) (* n n))))))))))))))

(defn enstrophy
  "`1/2 integral omega^2 dA`. The other inviscid invariant of two
  dimensional flow, and the one whose conservation is why energy moves to
  *large* scales here rather than small ones."
  ^double [{:keys [grid ^doubles vorticity]}]
  (let [{:keys [nlat nlon ^doubles area]} grid
        nlat (long nlat) nlon (long nlon)]
    (loop [j 0 acc 0.0]
      (if (= j nlat)
        (* 0.5 acc)
        (let [a (aget area j) base (* j nlon)]
          (recur (inc j)
                 (+ acc (* a (loop [i 0 s 0.0]
                               (if (= i nlon)
                                 s
                                 (let [w (aget vorticity (+ base i))]
                                   (recur (inc i) (+ s (* w w))))))))))))))

(defn divergence
  "The divergence of the velocity field. Identically zero, by construction.

  Worth having anyway, and worth explaining, because *which* discrete
  divergence is meant is the whole content of the claim.

  Write the meridional transport as `G = v cos(phi)`. Then

      u_east = -(1/R) D_phi(psi)        G = (1/R) D_lam(psi)
      div u  = (1/(R cos phi)) [ D_lam(u_east) + D_phi(G) ]
             = (1/R^2 cos phi) [ -D_lam D_phi psi + D_phi D_lam psi ]

  and centered differences in two independent directions commute exactly,
  so the two terms cancel to the last bit -- not to truncation error, and
  not to a solver tolerance. That is what the streamfunction buys, and it
  is why there is no pressure anywhere in this namespace.

  It holds at the poles too, but only because `G` and not `v` is what gets
  differenced across them. `G` is a longitude derivative of a scalar, so
  the scalar pole rule applies to it unchanged; `v` is a vector component
  whose frame turns over at the pole, and differencing that instead leaves
  a first-order error in the two polar rows and nowhere else."
  ^doubles [{:keys [nlat nlon ^double radius ^double dphi ^double dlam
                    ^doubles cos-phi]}
            ^doubles east ^doubles north]
  (let [nlat (long nlat) nlon (long nlon)
        cells (* nlat nlon)
        out (a/f64 cells)
        ;; The transport, formed once so that the pole rule below has a
        ;; scalar to work on.
        transport (a/f64 cells)]
    (dotimes [j nlat]
      (let [cj (aget cos-phi j) base (* j nlon)]
        (dotimes [i nlon]
          (aset transport (+ base i) (* cj (aget north (+ base i)))))))
    (dotimes [j nlat]
      (let [cj (aget cos-phi j) base (* j nlon)]
        (dotimes [i nlon]
          (let [de (/ (- (aget east (+ base (mod (inc i) nlon)))
                         (aget east (+ base (mod (+ i nlon -1) nlon))))
                      (* 2.0 dlam))
                dg (/ (- (pole-wrapped transport nlat nlon i (inc j))
                         (pole-wrapped transport nlat nlon i (dec j)))
                      (* 2.0 dphi))]
            (aset out (+ base i) (/ (+ de dg) (* radius cj)))))))
    out))
