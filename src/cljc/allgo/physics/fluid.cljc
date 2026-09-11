(ns allgo.physics.fluid
  "Eulerian fluid on a staggered grid, after Ten Minute Physics 17.

  Where `allgo.physics.xpbd` tracks particles and asks where each one
  goes, this tracks a fixed grid and asks what flows through each cell.
  Four stages a step:

    integrate    apply gravity to the velocity field
    project      make the field divergence-free, which is the whole
                 simulation: an incompressible fluid is one that neither
                 gains nor loses material anywhere
    extrapolate  carry velocities into the border cells
    advect       move the velocity field, and anything painted on it,
                 along itself

  The grid is *staggered*: pressure and dye live at cell centres, while
  the horizontal velocity `u` lives on vertical faces and the vertical
  velocity `v` on horizontal faces. That is not fussiness. With everything
  at the centre, the natural divergence of a cell does not involve that
  cell's own value, so a checkerboard of alternating pressures is
  invisible to the solver and grows without bound. On a staggered grid the
  divergence of a cell is exactly the four face velocities around it, and
  the mode cannot form.

  `s` marks which cells are open: 1 for fluid, 0 for solid. Every stage
  consults it, so obstacles are painted rather than modelled, and a moving
  obstacle is just a repainting.

  Fields are flat arrays indexed `i*ny + j`, with a one-cell border of
  solid on every side."
  (:require [clojure.math :as math]))

(defn- f32
  ([n] #?(:clj (float-array n) :cljs (js/Float32Array. n)))
  ([n fill] (let [^floats a (f32 n)] (dotimes [i n] (aset a i (float fill))) a)))

(defn fluid
  "A grid of `nx` by `ny` interior cells, each `h` across.

  A border of one cell is added on every side and marked solid, so the
  interior never has to check whether its neighbours exist."
  ([nx ny] (fluid nx ny 1.0 1000.0))
  ([nx ny h] (fluid nx ny h 1000.0))
  ([nx ny h density]
   (let [nx (+ nx 2)
         ny (+ ny 2)
         n  (* nx ny)]
     {:nx nx :ny ny :n n :h (double h) :density (double density)
      :u (f32 n) :v (f32 n)
      :u' (f32 n) :v' (f32 n)
      :p (f32 n)
      ;; Everything open by default; the border is closed by `close-border!`.
      :s (f32 n 1.0)
      :smoke (f32 n 1.0) :smoke' (f32 n 1.0)})))

(defn idx ^long [{:keys [ny]} i j] (+ (* (long i) (long ny)) (long j)))

(defn solid? [{:keys [^floats s] :as f} i j] (zero? (aget s (idx f i j))))

(defn set-solid!
  "Marks a cell solid or open. Solid cells take no velocity and pass none."
  [{:keys [^floats s] :as f} i j solid?]
  (aset s (idx f i j) (float (if solid? 0.0 1.0)))
  f)

(defn close-border!
  "Marks the outer ring solid, leaving the left column open if `inflow?`
  so a wind tunnel has somewhere to blow from."
  ([f] (close-border! f false))
  ([{:keys [nx ny] :as f} inflow?]
   (dotimes [i nx]
     (set-solid! f i 0 true)
     (set-solid! f i (dec ny) true))
   (dotimes [j ny]
     (set-solid! f 0 j (not inflow?))
     (set-solid! f (dec nx) j true))
   f))

(defn smoke-at [{:keys [^floats smoke] :as f} i j] (aget smoke (idx f i j)))
(defn set-smoke! [{:keys [^floats smoke] :as f} i j m]
  (aset smoke (idx f i j) (float m)) f)

(defn velocity-at [{:keys [^floats u ^floats v] :as f} i j]
  [(aget u (idx f i j)) (aget v (idx f i j))])

(defn set-velocity! [{:keys [^floats u ^floats v] :as f} i j ux vy]
  (aset u (idx f i j) (float ux))
  (aset v (idx f i j) (float vy))
  f)

;; ---------------------------------------------------------------------------

(defn integrate!
  "Gravity, applied only to the vertical faces between two open cells."
  [{:keys [nx ny ^floats v ^floats s]} dt gravity]
  (let [ny (long ny)]
    (dotimes [i nx]
      (when (pos? i)
        (loop [j 1]
          (when (< j (dec ny))
            (let [k (+ (* i ny) j)]
              (when (and (pos? (aget s k)) (pos? (aget s (dec k))))
                (aset v k (float (+ (aget v k) (* gravity dt))))))
            (recur (inc j))))))))

(defn project!
  "Gauss-Seidel pressure projection: the heart of it.

  For each cell, measure how much more is flowing out than in -- the
  divergence -- and push the surrounding face velocities to cancel it. Do
  that everywhere, repeatedly, and the whole field becomes divergence
  free. Nothing solves a global system here; the correction is local and
  the iteration carries it outward.

  `over-relaxation` above 1 overshoots each local correction deliberately.
  It has no physical meaning and it converges several times faster; 1.9 is
  the usual choice."
  [{:keys [nx ny h density ^floats u ^floats v ^floats p ^floats s]}
   iterations dt over-relaxation]
  (let [nx (long nx) ny (long ny)
        cp (/ (* density h) dt)]
    (dotimes [_ iterations]
      (loop [i 1]
        (when (< i (dec nx))
          (loop [j 1]
            (when (< j (dec ny))
              (let [k (+ (* i ny) j)]
                (when (pos? (aget s k))
                  (let [sx0 (aget s (- k ny))
                        sx1 (aget s (+ k ny))
                        sy0 (aget s (dec k))
                        sy1 (aget s (inc k))
                        total (+ sx0 sx1 sy0 sy1)]
                    (when (pos? total)
                      (let [div (+ (- (aget u (+ k ny)) (aget u k))
                                   (- (aget v (inc k)) (aget v k)))
                            corr (* over-relaxation (/ (- div) total))]
                        (aset p k (float (+ (aget p k) (* cp corr))))
                        (aset u k (float (- (aget u k) (* sx0 corr))))
                        (aset u (+ k ny) (float (+ (aget u (+ k ny)) (* sx1 corr))))
                        (aset v k (float (- (aget v k) (* sy0 corr))))
                        (aset v (inc k) (float (+ (aget v (inc k)) (* sy1 corr)))))))))
              (recur (inc j))))
          (recur (inc i)))))))

(defn extrapolate!
  "Copies velocities into the border, so sampling near an edge has
  something to read."
  [{:keys [nx ny ^floats u ^floats v]}]
  (let [nx (long nx) ny (long ny)]
    (dotimes [i nx]
      (aset u (* i ny) (aget u (inc (* i ny))))
      (aset u (+ (* i ny) (dec ny)) (aget u (+ (* i ny) (- ny 2)))))
    (dotimes [j ny]
      (aset v j (aget v (+ ny j)))
      (aset v (+ (* (dec nx) ny) j) (aget v (+ (* (- nx 2) ny) j))))))

(defn sample
  "Bilinear sample of a field at a point in world units.

  `field` is `:u`, `:v` or `:smoke`, and the offsets differ for each
  because they live in different places: `u` on vertical faces, `v` on
  horizontal faces, dye at centres. Getting those half-cell offsets wrong
  is the classic way to end up with a fluid that drifts diagonally."
  [{:keys [nx ny h ^floats u ^floats v ^floats smoke]} field x y]
  (let [nx (long nx) ny (long ny)
        h  (double h)
        h1 (/ 1.0 h)
        h2 (* 0.5 h)
        x  (max (min (double x) (* nx h)) h)
        y  (max (min (double y) (* ny h)) h)
        [^floats f dx dy] (case field
                            :u     [u 0.0 h2]
                            :v     [v h2 0.0]
                            :smoke [smoke h2 h2])
        x0 (min (long (math/floor (* (- x dx) h1))) (dec nx))
        tx (* (- (- x dx) (* x0 h)) h1)
        x1 (min (inc x0) (dec nx))
        y0 (min (long (math/floor (* (- y dy) h1))) (dec ny))
        ty (* (- (- y dy) (* y0 h)) h1)
        y1 (min (inc y0) (dec ny))
        sx (- 1.0 tx)
        sy (- 1.0 ty)]
    (+ (* sx sy (aget f (+ (* x0 ny) y0)))
       (* tx sy (aget f (+ (* x1 ny) y0)))
       (* tx ty (aget f (+ (* x1 ny) y1)))
       (* sx ty (aget f (+ (* x0 ny) y1))))))

(defn advect-velocity!
  "Semi-Lagrangian advection: for each face, ask where the material
  arriving there came from a step ago, and take its velocity.

  Tracing backwards rather than pushing forwards is what makes this
  unconditionally stable -- the answer is always an interpolation between
  values that already exist, so nothing can grow."
  [{:keys [nx ny h ^floats u ^floats v ^floats u' ^floats v' ^floats s] :as f} dt]
  (let [nx (long nx) ny (long ny) h (double h) h2 (* 0.5 h)]
    (dotimes [i (* nx ny)] (aset u' i (aget u i)) (aset v' i (aget v i)))
    (loop [i 1]
      (when (< i nx)
        (loop [j 1]
          (when (< j ny)
            (let [k (+ (* i ny) j)]
              (when (and (pos? (aget s k)) (pos? (aget s (- k ny))) (< j (dec ny)))
                (let [x (* i h)
                      y (+ (* j h) h2)
                      uu (aget u k)
                      ;; The other component is not stored here, so it is
                      ;; averaged from the four faces around this one.
                      ;; This face is at (i*h, (j+1/2)*h); the four v
                      ;; faces nearest it are the ones at columns i-1 and
                      ;; i, rows j and j+1. Taking the other four -- the
                      ;; stencil the v branch below wants -- reads the
                      ;; velocity from a cell diagonally away, and the
                      ;; whole fluid drifts along that diagonal.
                      vv (* 0.25 (+ (aget v (- k ny)) (aget v k)
                                    (aget v (inc (- k ny))) (aget v (inc k))))]
                  (aset u' k (float (sample f :u (- x (* dt uu)) (- y (* dt vv)))))))
              (when (and (pos? (aget s k)) (pos? (aget s (dec k))) (< i (dec nx)))
                (let [x (+ (* i h) h2)
                      y (* j h)
                      ;; And this face is at ((i+1/2)*h, j*h), so the
                      ;; four u faces nearest it are at columns i and
                      ;; i+1, rows j-1 and j.
                      uu (* 0.25 (+ (aget u (dec k)) (aget u k)
                                    (aget u (+ (dec k) ny)) (aget u (+ k ny))))
                      vv (aget v k)]
                  (aset v' k (float (sample f :v (- x (* dt uu)) (- y (* dt vv))))))))
            (recur (inc j))))
        (recur (inc i))))
    (dotimes [i (* nx ny)] (aset u i (aget u' i)) (aset v i (aget v' i)))))

(defn advect-smoke!
  "The same backward trace, carrying the dye."
  [{:keys [nx ny h ^floats u ^floats v ^floats s ^floats smoke ^floats smoke'] :as f} dt]
  (let [nx (long nx) ny (long ny) h (double h) h2 (* 0.5 h)]
    (dotimes [i (* nx ny)] (aset smoke' i (aget smoke i)))
    (loop [i 1]
      (when (< i (dec nx))
        (loop [j 1]
          (when (< j (dec ny))
            (let [k (+ (* i ny) j)]
              (when (pos? (aget s k))
                (let [uu (* 0.5 (+ (aget u k) (aget u (+ k ny))))
                      vv (* 0.5 (+ (aget v k) (aget v (inc k))))
                      x  (- (+ (* i h) h2) (* dt uu))
                      y  (- (+ (* j h) h2) (* dt vv))]
                  (aset smoke' k (float (sample f :smoke x y))))))
            (recur (inc j))))
        (recur (inc i))))
    (dotimes [i (* nx ny)] (aset smoke i (aget smoke' i)))))

(def default-world
  {:dt              (/ 1.0 60.0)
   :gravity         0.0
   :iterations      40
   ;; Deliberately overshooting each local correction. No physical
   ;; meaning, several times the convergence rate.
   :over-relaxation 1.9})

(defn step!
  "Advance the fluid one frame.

  The stages are public and run in this order because each depends on the
  last: the projection needs the forces already applied, advection needs a
  divergence-free field to carry things along, and the smoke has to be
  carried by the same velocities the fluid was left with. Compose them
  yourself when a simulation needs a pass of its own in between --
  `allgo.physics.fire` adds buoyancy, cooling and vortices that way."
  ([f] (step! f default-world))
  ([{:keys [^floats p] :as f} world]
   (let [{:keys [dt gravity iterations over-relaxation]} (merge default-world world)]
     (integrate! f dt gravity)
     (dotimes [i (alength p)] (aset p i (float 0.0)))
     (project! f iterations dt over-relaxation)
     (extrapolate! f)
     (advect-velocity! f dt)
     (advect-smoke! f dt)
     f)))

;; ---------------------------------------------------------------------------
;; Measures, mostly for checking that it is behaving

(defn divergence
  "How much more flows out of a cell than into it. Zero everywhere is what
  the projection is trying to achieve."
  [{:keys [ny ^floats u ^floats v] :as f} i j]
  (let [k (idx f i j)
        ny (long ny)]
    (+ (- (aget u (+ k ny)) (aget u k))
       (- (aget v (inc k)) (aget v k)))))

(defn max-divergence
  "The worst divergence anywhere in the open interior."
  [{:keys [nx ny ^floats s] :as f}]
  (reduce max 0.0
          (for [i (range 1 (dec (long nx)))
                j (range 1 (dec (long ny)))
                :when (pos? (aget s (idx f i j)))]
            (abs (divergence f i j)))))

(defn total-smoke
  "Summed dye over the open interior. Advection moves it about; it should
  not create or destroy much of it."
  [{:keys [nx ny ^floats s ^floats smoke] :as f}]
  (reduce + 0.0
          (for [i (range 1 (dec (long nx)))
                j (range 1 (dec (long ny)))
                :when (pos? (aget s (idx f i j)))]
            (aget smoke (idx f i j)))))

;; ---------------------------------------------------------------------------
;; Scenes

(defn wind-tunnel!
  "Steady flow from the left, with dye injected in a band.

  The left column is driven every frame rather than set once: the
  projection would otherwise bleed the inflow away within a few steps."
  [{:keys [ny] :as f} speed band]
  (let [ny (long ny)
        mid (quot ny 2)
        half (max 1 (quot (* band ny) 2))]
    (dotimes [j ny]
      (set-velocity! f 1 j speed 0.0))
    (dotimes [j ny]
      (set-smoke! f 0 j (if (< (abs (- j mid)) half) 0.0 1.0)))
    f))

(defn disc!
  "Paints a solid disc, clearing whatever was solid before except the
  border. Moving it between frames is how an obstacle is dragged around."
  [{:keys [nx ny h] :as f} cx cy r]
  (let [nx (long nx) ny (long ny) h (double h)]
    (loop [i 1]
      (when (< i (dec nx))
        (loop [j 1]
          (when (< j (dec ny))
            (set-solid! f i j false)
            (let [x (* (+ i 0.5) h)
                  y (* (+ j 0.5) h)
                  dx (- x cx) dy (- y cy)]
              (when (< (+ (* dx dx) (* dy dy)) (* r r))
                (set-solid! f i j true)
                (set-smoke! f i j 1.0)
                (set-velocity! f i j 0.0 0.0)
                (set-velocity! f (inc i) j 0.0 0.0)
                (aset ^floats (:v f) (idx f i (inc j)) (float 0.0))))
            (recur (inc j))))
        (recur (inc i))))
    f))
