(ns allgo.physics.flip
  "FLIP fluid: particles carrying the velocity, a grid enforcing
  incompressibility. After Ten Minute Physics 18.

  `allgo.physics.fluid` keeps everything on a grid, which makes the
  pressure solve easy and smears the fluid out a little more every time it
  advects. Particles have the opposite problem: they carry detail
  perfectly and have no way to stay incompressible. FLIP takes both. Each
  step the particle velocities are splatted onto a grid, the grid is made
  divergence free exactly as in 17, and the result is read back onto the
  particles -- so advection is the particles simply moving, and no
  interpolation error accumulates in it.

  How the result is read back is the whole character of the fluid:

    PIC  take the grid velocity. Stable, and every transfer averages the
         particles in a cell together, so the fluid loses its swirl and
         ends up looking like treacle.
    FLIP take the grid's *change* in velocity and add it to the velocity
         the particle already had. Nothing is averaged away, so it keeps
         every eddy -- and every bit of noise, until it looks like it is
         boiling.

  `:flip-ratio` blends them. Around 0.9 is the usual answer: lively, with
  just enough PIC to damp the noise.

  Three things beyond the transfer matter. Particles are pushed apart so
  they spread evenly rather than clumping, using the same counting sort as
  `allgo.spatial.hash`. Cells are classified fluid, air or solid each
  step, and only fluid cells are solved -- that is what gives a free
  surface, and what separates this from a fluid that fills its container.
  And drift compensation pushes back when a cell holds more particles than
  it should, which is what stops the whole body slowly compressing."
  (:require [allgo.array :as a]
            [clojure.math :as math]))

(def ^:private fluid-cell 0)
(def ^:private air-cell 1)
(def ^:private solid-cell 2)

(defn flip-fluid
  "A grid of `nx` by `ny` cells of size `h`, and room for `max-particles`.

  The particle grid used for pushing particles apart is separate and
  coarser -- sized to the particle diameter, which is what its queries
  actually need."
  [{:keys [nx ny h density particle-radius max-particles]
    :or   {density 1000.0}}]
  (let [nx (+ nx 2) ny (+ ny 2)
        n  (* nx ny)
        p-spacing (* 2.2 particle-radius)
        p-nx (inc (long (math/floor (/ (* nx h) p-spacing))))
        p-ny (inc (long (math/floor (/ (* ny h) p-spacing))))]
    {:nx nx :ny ny :n n :h (double h) :density (double density)
     :particle-radius (double particle-radius)
     :max-particles max-particles
     :count (volatile! 0)
     ;; Grid
     :u (a/f32 n) :v (a/f32 n) :u0 (a/f32 n) :v0 (a/f32 n)
     :du (a/f32 n) :dv (a/f32 n) :p (a/f32 n)
     :s (a/fill! (a/f32 n) 1.0)
     :cell-type (a/i32 n)
     :particle-density (a/f32 n)
     :rest-density (volatile! 0.0)
     ;; Particles
     :pos (a/f32 (* 2 max-particles))
     :vel (a/f32 (* 2 max-particles))
     ;; The counting-sort bins used to push particles apart
     :p-nx p-nx :p-ny p-ny :p-spacing p-spacing
     :cell-count (a/i32 (inc (* p-nx p-ny)))
     :cell-ids (a/i32 max-particles)}))

(defn particle-count [f] @(:count f))

(defn add-particle! [{:keys [^floats pos ^floats vel max-particles count]} x y]
  (let [i @count]
    (when (< i max-particles)
      (aset pos (* 2 i) (float x))
      (aset pos (inc (* 2 i)) (float y))
      (aset vel (* 2 i) (float 0.0))
      (aset vel (inc (* 2 i)) (float 0.0))
      (vreset! count (inc i))
      i)))

(defn particle [{:keys [^floats pos]} i]
  [(aget pos (* 2 i)) (aget pos (inc (* 2 i)))])

(defn set-solid! [{:keys [^floats s ny]} i j solid?]
  (aset s (+ (* (long i) (long ny)) (long j)) (float (if solid? 0.0 1.0))))

(defn close-border! [{:keys [nx ny] :as f}]
  (dotimes [i nx] (set-solid! f i 0 true) (set-solid! f i (dec ny) true))
  (dotimes [j ny] (set-solid! f 0 j true) (set-solid! f (dec nx) j true))
  f)

(defn fill-block!
  "Seeds a rectangle of particles on a staggered lattice, which packs
  more evenly than a square one."
  [{:keys [particle-radius] :as f} x0 y0 x1 y1]
  (let [d (* 2.0 particle-radius)
        dy (* d (/ (math/sqrt 3.0) 2.0))]
    (loop [y (+ y0 particle-radius) row 0]
      (when (< y y1)
        (loop [x (+ x0 particle-radius (if (odd? row) particle-radius 0.0))]
          (when (< x x1)
            (add-particle! f x y)
            (recur (+ x d))))
        (recur (+ y dy) (inc row))))
    f))

;; ---------------------------------------------------------------------------

(defn- integrate-particles! [{:keys [^floats pos ^floats vel count]} dt gravity]
  (dotimes [i @count]
    (let [b (* 2 i)]
      (aset vel (inc b) (float (+ (aget vel (inc b)) (* dt gravity))))
      (aset pos b (float (+ (aget pos b) (* dt (aget vel b)))))
      (aset pos (inc b) (float (+ (aget pos (inc b)) (* dt (aget vel (inc b)))))))))

(defn- push-apart!
  "Spreads particles that have piled up.

  The same counting sort as `allgo.spatial.hash`: bin the particles,
  prefix-sum the counts into offsets, fill backward. Without this the
  particles clump into the middle of cells and leave gaps the pressure
  solve reads as empty space."
  [{:keys [^floats pos count particle-radius p-nx p-ny p-spacing
           ^ints cell-count ^ints cell-ids]}
   iterations]
  (let [n    @count
        inv  (/ 1.0 p-spacing)
        p-nx (long p-nx) p-ny (long p-ny)
        bin  (fn ^long [i]
               (let [b  (* 2 i)
                     xi (min (max 0 (long (math/floor (* (aget pos b) inv)))) (dec p-nx))
                     yi (min (max 0 (long (math/floor (* (aget pos (inc b)) inv)))) (dec p-ny))]
                 (+ (* xi p-ny) yi)))
        min-d  (* 2.0 particle-radius)
        min-d2 (* min-d min-d)]
    (a/ifill! cell-count 0)
    (dotimes [i n] (let [c (bin i)] (aset cell-count c (inc (aget cell-count c)))))
    (loop [i 0 first 0]
      (when (<= i (* p-nx p-ny))
        (let [first (if (= i (* p-nx p-ny)) first (+ first (aget cell-count i)))]
          (aset cell-count i (int first))
          (recur (inc i) first))))
    (dotimes [i n]
      (let [c (bin i)
            slot (dec (aget cell-count c))]
        (aset cell-count c (int slot))
        (aset cell-ids slot (int i))))

    (dotimes [_ iterations]
      (dotimes [i n]
        (let [b  (* 2 i)
              px (aget pos b) py (aget pos (inc b))
              xi (min (max 0 (long (math/floor (* px inv)))) (dec p-nx))
              yi (min (max 0 (long (math/floor (* py inv)))) (dec p-ny))
              x-lo (max 0 (dec xi)) x-hi (min (dec p-nx) (inc xi))
              y-lo (max 0 (dec yi)) y-hi (min (dec p-ny) (inc yi))]
          (loop [x x-lo]
            (when (<= x x-hi)
              (loop [y y-lo]
                (when (<= y y-hi)
                  (let [c    (+ (* x p-ny) y)
                        stop (aget cell-count (inc c))]
                    (loop [k (aget cell-count c)]
                      (when (< k stop)
                        (let [j (aget cell-ids k)]
                          (when (not= j i)
                            (let [q  (* 2 j)
                                  dx (- (aget pos q) (aget pos b))
                                  dy (- (aget pos (inc q)) (aget pos (inc b)))
                                  d2 (+ (* dx dx) (* dy dy))]
                              (cond
                                ;; Exactly coincident, which the tutorial
                                ;; skips and so leaves stuck together for
                                ;; good. Clamping against a wall puts
                                ;; particles on the same point routinely --
                                ;; a corner collapses them onto one spot --
                                ;; so they are nudged apart along a fixed
                                ;; axis instead. Deterministic, so the
                                ;; simulation stays repeatable.
                                (zero? d2)
                                (let [nudge (* 0.5 min-d)
                                      ;; Alternating the axis keeps a pair
                                      ;; collapsed onto a wall from being
                                      ;; nudged straight back into it.
                                      axis (if (even? i) 0 1)]
                                  (aset pos (+ b axis) (float (- (aget pos (+ b axis)) nudge)))
                                  (aset pos (+ q axis) (float (+ (aget pos (+ q axis)) nudge))))

                                (< d2 min-d2)
                                (let [d (math/sqrt d2)
                                      sc (* 0.5 (/ (- min-d d) d))
                                      ox (* dx sc) oy (* dy sc)]
                                  (aset pos b (float (- (aget pos b) ox)))
                                  (aset pos (inc b) (float (- (aget pos (inc b)) oy)))
                                  (aset pos q (float (+ (aget pos q) ox)))
                                  (aset pos (inc q) (float (+ (aget pos (inc q)) oy))))

                                :else nil))))
                        (recur (inc k)))))
                  (recur (inc y))))
              (recur (inc x)))))))))

(defn- collide-walls!
  "Keeps particles inside the box, and outside the obstacle if there is
  one."
  [{:keys [^floats pos ^floats vel count nx ny h particle-radius]} obstacle]
  (let [r  (double particle-radius)
        lo (+ h r)
        hx (- (* (dec (long nx)) h) r)
        hy (- (* (dec (long ny)) h) r)]
    (dotimes [i @count]
      (let [b (* 2 i)
            x (aget pos b) y (aget pos (inc b))]
        (when obstacle
          (let [[ox oy orad] obstacle
                dx (- x ox) dy (- y oy)
                d2 (+ (* dx dx) (* dy dy))
                rr (+ orad r)]
            (when (< d2 (* rr rr))
              (let [d (math/sqrt (max d2 1e-12))
                    s (/ rr d)]
                (aset pos b (float (+ ox (* dx s))))
                (aset pos (inc b) (float (+ oy (* dy s))))
                (aset vel b (float 0.0))
                (aset vel (inc b) (float 0.0))))))
        (let [x (aget pos b) y (aget pos (inc b))]
          (when (< x lo) (aset pos b (float lo)) (aset vel b (float 0.0)))
          (when (> x hx) (aset pos b (float hx)) (aset vel b (float 0.0)))
          (when (< y lo) (aset pos (inc b) (float lo)) (aset vel (inc b) (float 0.0)))
          (when (> y hy) (aset pos (inc b) (float hy)) (aset vel (inc b) (float 0.0))))))))

(defn- update-density!
  "How many particles each cell holds, splatted bilinearly.

  Compared against the rest density -- measured once, when the fluid is
  first at rest -- this says which cells have been squeezed, and the
  projection pushes back on those. Without it the body of fluid slowly
  loses volume, because nothing else notices particles crowding."
  [{:keys [^floats pos count nx ny h ^floats particle-density
           ^ints cell-type rest-density]}]
  (let [nx (long nx) ny (long ny) h (double h)
        h1 (/ 1.0 h) h2 (* 0.5 h)]
    (a/fill! particle-density 0.0)
    (dotimes [i @count]
      (let [b (* 2 i)
            x (min (max (aget pos b) h) (* (dec nx) h))
            y (min (max (aget pos (inc b)) h) (* (dec ny) h))
            x0 (min (long (math/floor (* (- x h2) h1))) (- nx 2))
            tx (* (- (- x h2) (* x0 h)) h1)
            x1 (min (inc x0) (- nx 2))
            y0 (min (long (math/floor (* (- y h2) h1))) (- ny 2))
            ty (* (- (- y h2) (* y0 h)) h1)
            y1 (min (inc y0) (- ny 2))
            sx (- 1.0 tx) sy (- 1.0 ty)]
        (aset particle-density (+ (* x0 ny) y0)
              (float (+ (aget particle-density (+ (* x0 ny) y0)) (* sx sy))))
        (aset particle-density (+ (* x1 ny) y0)
              (float (+ (aget particle-density (+ (* x1 ny) y0)) (* tx sy))))
        (aset particle-density (+ (* x1 ny) y1)
              (float (+ (aget particle-density (+ (* x1 ny) y1)) (* tx ty))))
        (aset particle-density (+ (* x0 ny) y1)
              (float (+ (aget particle-density (+ (* x0 ny) y1)) (* sx ty))))))
    ;; The first time the fluid is settled, whatever a full cell holds is
    ;; what "full" means from then on.
    (when (zero? (double @rest-density))
      ;; `count` is the particle counter here, so clojure.core/count is
      ;; shadowed; summing both quantities in one pass avoids needing it.
      (let [[sum cells]
            (reduce (fn [[sum cells] i]
                      (if (= fluid-cell (aget cell-type i))
                        [(+ sum (aget particle-density i)) (inc cells)]
                        [sum cells]))
                    [0.0 0]
                    (range (* nx ny)))]
        (when (pos? cells)
          (vreset! rest-density (/ sum cells)))))))

(defn- transfer!
  "Particles to grid, or grid back to particles.

  Going to the grid it is a weighted scatter followed by a divide, which
  is a weighted average. Coming back it is `flip-ratio` between taking the
  grid velocity outright and taking only the grid's change."
  [{:keys [^floats pos ^floats vel count nx ny h ^floats u ^floats v
           ^floats u0 ^floats v0 ^floats du ^floats dv ^floats s ^ints cell-type]}
   to-grid? flip-ratio]
  (let [nx (long nx) ny (long ny) h (double h)
        h1 (/ 1.0 h) h2 (* 0.5 h)
        n  @count]
    (when to-grid?
      (dotimes [i (alength u)] (aset u0 i (aget u i)) (aset v0 i (aget v i)))
      (a/fill! du 0.0) (a/fill! dv 0.0) (a/fill! u 0.0) (a/fill! v 0.0)
      ;; Classify: solid where the mask says so, fluid where a particle
      ;; is, air everywhere else. Only fluid cells get solved, and that is
      ;; what gives the fluid a free surface.
      (dotimes [i (* nx ny)]
        (aset cell-type i (int (if (zero? (aget s i)) solid-cell air-cell))))
      (dotimes [i n]
        (let [b (* 2 i)
              xi (min (max 0 (long (math/floor (* (aget pos b) h1)))) (dec nx))
              yi (min (max 0 (long (math/floor (* (aget pos (inc b)) h1)))) (dec ny))
              c  (+ (* xi ny) yi)]
          (when (= air-cell (aget cell-type c))
            (aset cell-type c (int fluid-cell))))))

    (dotimes [component 2]
      (let [dx (if (zero? component) 0.0 h2)
            dy (if (zero? component) h2 0.0)
            ^floats f  (if (zero? component) u v)
            ^floats f0 (if (zero? component) u0 v0)
            ^floats d  (if (zero? component) du dv)
            offset (if (zero? component) ny 1)]
        (dotimes [i n]
          (let [b (* 2 i)
                x (min (max (aget pos b) h) (* (dec nx) h))
                y (min (max (aget pos (inc b)) h) (* (dec ny) h))
                x0 (min (long (math/floor (* (- x dx) h1))) (- nx 2))
                tx (* (- (- x dx) (* x0 h)) h1)
                x1 (min (inc x0) (- nx 2))
                y0 (min (long (math/floor (* (- y dy) h1))) (- ny 2))
                ty (* (- (- y dy) (* y0 h)) h1)
                y1 (min (inc y0) (- ny 2))
                sx (- 1.0 tx) sy (- 1.0 ty)
                d0 (* sx sy) d1 (* tx sy) d2 (* tx ty) d3 (* sx ty)
                n0 (+ (* x0 ny) y0) n1 (+ (* x1 ny) y0)
                n2 (+ (* x1 ny) y1) n3 (+ (* x0 ny) y1)]
            (if to-grid?
              (let [pv (aget vel (+ b component))]
                (aset f n0 (float (+ (aget f n0) (* pv d0)))) (aset d n0 (float (+ (aget d n0) d0)))
                (aset f n1 (float (+ (aget f n1) (* pv d1)))) (aset d n1 (float (+ (aget d n1) d1)))
                (aset f n2 (float (+ (aget f n2) (* pv d2)))) (aset d n2 (float (+ (aget d n2) d2)))
                (aset f n3 (float (+ (aget f n3) (* pv d3)))) (aset d n3 (float (+ (aget d n3) d3))))
              ;; A face is only usable if a cell on one side of it holds
              ;; fluid; reading air would drag the surface toward zero.
              (let [valid (fn ^double [nr]
                            (if (or (not= air-cell (aget cell-type nr))
                                    (not= air-cell (aget cell-type (- nr offset))))
                              1.0 0.0))
                    v0* (valid n0) v1* (valid n1) v2* (valid n2) v3* (valid n3)
                    dd  (+ (* v0* d0) (* v1* d1) (* v2* d2) (* v3* d3))]
                (when (pos? dd)
                  (let [pv  (aget vel (+ b component))
                        pic (/ (+ (* v0* d0 (aget f n0)) (* v1* d1 (aget f n1))
                                  (* v2* d2 (aget f n2)) (* v3* d3 (aget f n3))) dd)
                        corr (/ (+ (* v0* d0 (- (aget f n0) (aget f0 n0)))
                                   (* v1* d1 (- (aget f n1) (aget f0 n1)))
                                   (* v2* d2 (- (aget f n2) (aget f0 n2)))
                                   (* v3* d3 (- (aget f n3) (aget f0 n3)))) dd)]
                    (aset vel (+ b component)
                          (float (+ (* (- 1.0 flip-ratio) pic)
                                    (* flip-ratio (+ pv corr)))))))))))

        (when to-grid?
          (dotimes [i (alength f)]
            (when (pos? (aget d i)) (aset f i (float (/ (aget f i) (aget d i))))))
          ;; A face touching a solid keeps the velocity it had, so the wall
          ;; neither absorbs nor emits fluid.
          (dotimes [i nx]
            (dotimes [j ny]
              (let [k (+ (* i ny) j)
                    solid? (= solid-cell (aget cell-type k))]
                (when (or solid? (and (pos? i) (= solid-cell (aget cell-type (- k ny)))))
                  (aset u k (aget u0 k)))
                (when (or solid? (and (pos? j) (= solid-cell (aget cell-type (dec k)))))
                  (aset v k (aget v0 k)))))))))))

(defn- project!
  "The same local Gauss-Seidel projection as `allgo.physics.fluid`, with
  two differences: only fluid cells are solved, which is what lets the
  fluid have a surface, and a cell holding more than the rest density is
  pushed back, which is what stops the body slowly compressing."
  [{:keys [nx ny h density ^floats u ^floats v ^floats u0 ^floats v0
           ^floats p ^floats s
           ^ints cell-type ^floats particle-density rest-density]}
   iterations dt over-relaxation drift?]
  (let [nx (long nx) ny (long ny)
        cp (/ (* density h) dt)
        rest-d (double @rest-density)]
    (a/fill! p 0.0)
    ;; The FLIP correction read back in `transfer!` is the change the
    ;; *projection* made, so the reference it is measured against has to be
    ;; taken here -- after the particles have been splatted onto the grid.
    ;; `transfer!` also keeps a snapshot, but that one is taken before the
    ;; splat, because the solid-cell restore needs the velocities the walls
    ;; had last frame. Measuring the correction against that earlier
    ;; snapshot instead adds the whole particle-to-grid splat back onto
    ;; velocities that already contain it, and the simulation grows without
    ;; bound.
    (dotimes [i (alength u)]
      (aset u0 i (aget u i))
      (aset v0 i (aget v i)))
    (dotimes [_ iterations]
      (loop [i 1]
        (when (< i (dec nx))
          (loop [j 1]
            (when (< j (dec ny))
              (let [k (+ (* i ny) j)]
                (when (= fluid-cell (aget cell-type k))
                  (let [sx0 (aget s (- k ny)) sx1 (aget s (+ k ny))
                        sy0 (aget s (dec k)) sy1 (aget s (inc k))
                        total (+ sx0 sx1 sy0 sy1)]
                    (when (pos? total)
                      (let [div0 (+ (- (aget u (+ k ny)) (aget u k))
                                    (- (aget v (inc k)) (aget v k)))
                            div (if (and drift? (pos? rest-d))
                                  (let [excess (- (aget particle-density k) rest-d)]
                                    (if (pos? excess) (- div0 (* 1.0 excess)) div0))
                                  div0)
                            corr (* over-relaxation (/ (- div) total))]
                        (aset p k (float (+ (aget p k) (* cp corr))))
                        (aset u k (float (- (aget u k) (* sx0 corr))))
                        (aset u (+ k ny) (float (+ (aget u (+ k ny)) (* sx1 corr))))
                        (aset v k (float (- (aget v k) (* sy0 corr))))
                        (aset v (inc k) (float (+ (aget v (inc k)) (* sy1 corr)))))))))
              (recur (inc j))))
          (recur (inc i)))))))

(def default-world
  {:dt              (/ 1.0 60.0)
   :gravity         -9.81
   :iterations      50
   :separation      2
   :over-relaxation 1.9
   :flip-ratio      0.9
   :drift?          true
   :obstacle        nil})

(defn step!
  "One frame: move the particles, spread them out, bounce them off the
  walls, splat onto the grid, project, and read back."
  ([f] (step! f default-world))
  ([f world]
   (let [{:keys [dt gravity iterations separation over-relaxation
                 flip-ratio drift? obstacle]} (merge default-world world)]
     (integrate-particles! f dt gravity)
     (push-apart! f separation)
     (collide-walls! f obstacle)
     (transfer! f true flip-ratio)
     (update-density! f)
     (project! f iterations dt over-relaxation drift?)
     (transfer! f false flip-ratio)
     f)))

;; ---------------------------------------------------------------------------

(defn min-separation
  "The closest any two particles are, for checking that pushing apart
  works."
  [{:keys [^floats pos count]}]
  (let [n @count]
    (reduce min ##Inf
            (for [i (range n) j (range (inc i) n)
                  :let [a (* 2 i) b (* 2 j)
                        dx (- (aget pos a) (aget pos b))
                        dy (- (aget pos (inc a)) (aget pos (inc b)))]]
              (math/sqrt (+ (* dx dx) (* dy dy)))))))

(defn bounds-of-particles
  "`[[min-x min-y] [max-x max-y]]` over the particles."
  [{:keys [^floats pos count]}]
  (let [n @count]
    (reduce (fn [[[lx ly] [hx hy]] i]
              (let [b (* 2 i) x (aget pos b) y (aget pos (inc b))]
                [[(min lx x) (min ly y)] [(max hx x) (max hy y)]]))
            [[##Inf ##Inf] [##-Inf ##-Inf]]
            (range n))))
