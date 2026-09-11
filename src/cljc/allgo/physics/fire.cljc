(ns allgo.physics.fire
  "Fire and smoke, after Ten Minute Physics 21.

  This is `allgo.physics.fluid` with a temperature field and three things
  done to it, and it is worth seeing how little is new:

    buoyancy  hot fluid is pulled upward, toward a speed proportional to
              its temperature. This is the only reason anything rises;
              there is no gravity in a fire scene.
    cooling   temperature decays, and flame decays faster than smoke, so
              a plume turns from yellow to grey at a definite height
              rather than fading uniformly.
    vortices  `allgo.physics.vortex`, spawned at the source, because a
              grid at a playable resolution smooths its own eddies away
              and fire without eddies looks like a flat sheet.

  The temperature field *is* the fluid's scalar field -- the same array
  the smoke demo advects as dye. Nothing about advecting a scalar cares
  what the scalar means, so fire needs no new transport code at all; it
  reads the field as heat instead of as dye and adds passes that act on
  it. `temperature-at` is `allgo.physics.fluid/smoke-at` under a name that
  says what it is being used for here.

  Sources are data. An emitter is a map with a `:kind`, dispatched
  through the `emit!` multimethod, so a scene can hold a burning floor and
  a burning obstacle at once and anything else you care to add without
  touching the step."
  (:require [allgo.physics.fluid :as fluid]
            [allgo.physics.vortex :as vortex]
            [clojure.math :as math]))

(defn fire
  "A fluid sized `nx` by `ny`, cold everywhere, with room for `:vortices`
  spinners."
  [{:keys [nx ny h vortices] :or {h 0.01 vortices 200}}]
  ;; Deliberately not `close-border!`. A fire wants an open domain: a
  ;; solid ceiling over a floor-wide source makes a uniformly rising
  ;; column impossible -- there is nowhere for the fluid to go and the
  ;; projection cancels the buoyancy outright, leaving the heat pooled on
  ;; the floor. `extrapolate!` handles the open edges instead. Mark solids
  ;; yourself with `allgo.physics.fluid/disc!` or `set-solid!`.
  (let [f (fluid/fluid nx ny h)
        {:keys [^floats smoke n]} f]
    ;; The scalar field starts at 1.0 for dye, which as a temperature
    ;; would mean the whole grid is alight.
    (dotimes [i n] (aset smoke i (float 0.0)))
    (assoc f :vortices (vortex/pool vortices))))

(defn temperature-at
  "The temperature of cell `i j`, 0 cold and 1 alight."
  [f i j]
  (fluid/smoke-at f i j))

(defn set-temperature! [f i j t]
  (fluid/set-smoke! f i j t))

(defn max-temperature [{:keys [^floats smoke n]}]
  (loop [i 0 m 0.0]
    (if (= i (long n)) m (recur (inc i) (max m (aget smoke i))))))

(defn total-heat [{:keys [^floats smoke n]}]
  (loop [i 0 sum 0.0]
    (if (= i (long n)) sum (recur (inc i) (+ sum (aget smoke i))))))

;; ---------------------------------------------------------------------------
;; The passes fire adds

(def default-world
  {:dt              (/ 1.0 60.0)
   ;; No gravity: a fire rises because it is hot, not because the air
   ;; around it is heavy.
   :gravity         0.0
   :iterations      10
   :over-relaxation 1.9
   ;; The upward speed a cell at full temperature is pulled toward.
   :lift            3.0
   ;; How quickly it gets there, per second.
   :buoyancy        6.0
   ;; Flame burns out faster than the smoke it leaves behind, which is
   ;; what puts a boundary between the two.
   :fire-cooling    1.2
   :smoke-cooling   0.3
   :smoke-threshold 0.3
   :vortex-radius   0.05
   :vortex-omega    20.0
   :vortex-life     1.0
   ;; Vortices spawned per unit area per second at a source. Large
   ;; because the area of a source is small: a ring a few cells thick at
   ;; h = 0.01 is a few thousandths of a square unit, and it wants a
   ;; vortex every frame or two.
   :vortex-rate     2000.0
   :emitters        []})

(defn buoyancy-and-cooling!
  "Cools every cell and pulls the hot ones upward.

  Two cooling rates, not one: below `:smoke-threshold` the cell is smoke
  and fades slowly, above it the cell is flame and burns out quickly. A
  single rate gives a plume that dims all over instead of one that turns
  to smoke at a height."
  [{:keys [nx ny ^floats v ^floats smoke]} dt world]
  (let [{:keys [lift buoyancy fire-cooling smoke-cooling smoke-threshold]}
        (merge default-world world)
        nx (long nx) ny (long ny)
        dt (double dt)
        fire-drop  (* (double fire-cooling) dt)
        smoke-drop (* (double smoke-cooling) dt)
        accel (* (double buoyancy) dt)
        lift (double lift)
        threshold (double smoke-threshold)]
    (dotimes [i nx]
      (dotimes [j ny]
        (let [id (+ (* i ny) j)
              t  (aget smoke id)
              cooled (max 0.0 (- t (if (< t threshold) smoke-drop fire-drop)))]
          (aset smoke id (float cooled))
          (let [vel (aget v id)]
            (aset v id (float (+ vel (* accel (- (* t lift) vel)))))))))))

(defn smooth-hotspots!
  "Averages any cell that is fully alight with its diagonal neighbours.

  A source holds its cells pinned at 1.0, which is a plateau with a hard
  edge; advection carries that edge upward as a visible flat-topped slab.
  Smoothing only the pinned cells breaks the plateau without touching the
  plume above it."
  [{:keys [nx ny ^floats smoke]}]
  (let [nx (long nx) ny (long ny)]
    (dotimes [i' (- nx 2)]
      (dotimes [j' (- ny 2)]
        (let [i (inc i') j (inc j')
              id (+ (* i ny) j)]
          (when (== 1.0 (aget smoke id))
            (aset smoke id
                  (float (* 0.25 (+ (aget smoke (- id ny 1))
                                    (aget smoke (+ (- id ny) 1))
                                    (aget smoke (+ id ny 1))
                                    (aget smoke (- (+ id ny) 1))))))))))))

;; ---------------------------------------------------------------------------
;; Sources

(defmulti emit!
  "Heats the cells a source covers, and reports how many vortices to
  spawn there.

  Dispatches on `:kind`, so a scene is a vector of emitter maps and a new
  kind of source is a new method rather than a change to the step.
  Implementations should heat cells and return a seq of `[x y]` points at
  which to spawn a vortex."
  (fn [_f emitter _dt _world] (:kind emitter)))

(defmethod emit! :default [_ emitter _ _]
  (throw (ex-info "unknown emitter kind" {:emitter emitter})))

(defn- spawn-points
  "Poisson-ish: `rate` spawns per unit area per second over `cells`
  cells."
  [cells h rate dt rng]
  (let [p (* (double rate) (double h) (double h) (double dt))]
    (keep (fn [[x y]] (when (< (rng) p) [x y])) cells)))

(defmethod emit! :floor
  [{:keys [nx ny h ^floats u ^floats v ^floats smoke]}
   {:keys [rows] :or {rows 4}} dt world]
  (let [{:keys [vortex-rate rng]} (merge default-world {:rng rand} world)
        nx (long nx) ny (long ny) h (double h)
        rows (min (long rows) (dec ny))]
    (spawn-points
     (for [i (range 1 (dec nx)) j (range 1 (inc rows))]
       (do (let [id (+ (* i ny) j)]
             (aset smoke id (float 1.0))
             ;; A source that is also a wind is a jet, not a fire.
             (aset u id (float 0.0))
             (aset v id (float 0.0)))
           [(* i h) (* j h)]))
     h vortex-rate dt rng)))

(defmethod emit! :disc
  [{:keys [nx ny h ^floats smoke]}
   {:keys [x y radius] :or {radius 0.2}} dt world]
  (let [{:keys [vortex-rate rng]} (merge default-world {:rng rand} world)
        nx (long nx) ny (long ny) h (double h)
        cx (double x) cy (double y)
        outer (+ (double radius) h)
        ;; A shell, not a disc: the inside of a burning object is not
        ;; where it burns, and heating it too only makes a blob.
        inner (* 0.85 (double radius))
        i0 (max 1 (long (math/floor (/ (- cx outer) h))))
        i1 (min (- nx 2) (long (math/ceil (/ (+ cx outer) h))))
        j0 (max 1 (long (math/floor (/ (- cy outer) h))))
        j1 (min (- ny 2) (long (math/ceil (/ (+ cy outer) h))))]
    (spawn-points
     (for [i (range i0 (inc i1)) j (range j0 (inc j1))
           :let [dx (- (* (+ i 0.5) h) cx)
                 dy (- (* (+ j 0.5) h) cy)
                 d2 (+ (* dx dx) (* dy dy))]
           :when (and (<= (* inner inner) d2) (< d2 (* outer outer)))]
       (do (aset smoke (+ (* i ny) j) (float 1.0))
           [(* i h) (* j h)]))
     h vortex-rate dt rng)))

(defn ignite!
  "Runs every emitter and spawns the vortices they asked for."
  [{:keys [vortices] :as f} dt world]
  (let [{:keys [emitters vortex-omega vortex-life rng]}
        (merge default-world {:rng rand} world)]
    (doseq [e emitters
            [x y] (emit! f e dt world)]
      ;; Either direction, or the plume leans.
      (vortex/add! vortices x y (* (- (* 2.0 (rng)) 1.0) (double vortex-omega))
                   vortex-life))
    f))

;; ---------------------------------------------------------------------------

(defn step!
  "One frame: the fluid's own stages, then heat, sources and vortices.

  The fluid is solved first so that buoyancy acts on a divergence-free
  field; heating it afterwards leaves a little divergence for the next
  step's projection to clear, which is both cheaper and stabler than
  trying to heat and project at once."
  ([f] (step! f {}))
  ([{:keys [vortices] :as f} world]
   (let [{:keys [dt gravity iterations over-relaxation
                 vortex-radius] :as w}
         (merge default-world world)
         {:keys [^floats p]} f]
     (fluid/integrate! f dt gravity)
     (dotimes [i (alength p)] (aset p i (float 0.0)))
     (fluid/project! f iterations dt over-relaxation)
     (fluid/extrapolate! f)
     (fluid/advect-velocity! f dt)
     (fluid/advect-smoke! f dt)
     (buoyancy-and-cooling! f dt w)
     (ignite! f dt w)
     (smooth-hotspots! f)
     (vortex/step! vortices f dt {:radius vortex-radius})
     f)))
