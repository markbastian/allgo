(ns allgo.physics.vortex
  "Vortex particles: little spinners that stir a velocity field.

  A grid fluid at a playable resolution loses its small eddies. The
  projection is a smoothing operator and advection interpolates, so every
  step rounds the curl off a little more, and what is left rises in flat
  sheets. Real flames and smoke are not flat -- the detail is the point.

  Rather than raise the resolution until the eddies survive, put them back
  by hand. Each vortex is a point carrying an angular velocity; it drifts
  with the flow like anything else, and within its radius it pulls the
  velocity field toward solid-body rotation about itself. They cost
  nothing next to the solve, and they are the difference between smoke
  that curls and smoke that does not.

  Nothing here is specific to fire. Any `allgo.physics.fluid` can be
  stirred this way, and the pool is a flat array with the usual
  compaction, so vortices can be spawned and expired every frame without
  allocating.

  A vortex is spawned with a lifetime and fades as it ages, because a
  spinner that never expires becomes a permanent feature of the flow and
  reads as a mistake rather than turbulence."
  (:require [allgo.physics.fluid :as fluid]
            [clojure.math :as math]))

(defn- f32 [n] #?(:clj (float-array n) :cljs (js/Float32Array. n)))

(defn pool
  "Room for `capacity` vortices."
  [capacity]
  {:capacity (long capacity)
   :count (volatile! 0)
   :x     (f32 capacity)
   :y     (f32 capacity)
   :omega (f32 capacity)
   :life  (f32 capacity)})

(defn vortex-count [{:keys [count]}] @count)

(defn vortex
  "`{:x :y :omega :life}` for vortex `i`."
  [{:keys [^floats x ^floats y ^floats omega ^floats life]} i]
  {:x (aget x i) :y (aget y i) :omega (aget omega i) :life (aget life i)})

(defn add!
  "Spawns a vortex, or does nothing if the pool is full.

  The reference never checks: its capacity field is left undefined, so the
  guard compares against nothing and always passes. Writes past the end of
  a typed array are silently dropped, so the extra vortices exist as
  counted slots holding no data, and read back as NaN until the next
  sweep quietly culls them."
  [{:keys [capacity count ^floats x ^floats y ^floats omega ^floats life]} px py w lifetime]
  (let [i @count]
    (when (< i (long capacity))
      (aset x i (float px))
      (aset y i (float py))
      (aset omega i (float w))
      (aset life i (float lifetime))
      (vreset! count (inc i))
      i)))

(defn expire!
  "Ages every vortex by `dt` and compacts the survivors to the front."
  [{:keys [count ^floats x ^floats y ^floats omega ^floats life]} dt]
  (let [n (long @count)]
    (loop [i 0 kept 0]
      (if (= i n)
        (vreset! count kept)
        (let [remaining (- (aget life i) (double dt))]
          (if (pos? remaining)
            (do (aset life kept (float remaining))
                (aset x kept (aget x i))
                (aset y kept (aget y i))
                (aset omega kept (aget omega i))
                (recur (inc i) (inc kept)))
            (recur (inc i) kept)))))))

(def default-world
  {:radius  0.05
   ;; How fast a vortex is slowed by the flow it sits in, per second.
   :damping 10.0
   ;; The fraction of its radius over which its influence falls to zero.
   ;; A vortex that stops abruptly at its rim stamps a visible disc on the
   ;; flow.
   :falloff 0.2})

(defn- blend
  "1 inside the core, falling to 0 at the rim."
  ^double [^double r ^double radius ^double falloff]
  (let [core (* radius (- 1.0 falloff))]
    (if (<= r core)
      1.0
      (max 0.0 (/ (- radius r) (max 1e-12 (* radius falloff)))))))

(defn stir!
  "Drifts every vortex with the flow, then pulls the flow around it.

  Within its radius each staggered velocity sample is moved toward what
  solid-body rotation about the vortex would give it: a velocity
  perpendicular to the offset, proportional to the distance out, plus the
  speed the vortex itself is travelling at so the whole eddy is carried
  along rather than left behind.

  Mutates the fluid's velocity arrays and the pool."
  ([p f dt] (stir! p f dt {}))
  ([{:keys [count ^floats x ^floats y ^floats omega] :as p}
    {:keys [nx ny h ^floats u ^floats v] :as f} dt world]
   (let [{:keys [radius damping falloff]} (merge default-world world)
         nx (long nx) ny (long ny)
         h  (double h) dt (double dt)
         radius (double radius) falloff (double falloff)
         max-x (* (dec nx) h)
         max-y (* (dec ny) h)
         drag (max 0.0 (- 1.0 (* (double damping) dt)))]
     (dotimes [k (long @count)]
       (let [px (aget x k) py (aget y k)
             ;; The flow the vortex is sitting in, slowed a little so a
             ;; vortex does not accelerate itself through its own field.
             vu (* drag (fluid/sample f :u px py))
             vv (* drag (fluid/sample f :v px py))
             px (min (max (+ px (* vu dt)) h) max-x)
             py (min (max (+ py (* vv dt)) h) max-y)
             w  (aget omega k)]
         (aset x k (float px))
         (aset y k (float py))
         (let [i0 (max 0 (long (math/floor (/ (- px radius) h))))
               i1 (min (dec nx) (inc (long (math/floor (/ (+ px radius) h)))))
               j0 (max 0 (long (math/floor (/ (- py radius) h))))
               j1 (min (dec ny) (inc (long (math/floor (/ (+ py radius) h)))))]
           (loop [i i0]
             (when (<= i i1)
               (loop [j j0]
                 (when (<= j j1)
                   (let [id (+ (* i ny) j)
                         ;; u sits on the left face, v on the bottom one,
                         ;; so each is a different distance from the
                         ;; vortex and needs its own offset.
                         urx (- (* i h) px)
                         ury (- (* (+ j 0.5) h) py)
                         ur  (math/sqrt (+ (* urx urx) (* ury ury)))
                         vrx (- (* (+ i 0.5) h) px)
                         vry (- (* j h) py)
                         vr  (math/sqrt (+ (* vrx vrx) (* vry vry)))]
                     (when (< ur radius)
                       (let [s (blend ur radius falloff)
                             target (+ (* ury w) vu)]
                         ;; The reference assigns here and accumulates in
                         ;; the branch below. Assigning throws away the
                         ;; velocity already in the cell and leaves
                         ;; whatever is left scaled by the falloff, which
                         ;; is not a velocity at all -- at the rim, where
                         ;; the falloff goes to zero, it erases the flow.
                         (aset u id (float (+ (aget u id)
                                              (* s (- target (aget u id))))))))
                     (when (< vr radius)
                       (let [s (blend vr radius falloff)
                             target (+ (* (- vrx) w) vv)]
                         (aset v id (float (+ (aget v id)
                                              (* s (- target (aget v id)))))))))
                   (recur (inc j))))
               (recur (inc i)))))))
     p)))

(defn step!
  "Age the vortices, then stir the fluid with the survivors."
  ([p f dt] (step! p f dt {}))
  ([p f dt world]
   (expire! p dt)
   (stir! p f dt world)
   p))
