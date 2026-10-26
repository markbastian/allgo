(ns allgo.procedural.city.land
  "The ground a city is grown on: Parish and Mueller's input maps, as
  functions instead of images.

  Their CityEngine (\"Procedural Modeling of Cities\", SIGGRAPH 2001)
  takes everything it knows about a place as 2D image maps, in two
  kinds: geographical -- elevation, and where the land, water and parks
  are -- and sociostatistical -- population density, zoning, the street
  pattern to follow, the tallest a building may be. Here each of them is
  a function of `[x y]`, in meters, built from noise, so a city needs no
  input but a seed; anything that answers the same questions could stand
  in for this one, a scanned map included.

      :elevation  (fn [x y]) -> meters above the water
      :water?     (fn [x y]) -> true where nothing can be built
      :population (fn [x y]) -> density, 0 to 1
      :pattern    (fn [x y]) -> {:grid w :radial w :organic w :terrain w}
      :grid-angle (fn [x y]) -> the local grid's direction, radians
      :zone       (fn [x y]) -> :downtown, :commercial or :residential
      :max-height (fn [x y]) -> the tallest building allowed, meters

  The pattern weights are the paper's street-pattern control maps: where
  more than one is nonzero, the roads blend what each proposes, which
  is how a radial old town shades into a grid suburb."
  (:require [allgo.procedural.fractal :as fractal]
            [allgo.procedural.noise :as noise]
            [allgo.random :as random]
            [clojure.math :as math]))

(def defaults
  {:seed 1
   ;; Half the side of the square the city may use, meters.
   :size 1600.0
   ;; :mixed, :grid, :radial, :organic or :terrain -- New York, Paris,
   ;; the paper's basic rule, and San Francisco.
   :pattern :mixed
   :river? true
   :coast? true
   ;; Secondary centers, besides the main one.
   :centers 3
   ;; Meters of relief.
   :relief 45.0})

(defn- smoothstep ^double [^double a ^double b ^double x]
  (let [t (max 0.0 (min 1.0 (/ (- x a) (- b a))))]
    (* t t (- 3.0 (* 2.0 t)))))

(defn- field
  "fBm of `seed`, about [-1, 1], with features `scale` meters across."
  [seed scale octaves]
  (let [f (fractal/fbm (noise/gradient-basis {:seed seed}) {:octaves octaves})
        k (/ 1.0 (double scale))]
    (fn ^double [^double x ^double y] (* 1.4 (double (f (* k x) (* k y) 0.5))))))

(defn- river
  "A river across the square at a random angle, meandering: `[water?
  center-distance]`, the distance from its middle line in meters."
  [seed size]
  (let [rng (random/rng seed 71 0 0)
        a (random/uniform rng 0.0 math/PI)
        off (random/uniform rng -0.35 0.35)
        u [(math/cos a) (math/sin a)]
        v [(- (math/sin a)) (math/cos a)]
        wander (field (+ seed 7) 900.0 3)
        size (double size)]
    (fn ^double [^double x ^double y]
      (let [s (+ (* x (first u)) (* y (second u)))
            t (+ (* x (first v)) (* y (second v)))
            mid (+ (* off size) (* 0.22 size (wander s 0.0)))]
        (abs (- t mid))))))

(defn land
  "The input maps for a city of `opts` (see `defaults`), as a map of
  functions."
  ([] (land {}))
  ([opts]
   (let [{:keys [seed size pattern river? coast? centers relief]} (merge defaults opts)
         seed (long seed)
         size (double size)
         relief (double relief)
         rng (random/rng seed 13 0 0)
         ;; The coast: the sea comes in over one side of the square.
         coast-dir (random/uniform rng 0.0 (* 2.0 math/PI))
         cu [(math/cos coast-dir) (math/sin coast-dir)]
         coast-wobble (field (+ seed 3) 700.0 4)
         hills (field (+ seed 5) 1100.0 5)
         river-d (river seed size)
         river-width (fn ^double [^double x ^double y] (+ 28.0 (* 10.0 (coast-wobble (* 0.5 y) (* 0.5 x)))))
         ;; Distance inland from the coastline, meters; negative is sea.
         inland (fn ^double [^double x ^double y]
                  (if coast?
                    (- (* 0.78 size)
                       (+ (* x (first cu)) (* y (second cu)))
                       (* 0.25 size (coast-wobble x y)))
                    1e9))
         elevation (fn ^double [^double x ^double y]
                     (let [h (* relief (+ 0.55 (* 0.6 (hills x y))))
                           ;; Land rises from the shore and the river bank.
                           shore (smoothstep -40.0 500.0 (inland x y))
                           bank (if river? (smoothstep 0.0 260.0 (- (river-d x y) (river-width x y))) 1.0)]
                       (* (max 0.0 h) shore bank)))
         water? (fn [x y]
                  (let [x (double x) y (double y)]
                    (or (neg? (inland x y))
                        (and river? (< (river-d x y) (river-width x y))))))
         ;; The main center somewhere central and on land; the others
         ;; around it, closer to the middle than the edges.
         on-land (fn [p] (not (water? (first p) (second p))))
         pick (fn [r]
                (let [cands (for [_ (range 40)]
                              (let [a (random/uniform rng 0.0 (* 2.0 math/PI))
                                    d (* r (math/sqrt (double (rng))))]
                                [(* d (math/cos a)) (* d (math/sin a))]))]
                  (or (first (filter #(and (on-land %)
                                           (> (inland (first %) (second %)) 250.0)
                                           (or (not river?) (> (river-d (first %) (second %)) 120.0)))
                                     cands))
                      (first (filter on-land cands))
                      [0.0 0.0])))
         main (pick (* 0.25 size))
         others (vec (repeatedly centers #(pick (* 0.7 size))))
         hubs (into [{:at main :radius (* 0.42 size) :weight 1.0
                      :angle (random/uniform rng 0.0 (* 0.5 math/PI))}]
                    (for [p others]
                      {:at p :radius (random/uniform rng (* 0.16 size) (* 0.26 size))
                       :weight (random/uniform rng 0.45 0.8)
                       :angle (random/uniform rng 0.0 (* 0.5 math/PI))}))
         grain (field (+ seed 9) 260.0 3)
         population (fn ^double [^double x ^double y]
                      (if (water? x y)
                        0.0
                        (let [p (reduce (fn [acc {[cx cy] :at r :radius w :weight}]
                                          (let [dx (- x cx) dy (- y cy)
                                                r (double r)
                                                d2 (/ (+ (* dx dx) (* dy dy)) (* r r))]
                                            (max acc (* (double w) (math/exp (- d2))))))
                                        0.0 hubs)
                              ;; Steep ground is built on less.
                              e0 (elevation x y)
                              slope (/ (+ (abs (- (elevation (+ x 15.0) y) e0))
                                          (abs (- (elevation x (+ y 15.0)) e0)))
                                       15.0)]
                          (max 0.0 (min 1.0 (* p (+ 0.85 (* 0.3 (grain x y)))
                                               (- 1.0 (* 0.7 (smoothstep 0.15 0.5 slope)))))))))
         ;; The grid of whichever center has most say here.
         grid-angle (fn ^double [^double x ^double y]
                      (double (:angle (apply max-key
                                             (fn [{[cx cy] :at r :radius w :weight}]
                                               (let [dx (- x cx) dy (- y cy) r (double r)]
                                                 (* (double w) (math/exp (- (/ (+ (* dx dx) (* dy dy)) (* r r)))))))
                                             hubs))))
         [mx my] main
         pattern-fn (case pattern
                      :grid (constantly {:grid 1.0})
                      :radial (constantly {:radial 1.0})
                      :organic (constantly {:organic 1.0})
                      :terrain (constantly {:terrain 1.0})
                      :mixed (fn [x y]
                               (let [dx (- (double x) mx) dy (- (double y) my)
                                     ;; The old town around the main center
                                     ;; is radial, the newer districts a grid,
                                     ;; the thin outskirts grown by hand.
                                     radial (math/exp (- (/ (+ (* dx dx) (* dy dy)) (* 0.09 size size))))
                                     dense (smoothstep 0.12 0.35 (population x y))
                                     rest (- 1.0 radial)]
                                 {:radial radial
                                  :grid (* rest dense)
                                  :organic (* rest (- 1.0 dense))})))
         towers (field (+ seed 11) 400.0 2)
         zone (fn [x y]
                (let [p (population x y)
                      t (towers x y)]
                  (cond (> (+ p (* 0.12 t)) 0.82) :downtown
                        (> (+ p (* 0.15 t)) 0.5) :commercial
                        :else :residential)))
         max-height (fn ^double [^double x ^double y]
                      (let [p (population x y)]
                        (+ 9.0 (* 230.0 (math/pow (smoothstep 0.45 1.0 p) 2.0)))))]
     {:size size
      :center main
      :hubs hubs
      :elevation elevation
      :water? water?
      :population population
      :pattern pattern-fn
      :grid-angle grid-angle
      :zone zone
      :max-height max-height})))

(defn legal?
  "Whether a road may be built at `[x y]`: inside the square and not in
  the water."
  [{:keys [size water?]} [x y]]
  (let [s (double size)]
    (and (<= (- s) (double x) s) (<= (- s) (double y) s)
         (not (water? x y)))))
