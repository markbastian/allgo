(ns allgo.procedural.city.buildings
  "A building for every lot, grown by a stochastic, parametric L-system:
  Parish and Mueller's section 4.2, with Greuter's stacked floor plans
  for one kind of office.

  The axiom is the building's bounding box -- the lot, set back, and
  extruded to the height its zone allows -- and every production refines
  a box into what is inside it. So each round of rewriting is a finer
  model of the same building, and any round can be drawn: the paper's
  level of detail for free, which `rounds` exposes.

  The modules are maps. A module with a `:sym` is still to be rewritten
  and is drawn as the box it bounds; the rest are terminals:

      {:op :prism   :footprint pts :z z :h h :style s}   walls, extruded
      {:op :roof    :kind :gable/:hip :rect pts :z z :h h}
      {:op :antenna :at [x y] :z z :h h}

  Three kinds of building, by zone, as the paper has, each its own set of
  productions:

      :downtown     towers that step in as they rise -- a base, then
                    setback after setback, and an antenna on top
      :commercial   offices after Greuter et al. (from the survey of
                    Kelly and McCabe): a plan of a few primitive shapes,
                    the top section extruded from the simplest plan and
                    each section below it adding one more shape; or a
                    plain low block with plant on the roof
      :residential  a house set back from the street, with a gable or a
                    hipped roof"
  (:require [allgo.geometry.polygon :as poly]
            [allgo.random :as random]
            [clojure.math :as math]))

(def defaults
  {:seed 1
   ;; Distance from the lot line to the walls, meters.
   :setback {:downtown 2.0 :commercial 2.5 :residential 3.0}
   :floor {:downtown 3.8 :commercial 3.6 :residential 3.0}
   :max-rounds 8})

;; ---------------------------------------------------------------------------
;; Footprints

(defn- scale-about [pts c k]
  (mapv #(poly/add c (poly/scale (poly/sub % c) k)) pts))

(defn- street-frame
  "The lot's frame: `u` along its first street side, `v` inward from it,
  and that side's midpoint."
  [{:keys [points street]}]
  (let [[a b] (nth (poly/edges points) (first street))
        u (poly/unit (poly/sub b a))]
    {:u u :v [(- (double (second u))) (double (first u))] :mid (poly/lerp a b 0.5)}))

(defn- rect
  "The corners of the rectangle `[u0 u1] x [v0 v1]` in the frame at
  `origin`, counterclockwise."
  [origin u v [u0 u1 v0 v1]]
  (let [at (fn [s t] (poly/add origin (poly/add (poly/scale u s) (poly/scale v t))))]
    [(at u0 v0) (at u1 v0) (at u1 v1) (at u0 v1)]))

(defn- fit-rect
  "The rectangle `[u0 u1 v0 v1]` shrunk toward its middle until it lies
  in the convex polygon `pts`, or nil."
  [pts origin u v [u0 u1 v0 v1]]
  (let [cu (* 0.5 (+ u0 u1)) cv (* 0.5 (+ v0 v1))]
    (some (fn [k]
            (let [r [(+ cu (* k (- u0 cu))) (+ cu (* k (- u1 cu)))
                     (+ cv (* k (- v0 cv))) (+ cv (* k (- v1 cv)))]]
              (when (every? #(poly/contains-point? pts %) (rect origin u v r)) r)))
          (range 1.0 0.25 -0.05))))

;; ---------------------------------------------------------------------------
;; Productions

(defn- box [fp z h] {:footprint fp :z (double z) :h (double h)})

(def ^:private plant-room
  "Meters above a low block's roof for the plant on it."
  4.0)

(defn- tower
  "A tower's productions. `:B` is the whole building: a base of a few
  floors, then the shaft. `:T` is a shaft: a stretch of floors, then a
  narrower shaft above it, until it is too short or too thin to step in
  again, when it takes its antenna."
  [{:keys [sym footprint z h floor] :as m} rng]
  (let [c (poly/centroid footprint)]
    (case sym
      :B (let [base (* floor (math/round (random/uniform rng 2.0 5.0)))]
           (if (< h (* 3.0 base))
             [(assoc (box footprint z h) :op :prism :style :tower)]
             [(assoc (box footprint z base) :op :prism :style :base)
              (assoc m :sym :T :footprint (scale-about footprint c 0.9) :z (+ z base) :h (- h base))]))
      :T (let [tier (* floor (math/round (/ (* h (random/uniform rng 0.35 0.6)) floor)))]
           (if (or (< h (* 8.0 floor)) (< (poly/area footprint) 250.0) (< tier floor))
             [(assoc (box footprint z h) :op :prism :style :tower)
              {:op :antenna :at c :z (+ z h) :h (* h (random/uniform rng 0.1 0.25))}]
             [(assoc (box footprint z tier) :op :prism :style :tower)
              (assoc m :footprint (scale-about footprint c (random/uniform rng 0.7 0.85))
                     :z (+ z tier) :h (- h tier))])))))

(defn- office
  "Greuter's floor plans. `:B` picks the primitive shapes -- rectangles in
  the lot's street frame, and now and then an octagon -- and `:S` is the
  stack of sections: the simplest plan rises the full height, and each
  further shape stops lower, so that the plan grows as you go down."
  [{:keys [sym footprint z h lot floor] :as m} rng]
  (case sym
    :B (let [{:keys [u v mid]} (street-frame lot)
             [u0 u1 v0 v1] (poly/oriented-extent footprint mid u)
             shapes (vec (keep (fn [_]
                                 (let [w (random/uniform rng 0.35 0.9) d (random/uniform rng 0.35 0.9)
                                       su (random/uniform rng u0 (- u1 (* w (- u1 u0))))
                                       sv (random/uniform rng v0 (- v1 (* d (- v1 v0))))
                                       r [su (+ su (* w (- u1 u0))) sv (+ sv (* d (- v1 v0)))]]
                                   (when-let [r (fit-rect footprint mid u v r)]
                                     (if (< (double (rng)) 0.2)
                                       ;; An octagon in that rectangle.
                                       (let [[a b c d] r
                                             cx (* 0.5 (+ a b)) cy (* 0.5 (+ c d))
                                             rr (* 0.5 (min (- b a) (- d c)))]
                                         (mapv (fn [k]
                                                 (let [t (+ (/ math/PI 8) (* k (/ math/PI 4)))]
                                                   (poly/add mid (poly/add (poly/scale u (+ cx (* rr (math/cos t))))
                                                                           (poly/scale v (+ cy (* rr (math/sin t))))))))
                                               (range 8)))
                                       (rect mid u v r)))))
                               (range (+ 2 (random/below rng 3)))))]
         (if (empty? shapes)
           [(assoc (box footprint z h) :op :prism :style :office)]
           [(assoc m :sym :S :shapes shapes)]))
    :S (let [n (count (:shapes m))]
         (into []
               (map-indexed (fn [i fp]
                              ;; The first shape the full height; the
                              ;; others each a section lower.
                              (let [hi (if (zero? i) h (* floor (max 1.0 (math/round (/ (* h (- 1.0 (/ i (+ n 0.5)))) floor)))))]
                                (assoc (box fp z hi) :op :prism :style :office))))
               (:shapes m)))))

(defn- block
  "A low commercial block on the whole lot, and plant on its roof: the
  box it was given has room for both."
  [{:keys [footprint z h]} rng]
  (let [c (poly/centroid footprint)
        walls (- h plant-room)]
    (cond-> [(assoc (box footprint z walls) :op :prism :style :shops)]
      (< (double (rng)) 0.7)
      (conj (assoc (box (scale-about footprint c (random/uniform rng 0.2 0.4)) (+ z walls)
                        (random/uniform rng 2.0 plant-room))
                   :op :prism :style :plant)))))

(defn- house
  "A house: a rectangle set back from the street, as wide as the lot
  allows and no deeper than a house is, walls, and a roof whose ridge
  runs along its longer side."
  [{:keys [footprint z h lot floor]} rng]
  (let [{:keys [u v mid]} (street-frame lot)
        [u0 u1 v0 v1] (poly/oriented-extent footprint mid u)
        depth (min (- v1 v0) (random/uniform rng 8.0 12.0))
        width (min (- u1 u0) (random/uniform rng 8.0 15.0))
        cu (* 0.5 (+ u0 u1))
        r (fit-rect footprint mid u v [(- cu (* 0.5 width)) (+ cu (* 0.5 width)) v0 (+ v0 depth)])]
    (if-not r
      [(assoc (box footprint z (min h (* 2 floor))) :op :prism :style :house)]
      (let [pts (rect mid u v r)
            walls (* floor (if (< (double (rng)) 0.6) 1 2))
            [a b c d] r
            ;; The ridge along the longer side: rotate the corners so the
            ;; first side is the long one.
            pts (if (>= (- b a) (- d c)) pts (into (subvec pts 1) [(first pts)]))]
        [(assoc (box pts z walls) :op :prism :style :house)
         {:op :roof :kind (if (< (double (rng)) 0.65) :gable :hip)
          :rect pts :z (+ z walls) :h (random/uniform rng 2.2 4.0)}]))))

(defn- produce
  "The successors of a module that is still to be rewritten."
  [{:keys [kind] :as m} rng]
  (case kind
    :tower (tower m rng)
    :office (office m rng)
    :block (block m rng)
    :house (house m rng)))

;; ---------------------------------------------------------------------------
;; Derivation

(defn axiom
  "A lot's building before any rewriting: its bounding box. The kind of
  building and its height come from the lot's zone and the land's
  height limit there."
  [{:keys [points zone] :as lot} max-height opts rng]
  (let [{:keys [setback floor]} (merge defaults opts)
        fp (or (poly/inset points (get setback zone 2.0)) points)
        [x y] (poly/centroid points)
        cap (double (max-height x y))
        f (double (get floor zone 3.5))
        kind (case zone
               :downtown :tower
               :commercial (if (< (double (rng)) 0.55) :office :block)
               :house)
        h (case kind
            :tower (max (* 10 f) (* cap (random/uniform rng 0.45 1.0)))
            :office (* f (math/round (random/uniform rng 4.0 (max 5.0 (min 16.0 (/ cap f))))))
            :block (+ plant-room (* f (math/round (random/uniform rng 1.0 4.0))))
            ;; A house's box: two storeys and a roof.
            (+ (* 2 f) 4.0))]
    (assoc (box fp 0.0 h) :sym :B :kind kind :lot lot :floor f)))

(defn rounds
  "Every stage of a lot's building: the axiom, then the string after each
  round of rewriting, until nothing is left to rewrite. The last is the
  building; any earlier one is a coarser model of it."
  ([lot max-height] (rounds lot max-height {}))
  ([lot max-height opts]
   (let [{:keys [seed max-rounds]} (merge defaults opts)
         [x y] (poly/centroid (:points lot))
         ;; Each lot its own stream, from where it is.
         rng (random/rng seed (long (math/floor x)) (long (math/floor y)) 3)]
     (loop [string [(axiom lot max-height opts rng)] out [] k 0]
       (let [out (conj out string)]
         (if (or (not-any? :sym string) (>= k (long max-rounds)))
           out
           (recur (into [] (mapcat #(if (:sym %) (produce % rng) [%])) string)
                  out (inc k))))))))

(defn building
  "A lot's finished building: its terminal modules."
  ([lot max-height] (building lot max-height {}))
  ([lot max-height opts] (peek (rounds lot max-height opts))))

(defn at-detail
  "A lot's building after at most `k` rounds of rewriting."
  [stages k]
  (nth stages (min (long k) (dec (count stages)))))
