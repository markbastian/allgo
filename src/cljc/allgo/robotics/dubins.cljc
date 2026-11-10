(ns allgo.robotics.dubins
  "Shortest paths for cars: LaValle, *Planning Algorithms*, Cambridge
  University Press, 2006, section 15.3.

    Dubins       a car that drives only forward and turns no tighter than
                 radius rho: the shortest path between two poses is one of
                 six words -- LSL, RSR, LSR, RSL, RLR, LRL -- each a turn,
                 a straight run or another turn, and a turn (Dubins,
                 American J. Mathematics 79, 1957); the closed forms here
                 are Shkel and Lumelsky's (Computational Geometry 20, 2001)
    Reeds-Shepp  a car that may also reverse: up to five segments, each
                 forward or backward, in 48 words -- Reeds and Shepp's
                 closed forms for a few, and the rest from them by
                 driving backward, swapping left for right, and reversing
                 the path (Pacific J. Mathematics 145, 1990)

  A pose is `[x y heading]`. A path is `{:word :lengths :length}`, the
  lengths of its segments in the turning radius's units, each signed for
  Reeds-Shepp (negative backward), and `drive` follows it."
  (:require [allgo.math :as am]
            [clojure.math :as math]))

(defn- mod2pi [x] (am/wrap-2pi x))

(defn- words
  "The six Dubins words' segment lengths for the normalized problem: the
  start at the origin heading `a`, the goal at (`d`, 0) heading `b`,
  radius 1. A word whose lengths do not exist is left out."
  [a b d]
  (let [sa (math/sin a) sb (math/sin b) ca (math/cos a) cb (math/cos b) cab (math/cos (- a b))
        when-real (fn [p2 f] (when (>= p2 0.0) (f (math/sqrt p2))))]
    (keep identity
          [(when-real (+ 2.0 (* d d) (* -2.0 cab) (* 2.0 d (- sa sb)))
                      (fn [p] (let [t (math/atan2 (- cb ca) (+ d sa (- sb)))]
                                {:word "LSL" :lengths [(mod2pi (- t a)) p (mod2pi (- b t))]})))
           (when-real (+ 2.0 (* d d) (* -2.0 cab) (* 2.0 d (- sb sa)))
                      (fn [p] (let [t (math/atan2 (- ca cb) (+ d (- sa) sb))]
                                {:word "RSR" :lengths [(mod2pi (- a t)) p (mod2pi (- t b))]})))
           (when-real (+ -2.0 (* d d) (* 2.0 cab) (* 2.0 d (+ sa sb)))
                      (fn [p] (let [t (- (math/atan2 (- (- ca) cb) (+ d sa sb)) (math/atan2 -2.0 p))]
                                {:word "LSR" :lengths [(mod2pi (- t a)) p (mod2pi (- t b))]})))
           (when-real (+ -2.0 (* d d) (* 2.0 cab) (* -2.0 d (+ sa sb)))
                      (fn [p] (let [t (- (math/atan2 (+ ca cb) (- d sa sb)) (math/atan2 2.0 p))]
                                {:word "RSL" :lengths [(mod2pi (- a t)) p (mod2pi (- b t))]})))
           (let [c (/ (+ 6.0 (- (* d d)) (* 2.0 cab) (* 2.0 d (- sa sb))) 8.0)]
             (when (<= (abs c) 1.0)
               (let [p (mod2pi (- (* 2 math/PI) (math/acos c)))
                     t (mod2pi (+ (- a (math/atan2 (- ca cb) (+ d (- sa) sb))) (/ p 2.0)))]
                 {:word "RLR" :lengths [t p (mod2pi (+ (- a b t) p))]})))
           (let [c (/ (+ 6.0 (- (* d d)) (* 2.0 cab) (* 2.0 d (- sb sa))) 8.0)]
             (when (<= (abs c) 1.0)
               (let [p (mod2pi (- (* 2 math/PI) (math/acos c)))
                     t (mod2pi (+ (- (- a) (math/atan2 (- ca cb) (+ d sa (- sb)))) (/ p 2.0)))]
                 {:word "LRL" :lengths [t p (mod2pi (+ (- b a t) p))]})))])))

(defn segment
  "Where `pose` ends after a segment of `kind` -- \\L, \\R or \\S -- of
  length `s` (negative backward), turning radius `rho`."
  [[x y th] kind s rho]
  (case kind
    \S [(+ x (* s (math/cos th))) (+ y (* s (math/sin th))) th]
    \L (let [th' (+ th (/ s rho))]
         [(+ x (* rho (- (math/sin th') (math/sin th)))) (- y (* rho (- (math/cos th') (math/cos th)))) th'])
    \R (let [th' (- th (/ s rho))]
         [(- x (* rho (- (math/sin th') (math/sin th)))) (+ y (* rho (- (math/cos th') (math/cos th)))) th'])))

(defn drive
  "The poses along `path` from `start`, `n` a segment, radius `rho`; the
  last is where the path ends."
  [start {:keys [word lengths]} rho n]
  (reduce (fn [poses [kind len]]
            (let [p (peek poses)]
              (into poses (map #(segment p kind (* rho len (/ % n)) rho) (range 1 (inc n))))))
          [start] (map vector word lengths)))

(defn dubins-paths
  "Every Dubins word that joins `start` to `goal` for radius `rho`, with
  its total length in metres, shortest first."
  [[x0 y0 th0] [x1 y1 th1] rho]
  (let [dx (- x1 x0) dy (- y1 y0)
        d (/ (math/hypot dx dy) rho)
        theta (mod2pi (math/atan2 dy dx))
        a (mod2pi (- th0 theta)) b (mod2pi (- th1 theta))]
    (sort-by :length (map #(assoc % :length (* rho (reduce + (:lengths %)))) (words a b d)))))

(defn dubins
  "The shortest path for a forward-only car of turning radius `rho` from
  pose `start` to pose `goal`."
  [start goal rho]
  (first (dubins-paths start goal rho)))

;; ---------------------------------------------------------------------------
;; Reeds-Shepp

(defn- polar [x y] [(math/hypot x y) (math/atan2 y x)])

(defn- tau-omega [u v xi eta phi]
  (let [delta (am/wrap-angle (- u v))
        a (- (math/sin u) (math/sin delta))
        b (- (math/cos u) (math/cos delta) 1.0)
        t1 (math/atan2 (- (* eta a) (* xi b)) (+ (* xi a) (* eta b)))
        t2 (+ (* 2.0 (- (math/cos delta) (math/cos v) (math/cos u))) 3.0)
        tau (if (neg? t2) (am/wrap-angle (+ t1 math/PI)) (am/wrap-angle t1))]
    [tau (am/wrap-angle (- (+ (- tau u) v) phi))]))

(defn- base-words
  "Reeds and Shepp's base formulas (\"Optimal Paths for a Car That Goes
  Both Forwards and Backwards\", Pacific J. Mathematics 145, 1990, section
  8) for the goal `[x y phi]` with the start at the origin heading 0,
  radius 1: `[word lengths]` for each that exists, lengths signed."
  [x y phi]
  (let [m am/wrap-angle
        hp (/ math/PI 2)
        sp (math/sin phi) cp (math/cos phi)]
    (keep identity
          [;; 8.1 C S C, the turns the same way
           (let [[u t] (polar (- x sp) (+ y -1.0 cp)) v (m (- phi t))]
             (when (and (>= t 0) (>= v 0)) ["LSL" [t u v]]))
           ;; 8.2 C S C, the turns opposite
           (let [[u1 t1] (polar (+ x sp) (- y 1.0 cp))]
             (when (>= (* u1 u1) 4.0)
               (let [u (math/sqrt (- (* u1 u1) 4.0)) t (m (+ t1 (math/atan2 2.0 u))) v (m (- t phi))]
                 (when (and (>= t 0) (>= v 0)) ["LSR" [t u v]]))))
           ;; 8.3 C|C C
           (let [[u1 theta] (polar (- x sp) (+ y -1.0 cp))]
             (when (<= u1 4.0)
               (let [u (* -2.0 (math/asin (/ u1 4.0))) t (m (+ theta (/ u 2.0) math/PI)) v (m (+ (- phi t) u))]
                 (when (and (>= t 0) (<= u 0)) ["LRL" [t u v]]))))
           ;; 8.7 C Cu|Cu C
           (let [xi (+ x sp) eta (- y 1.0 cp) rho (/ (+ 2.0 (math/hypot xi eta)) 4.0)]
             (when (<= rho 1.0)
               (let [u (math/acos rho) [t v] (tau-omega u (- u) xi eta phi)]
                 (when (and (>= t 0) (<= v 0)) ["LRLR" [t u (- u) v]]))))
           ;; 8.8 C|Cu Cu|C
           (let [xi (+ x sp) eta (- y 1.0 cp) rho (/ (- 20.0 (* xi xi) (* eta eta)) 16.0)]
             (when (<= 0.0 rho 1.0)
               (let [u (- (math/acos rho))]
                 (when (>= u (- hp))
                   (let [[t v] (tau-omega u u xi eta phi)]
                     (when (and (>= t 0) (>= v 0)) ["LRLR" [t u u v]]))))))
           ;; 8.9 C|C(pi/2) S C, the last turn the first's way
           (let [xi (- x sp) eta (+ y -1.0 cp) [rho theta] (polar (- eta) xi)]
             (when (>= rho 2.0)
               (let [r (math/sqrt (- (* rho rho) 4.0)) u (- 2.0 r) t (m (+ theta (math/atan2 2.0 r))) v (m (- phi hp t))]
                 (when (and (>= t 0) (<= u 0) (<= v 0)) ["LRSL" [t (- hp) u v]]))))
           ;; 8.10 C|C(pi/2) S C, the last turn the other way
           (let [xi (+ x sp) eta (- y 1.0 cp) [rho theta] (polar (- eta) xi)]
             (when (>= rho 2.0)
               (let [t theta u (- 2.0 rho) v (m (- (+ t hp) phi))]
                 (when (and (>= t 0) (<= u 0) (<= v 0)) ["LRSR" [t (- hp) u v]]))))
           ;; 8.11 C|C(pi/2) S C(pi/2)|C
           (let [xi (+ x sp) eta (- y 1.0 cp) [rho _] (polar xi eta)]
             (when (>= rho 2.0)
               (let [u (- 4.0 (math/sqrt (- (* rho rho) 4.0)))]
                 (when (<= u 0)
                   (let [t (m (math/atan2 (- (* (- 4.0 u) xi) (* 2.0 eta)) (+ (* -2.0 xi) (* (- u 4.0) eta)))) v (m (- t phi))]
                     (when (and (>= t 0) (>= v 0)) ["LRSLR" [t (- hp) u (- hp) v]]))))))])))

(defn- variants
  "Every word the base formulas give for `goal`, through the symmetries:
  time-flip (drive the word backward), reflection (swap left and right),
  both, and, for the words that are not their own reverse, the reversed
  path to the backward goal."
  [[x y phi]]
  (let [flip-turns (fn [w] (apply str (map {\L \R \R \L \S \S} w)))
        ;; the reversed problem: the goal seen from itself, driven in reverse order
        back [(+ (* x (math/cos phi)) (* y (math/sin phi))) (- (* x (math/sin phi)) (* y (math/cos phi))) phi]
        from (fn [[gx gy gp] reverse?]
               (for [[tf rf] [[false false] [true false] [false true] [true true]]
                     :let [g [(if tf (- gx) gx) (if rf (- gy) gy) (if (not= tf rf) (- gp) gp)]]
                     [w ls] (base-words (g 0) (g 1) (g 2))]
                 (let [w (if rf (flip-turns w) w)
                       ls (if tf (mapv - ls) ls)]
                   (if reverse? [(apply str (reverse w)) (vec (rseq ls))] [w ls]))))]
    (concat (from [x y phi] false) (from back true))))

(defn reeds-shepp-paths
  "Every Reeds-Shepp path from `start` to `goal` for turning radius `rho`
  that the formulas and their symmetries give: `{:word :lengths :length}`,
  lengths signed (negative backward) in units of `rho`, shortest first."
  [[x0 y0 th0] [x1 y1 th1] rho]
  (let [c (math/cos th0) s (math/sin th0)
        dx (/ (- x1 x0) rho) dy (/ (- y1 y0) rho)
        goal [(+ (* c dx) (* s dy)) (- (* c dy) (* s dx)) (am/wrap-angle (- th1 th0))]]
    (sort-by :length (for [[w ls] (variants goal)]
                       {:word w :lengths ls :length (* rho (reduce + (map abs ls)))}))))

(defn reeds-shepp
  "The shortest path for a car that may reverse, turning radius `rho`,
  from `start` to `goal`."
  [start goal rho]
  (first (reeds-shepp-paths start goal rho)))
