(ns allgo.vallado-maneuvers-test
  "Vallado chapter 6 -- transfers, plane changes, rendezvous, low thrust and
  Hill's equations -- each checked by flying it: the burns are applied to
  state vectors and the orbits propagated with the two-body propagator,
  plane changes are measured between the planes' own vectors, and the
  closed forms are held against the equations they solve, integrated."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.kepler :as kep]
            [allgo.astro.maneuvers :as mv]
            [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu c/GM-earth)
(defn- close? [a b tol] (<= (abs (- a b)) (* tol (max 1.0 (abs b)))))
(defn- deg [x] (math/to-radians x))
(defn- radius [[r]] (v3/length r))

(defn- circular
  "The state at argument of latitude `u` on a circle of radius `a`."
  ([a u] (circular a 0.0 0.0 u))
  ([a i raan u] (kep/elements->state mu {:a a :e 0.0 :i i :raan raan :argp 0.0 :nu u})))

(defn- burn
  "Add `dv` along the velocity."
  [[r v] dv]
  [r (v3/add v (v3/scale (v3/normalize v) dv))])

(defn- fly [[r v] dt] (kep/propagate mu r v dt))

(defn- circular-there?
  "Is the state on a circular orbit -- speed sqrt(mu/r), no radial part?"
  [[r v] tol]
  (and (close? (v3/length v) (math/sqrt (/ mu (v3/length r))) tol)
       (< (abs (v3/dot (v3/normalize r) (v3/normalize v))) tol)))

(def ^:private leo 6678.0)
(def ^:private geo 42164.0)

(deftest transfers
  (testing "Hohmann: one burn out, half an orbit, one burn in, and circular at the end"
    (let [{:keys [dva dvb tof]} (mv/hohmann leo geo)
          there (fly (burn (circular leo 0.0) dva) tof)]
      (is (close? (radius there) geo 1e-9))
      (is (circular-there? (burn there dvb) 1e-9))))
  (testing "Hohmann from and to ellipses, burning at their periapsis and apoapsis"
    (let [e1 0.01 e2 0.02
          a1 (/ leo (- 1.0 e1)) a2 (/ geo (+ 1.0 e2))
          {:keys [dva dvb tof]} (mv/hohmann leo geo e1 e2 0.0 math/PI)
          start (kep/elements->state mu {:a a1 :e e1 :i 0.0 :raan 0.0 :argp 0.0 :nu 0.0})
          there (burn (fly (burn start dva) tof) dvb)
          {:keys [a e]} (kep/state->elements mu (first there) (second there))]
      (is (close? (radius there) geo 1e-9))
      (is (close? a a2 1e-9))
      (is (close? e e2 1e-7))))
  (testing "bi-elliptic: three burns, and cheaper than Hohmann for a large ratio"
    (let [rf (* 15.0 leo) rb (* 40.0 leo)
          {:keys [dva dvb dvc tof]} (mv/bi-elliptic leo rb rf)
          t1 (* math/PI (math/sqrt (/ (math/pow (* 0.5 (+ leo rb)) 3) mu)))
          at-b (fly (burn (circular leo 0.0) dva) t1)
          there (fly (burn at-b dvb) (- tof t1))]
      (is (close? (radius at-b) rb 1e-9))
      (is (close? (radius there) rf 1e-9))
      (is (circular-there? (burn there (- dvc)) 1e-9))
      (let [{h1 :dva h2 :dvb} (mv/hohmann leo rf)]
        (is (< (+ dva dvb dvc) (+ h1 h2))))))
  (testing "one-tangent: reaches the final radius at the transfer anomaly asked for"
    (let [nu (deg 160)
          {:keys [dva dvb tof e a va vb]} (mv/one-tangent leo geo 0.0 0.0 nu)
          start (burn (circular leo 0.0) dva)
          there (fly start tof)
          el (kep/state->elements mu (first there) (second there))
          [r v] there
          circ (v3/scale (v3/normalize (v3/cross [0.0 0.0 1.0] r)) (math/sqrt (/ mu geo)))]
      (is (close? (radius there) geo 1e-9))
      (is (close? (:nu el) nu 1e-9))
      (is (close? (:a el) a 1e-9))
      (is (close? (:e el) e 1e-9))
      (is (close? (v3/length (second start)) va 1e-12))
      (is (close? (v3/length v) vb 1e-9))
      (is (close? (v3/distance v circ) dvb 1e-9) "the burn onto the circle")
      (let [{h1 :dva h2 :dvb h-tof :tof} (mv/hohmann leo geo)]
        (is (< tof h-tof))
        (is (> (+ dva dvb) (+ h1 h2)) "faster, and dearer"))))
  (is (nil? (mv/one-tangent leo geo 0.0 0.0 0.0)) "no tangent conic reaches it"))

(defn- plane-normal [i raan] (v3/normalize (apply v3/cross (circular leo i raan 0.0))))

(deftest plane-changes
  (testing "inclination only: the velocity turned about the radius"
    (let [v 7.5 di (deg 15)
          before [v 0.0 0.0]
          after [(* v (math/cos di)) (* v (math/sin di)) 0.0]]
      (is (close? (mv/inclination-only di v 0.0) (v3/distance before after) 1e-12)))
    (testing "and only the horizontal part turns"
      (is (close? (mv/inclination-only (deg 15) 7.5 (deg 6))
                  (* (math/cos (deg 6)) (mv/inclination-only (deg 15) 7.5 0.0)) 1e-12))))
  (testing "node only, on a circle: the burn where the planes cross"
    (let [i (deg 55) draan (deg 45) v (math/sqrt (/ mu leo))
          {:keys [dv u-init u-final] i' :i} (mv/node-only i 0.0 draan v 0.0)
          [r1 v1] (circular leo i 0.0 u-init)
          [r2 v2] (circular leo i draan u-final)]
      (is (close? i' i 1e-12))
      (is (< (v3/distance r1 r2) 1e-6) "the same point on both orbits")
      (is (close? dv (v3/distance v1 v2) 1e-9))))
  (testing "node only, on an ellipse: burn at the apex, which tilts the plane"
    (let [i (deg 55) draan (deg 45)
          {:keys [i u-final] i-new :i} (mv/node-only i 0.0123 draan 7.5 0.0)
          apex (first (circular leo (deg 55) 0.0 (/ math/PI 2)))]
      (is (> i-new (deg 55)))
      (is (< (abs (v3/dot apex (plane-normal i-new draan))) 1e-6)
          "the apex lies in the new plane, whose node has moved by draan")
      (is (close? (math/sin u-final) (/ (math/sin (deg 55)) (math/sin i)) 1e-12))))
  (testing "inclination and node together"
    (let [v (math/sqrt (/ mu leo))
          {:keys [dv u-init u-final]} (mv/inclination-and-node (deg 55) (deg 40) (deg 45) v 0.0)
          [r1 v1] (circular leo (deg 55) 0.0 u-init)
          [r2 v2] (circular leo (deg 40) (deg 45) u-final)]
      (is (< (v3/distance r1 r2) 1e-6))
      (is (close? dv (v3/distance v1 v2) 1e-9))))
  (testing "a Hohmann transfer that also turns the plane"
    (let [di (deg -28.5)
          {:keys [di1 di2 dva dvb tof] :as plan} (mv/combined leo geo di)
          atran (* 0.5 (+ leo geo))
          total (fn [s]
                  ;; the burns as vector differences, turning s of it at the start
                  (let [dv (fn [v1 v2 t] (v3/distance [v1 0.0 0.0] [(* v2 (math/cos t)) (* v2 (math/sin t)) 0.0]))]
                    (+ (dv (mv/speed leo leo) (mv/speed leo atran) (* s di))
                       (dv (mv/speed geo atran) (mv/speed geo geo) (* (- 1.0 s) di)))))
          s (/ di1 di)]
      (is (close? (+ di1 di2) di 1e-12))
      (is (close? (+ dva dvb) (total s) 1e-12))
      (is (close? tof (:tof (mv/hohmann leo geo)) 1e-12))
      (is (< (abs di1) (abs di2)) "most of the turn where the orbit is slow")
      (testing "and the numerical optimum is a minimum the closed form nearly finds"
        (let [o (mv/combined leo geo di {:optimal true})
              so (/ (:di1 o) di)]
          (is (<= (+ (:dva o) (:dvb o)) (+ 1e-9 (+ dva dvb))))
          (is (< (+ (:dva o) (:dvb o)) (+ 1e-9 (total (+ so 1e-3)))))
          (is (< (+ (:dva o) (:dvb o)) (+ 1e-9 (total (- so 1e-3)))))
          (is (< (- (total s) (total so)) (* 0.01 (total so))) (str plan)))))))

(defn- after
  "Where a body at radius `a` and angle `theta` on a circle is `t` later."
  [a theta t]
  (+ theta (* t (math/sqrt (/ mu (* a a a))))))

(deftest rendezvous
  (testing "on the same orbit: one lap of a phasing orbit and the target arrives too"
    (doseq [[k-int k-tgt] [[1 1] [2 2]]]
      (let [a 12756.0 lead (deg 20)
            {:keys [tof a-phase dv]} (mv/rendezvous-same-orbit a lead k-int k-tgt)
            laps (/ tof (kep/period mu a-phase))]
        (is (close? laps k-int 1e-12))
        (is (close? (math/sin (after a lead tof)) 0.0 1e-9))
        (is (close? (math/cos (after a lead tof)) 1.0 1e-9))
        (is (close? dv (* 2.0 (- (mv/speed a a-phase) (mv/speed a a))) 1e-12)
            "two equal burns, into the phasing orbit and back out"))))
  (testing "between circles: wait, transfer, and meet"
    (let [r1 12756.0 r2 42164.0 lead (deg 20)
          {:keys [phase-final wait tof dv]} (mv/rendezvous-coplanar r1 r2 lead 0)
          chaser (+ (after r1 0.0 wait) math/PI)
          target (after r2 lead (+ wait tof))]
      (is (>= wait 0.0))
      (is (close? (math/cos (- chaser target)) 1.0 1e-9))
      (is (close? (math/cos (- (after r2 lead wait) (after r1 0.0 wait))) (math/cos phase-final) 1e-9)
          "the lead when the transfer starts is the one it needs")
      (is (close? dv (let [{:keys [dva dvb]} (mv/hohmann r1 r2)] (+ dva dvb)) 1e-12)))))

(deftest rendezvous-across-planes
  (let [a-int 7143.0 a-tgt 42164.0 di (deg -28.5)
        u-int (deg 200) u-tgt (deg 5)
        {:keys [t-coast t-phase t-trans a-phase dv-phase dv-trans1 dv-trans2]}
        (mv/rendezvous-noncoplanar a-int a-tgt di u-int u-tgt 1 0)]
    (testing "coast to the node, then half a Hohmann transfer"
      (is (close? (after a-int u-int t-coast) (* 2 math/PI) 1e-9))
      (is (close? t-trans (:tof (mv/hohmann a-int a-tgt)) 1e-12)))
    (testing "the burns into and out of the phasing orbit"
      (is (close? dv-phase (- (mv/speed a-int a-phase) (mv/speed a-int a-int)) 1e-12))
      (is (close? dv-trans1 (abs (- (mv/speed a-int (* 0.5 (+ a-int a-tgt))) (mv/speed a-int a-phase))) 1e-12)))
    (testing "the last burn circularizes and turns the plane at once"
      (let [va (mv/speed a-tgt (* 0.5 (+ a-int a-tgt))) vc (mv/speed a-tgt a-tgt)]
        (is (close? dv-trans2 (v3/distance [va 0.0 0.0] [(* vc (math/cos di)) (* vc (math/sin di)) 0.0]) 1e-12))))
    (testing "and the two arrive at the far node together"
      (let [u-target (after a-tgt u-tgt (+ t-coast t-phase t-trans))]
        (is (close? t-phase (kep/period mu a-phase) 1e-12) "one revolution of the phasing orbit")
        (is (close? (math/cos (- u-target math/PI)) 1.0 1e-9))))))

(deftest low-thrust
  (let [r1 6678.0 r2 42164.0 v1 (math/sqrt (/ mu r1)) v2 (math/sqrt (/ mu r2))]
    (testing "in the plane, Edelbaum's delta-v is the difference of the circular speeds"
      (is (close? (:dv (mv/low-thrust r1 r2 0.0 -2.5e-7 4e-6)) (- v1 v2) 1e-12)))
    (testing "turning the plane as well costs more"
      (is (> (:dv (mv/low-thrust r1 r2 (deg 28.5) -2.5e-7 4e-6)) (- v1 v2))))
    (testing "the time: thrust at a constant rate of mass loss, integrated"
      (let [accel 4e-6 mdot -2.5e-7
            {:keys [dv tof]} (mv/low-thrust r1 r2 0.0 mdot accel)
            ;; the acceleration grows as the mass falls: a(t) = a0 / (1 + mdot t)
            n 200000 h (/ tof n)
            gained (* h (reduce + (for [k (range n) :let [t (* (+ k 0.5) h)]]
                                    (/ accel (+ 1.0 (* mdot t))))))]
        (is (close? gained dv 1e-9))))))

(defn- cw-rk4
  "The Clohessy-Wiltshire equations -- x'' = 3n^2 x + 2n y', y'' = -2n x',
  z'' = -n^2 z -- integrated numerically."
  [n r v t steps]
  (let [f (fn [[x _ z] [vx vy _]]
            [(+ (* 3 n n x) (* 2 n vy)) (* -2 n vx) (* (- n) n z)])
        h (/ t steps)]
    (loop [r r v v k 0]
      (if (= k steps)
        [r v]
        (let [a1 (f r v)
              r2 (v3/add-scaled r v (* 0.5 h)) v2 (v3/add-scaled v a1 (* 0.5 h)) a2 (f r2 v2)
              r3 (v3/add-scaled r v2 (* 0.5 h)) v3' (v3/add-scaled v a2 (* 0.5 h)) a3 (f r3 v3')
              r4 (v3/add-scaled r v3' h) v4 (v3/add-scaled v a3 h) a4 (f r4 v4)]
          (recur (v3/add r (v3/scale (v3/add (v3/add v (v3/scale v2 2.0)) (v3/add (v3/scale v3' 2.0) v4)) (/ h 6.0)))
                 (v3/add v (v3/scale (v3/add (v3/add a1 (v3/scale a2 2.0)) (v3/add (v3/scale a3 2.0) a4)) (/ h 6.0)))
                 (inc k)))))))

(deftest hills-equations
  (let [a 6968.0 t 300.0 n (kep/mean-motion mu a)]
    (testing "the closed form solves the equations it comes from"
      (let [r0 [0.5 -1.0 0.3] v0 [-0.1 -0.04 -0.02]
            [r v] (mv/hill a r0 v0 t)
            [r' v'] (cw-rk4 n r0 v0 t 3000)]
        (is (< (v3/distance r r') 1e-9))
        (is (< (v3/distance v v') 1e-12))))
    (testing "and, for a small separation, follows two-body motion"
      (let [[rt vt] (circular a 0.0)
            ;; the chaser 1 km behind and 100 m below, drifting
            rel-r [-0.1 -1.0 0.0] rel-v [0.0 0.0003 0.0]
            ;; local axes on the x axis: radial x, along-track y; the
            ;; frame turns at n, which the inertial velocity must include
            rc (v3/add rt rel-r)
            vc (v3/add vt (v3/add rel-v [(* (- n) (second rel-r)) (* n (first rel-r)) 0.0]))
            [r _] (mv/hill a rel-r rel-v t)
            [rt' vt'] (kep/propagate mu rt vt t)
            [rc' _] (kep/propagate mu rc vc t)
            d (v3/sub rc' rt')
            x-hat (v3/normalize rt') y-hat (v3/normalize vt')]
        (is (< (abs (- (v3/dot d x-hat) (first r))) 1e-3))
        (is (< (abs (- (v3/dot d y-hat) (second r))) 1e-3))))
    (testing "the velocity to intercept, and that it does"
      (let [r0 [-70.0 20.0 -11.0]
            v0 (mv/hill-intercept a r0 t)
            [r _] (mv/hill a r0 v0 t)]
        (is (every? #(< (abs %) 1e-9) r))))))

(defn- burn-apsides
  "The apsides, found from the state vector, after a burn of `dv` at
  angle `alpha` from the velocity at radius `r`, speed `v` and flight-path
  angle `fpa`, all in the x-y plane."
  [r v fpa dv alpha]
  (let [pos [r 0.0 0.0]
        vel [(* v (math/sin fpa)) (* v (math/cos fpa)) 0.0]
        dvv [(* dv (math/sin (+ fpa alpha))) (* dv (math/cos (+ fpa alpha))) 0.0]
        {:keys [a e]} (kep/state->elements c/GM-earth pos (v3/add vel dvv))]
    {:rp (* a (- 1.0 e)) :ra (* a (+ 1.0 e))}))

(deftest fixed-dv
  (let [mu c/GM-earth
        r 7000.0 v (math/sqrt (/ mu r))]
    (testing "fixed-dv-orbit is the orbit the burned state vector has"
      (doseq [[fpa dv alpha] [[0.0 0.5 0.0] [0.1 0.3 1.2] [-0.05 0.8 -2.0]]]
        (let [{:keys [ra rp]} (mv/fixed-dv-orbit r v fpa dv alpha)
              flown (burn-apsides r v fpa dv alpha)]
          (is (< (abs (- ra (:ra flown))) 1e-6) (str alpha))
          (is (< (abs (- rp (:rp flown))) 1e-6) (str alpha)))))
    (testing "a burn along the velocity is the tangent burn, its apoapsis where vis-viva puts it"
      (let [dv 0.5
            {:keys [ra rp]} (mv/fixed-dv-orbit r v 0.0 dv 0.0)
            a (/ 1.0 (- (/ 2.0 r) (/ (math/pow (+ v dv) 2) mu)))]
        (is (< (abs (- rp r)) 1e-8))
        (is (< (abs (- ra (- (* 2.0 a) r))) 1e-8))))
    (testing "the tangent burn reaches farthest: just short of its apoapsis two directions, either side of it; just beyond, none"
      (let [farthest (:ra (mv/fixed-dv-orbit r v 0.0 0.5 0.0))
            alphas (mv/fixed-dv-to-radius r v 0.0 0.5 (* 0.999 farthest))]
        (is (= 2 (count alphas)))
        (is (< (abs (+ (first alphas) (second alphas))) 1e-9) "symmetric about the velocity")
        (doseq [alpha alphas]
          (is (< (abs (- (:ra (burn-apsides r v 0.0 0.5 alpha)) (* 0.999 farthest))) 1e-6) (str alpha)))
        (is (empty? (mv/fixed-dv-to-radius r v 0.0 0.5 (* 1.001 farthest))))))
    (testing "a lower target, reachable two ways, each checked from the state vector"
      (let [alphas (mv/fixed-dv-to-radius r v 0.1 0.3 7500.0)]
        (is (= 2 (count alphas)))
        (doseq [alpha alphas]
          (is (< (abs (- (:ra (burn-apsides r v 0.1 0.3 alpha)) 7500.0)) 1e-6) (str alpha)))))
    (testing "out of reach"
      (is (empty? (mv/fixed-dv-to-radius r v 0.0 0.1 20000.0))))
    (testing "the largest turn: the new velocity perpendicular to the burn, and no direction turns further"
      (let [{:keys [turn alpha] v' :v} (mv/max-turn 7.5 1.0)
            new [(+ 7.5 (math/cos alpha)) (math/sin alpha)]]
        (is (< (abs (- (math/atan2 (second new) (first new)) turn)) 1e-12))
        (is (< (abs (- (math/hypot (first new) (second new)) v')) 1e-12))
        (is (< (abs (+ (* (math/cos alpha) (first new)) (* (math/sin alpha) (second new)))) 1e-12))
        (is (every? (fn [b] (<= (math/atan2 (math/sin b) (+ 7.5 (math/cos b))) (+ turn 1e-12)))
                    (map #(* % 0.01) (range 315))))))))
