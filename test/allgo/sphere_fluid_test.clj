(ns allgo.sphere-fluid-test
  (:require [allgo.physics.sphere-fluid :as sf]
            [clojure.test :refer [deftest is testing]]))

(defn- field
  "A field from a function of latitude and longitude."
  ^doubles [{:keys [nlat nlon ^doubles phi ^double dlam]} f]
  (let [nlat (long nlat) nlon (long nlon)
        a (double-array (* nlat nlon))]
    (dotimes [j nlat]
      (dotimes [i nlon]
        (aset a (+ (* j nlon) i) (double (f (aget phi j) (* i dlam))))))
    a))

(defn- demean ^doubles [^doubles x]
  (let [m (/ (areduce x i s 0.0 (+ s (aget x i))) (alength x))]
    (amap x i _r (- (aget x i) m))))

(defn- rel-error
  "Relative L2 difference, with the mean taken out of both -- a
  streamfunction is only defined up to a constant."
  ^double [a b]
  (let [^doubles a (demean a) ^doubles b (demean b)
        num (Math/sqrt (areduce a i s 0.0 (+ s (let [d (- (aget a i) (aget b i))] (* d d)))))
        den (Math/sqrt (areduce b i s 0.0 (+ s (* (aget b i) (aget b i)))))]
    (/ num (max 1e-300 den))))

(defn- max-abs ^double [^doubles a]
  (areduce a i s 0.0 (max s (abs (aget a i)))))

;; The spherical harmonics are eigenfunctions of the Laplace-Beltrami
;; operator with eigenvalue -l(l+1)/R^2, which gives this solver something
;; exact to be graded against rather than a finer version of itself.
(def ^:private harmonics
  [["Y(1,0)" 1 (fn [p _] (Math/sin p))]
   ["Y(1,1)" 1 (fn [p l] (* (Math/cos p) (Math/cos l)))]
   ["Y(2,0)" 2 (fn [p _] (- (* 3 (Math/sin p) (Math/sin p)) 1.0))]
   ["Y(3,2)" 3 (fn [p l] (* (Math/sin p) (Math/pow (Math/cos p) 2) (Math/cos (* 2 l))))]])

(deftest grid-test
  (let [g (sf/grid {:nlat 16 :nlon 32 :radius 2.0})]
    (testing "no sample sits on a pole, so nothing divides by zero"
      (is (every? pos? (seq ^doubles (:cos-phi g))))
      (is (< (abs (+ (* 0.5 Math/PI) (aget ^doubles (:phi g) 0))) (:dphi g)))
      (is (< (abs (- (* 0.5 Math/PI) (aget ^doubles (:phi g) 15))) (:dphi g))))

    (testing "but the poles are cell faces, where cos is exactly zero"
      (is (zero? (aget ^doubles (:cos-face g) 0)))
      (is (zero? (aget ^doubles (:cos-face g) 16))))

    (testing "the cells cover the sphere and nothing else"
      ;; Not exactly: the cell areas are a midpoint rule for the integral
      ;; of cos, so they are short by O(h^2) and converge as the grid does.
      (doseq [[n tol] [[16 2e-3] [64 2e-4] [256 1e-5]]]
        (let [g (sf/grid {:nlat n :nlon 32 :radius 2.0})
              total (* 32 (reduce + (seq ^doubles (:area g))))]
          (is (< (abs (- 1.0 (/ total (* 4 Math/PI 4.0)))) tol)))))

    (testing "and a longitude count that is not a power of two is refused"
      (is (thrown? clojure.lang.ExceptionInfo (sf/grid {:nlat 16 :nlon 30}))))))

(deftest poisson-test
  (testing "the Poisson solve reproduces the exact harmonics"
    (doseq [[label deg f] harmonics]
      (let [g (sf/grid {:nlat 64 :nlon 128})
            sim (sf/simulation {:grid g :tracer? false})
            psi (field g f)
            ev (double (- (* deg (inc deg))))
            omega (amap psi i _r (* ev (aget psi i)))
            out (double-array (:cells g))]
        (sf/solve-helmholtz! g (:workspace sim) omega 0.0 out)
        (is (< (rel-error out psi) 2e-3) label))))

  (testing "and converges at second order, which is what the differencing
            promises and the only way to tell a discretisation from a bug"
    (let [[_ deg f] (first harmonics)
          err (fn [n]
                (let [g (sf/grid {:nlat n :nlon (* 2 n)})
                      sim (sf/simulation {:grid g :tracer? false})
                      psi (field g f)
                      ev (double (- (* deg (inc deg))))
                      omega (amap psi i _r (* ev (aget psi i)))
                      out (double-array (:cells g))]
                  (sf/solve-helmholtz! g (:workspace sim) omega 0.0 out)
                  (rel-error out psi)))
          coarse (err 16)
          fine (err 32)
          finer (err 64)]
      (is (< (abs (- 4.0 (/ coarse fine))) 0.5))
      (is (< (abs (- 4.0 (/ fine finer))) 0.5))))

  (testing "a positive alpha gives the Helmholtz operator, also exactly"
    ;; (laplacian - alpha) psi = (-l(l+1)/R^2 - alpha) psi, and unlike
    ;; Poisson this one is nonsingular, so nothing has to be pinned.
    (doseq [alpha [0.5 10.0]]
      (let [g (sf/grid {:nlat 64 :nlon 128})
            sim (sf/simulation {:grid g :tracer? false})
            psi (field g (fn [p l] (* (Math/cos p) (Math/cos l))))
            ev (- -2.0 (double alpha))
            rhs (amap psi i _r (* ev (aget psi i)))
            out (double-array (:cells g))]
        (sf/solve-helmholtz! g (:workspace sim) rhs alpha out)
        (is (< (rel-error out psi) 2e-3))))))

(deftest velocity-test
  (let [g (sf/grid {:nlat 64 :nlon 128 :radius 2.0})
        omega 0.7
        psi (field g (fn [p _] (* -1.0 omega 4.0 (Math/sin p))))
        east (double-array (:cells g))
        north (double-array (:cells g))]
    (sf/velocity! g psi east north)

    (testing "solid-body rotation comes out of the streamfunction that
              describes it, and has no meridional component at all"
      (is (< (rel-error east (field g (fn [p _] (* omega 2.0 (Math/cos p))))) 1e-3))
      (is (zero? (max-abs north))))

    (testing "the velocity is divergence-free to the last bit, everywhere,
              which is the whole point of carrying a streamfunction"
      (is (< (max-abs (sf/divergence g east north)) 1e-12)))))

(deftest divergence-test
  (testing "an arbitrary vorticity field still gives an exactly
            non-divergent velocity, poles included"
    (doseq [n [16 32 64]]
      (let [sim (sf/simulation {:grid {:nlat n :nlon (* 2 n)} :tracer? false})
            g (:grid sim)]
        (sf/add-vortex! g (:vorticity sim) [1 0.3 0.2] 6.0 0.25)
        (sf/add-vortex! g (:vorticity sim) [-0.4 0.8 0.1] -4.0 0.3)
        (sf/update-velocity! sim)
        (let [d (sf/divergence g (:east sim) (:north sim))]
          (is (< (max-abs d) (* 1e-10 (max-abs (:east sim))))
              (str "at nlat " n)))))))

(deftest advection-test
  (let [g (sf/grid {:nlat 32 :nlon 64})
        cells (:cells g)]
    (testing "a constant field is carried unchanged"
      (let [sim (sf/simulation {:grid g :tracer? false})
            f (double-array cells 3.0)
            out (double-array cells)]
        (sf/add-vortex! g (:vorticity sim) [0 1 0] 5.0 0.4)
        (sf/update-velocity! sim)
        (sf/advect! g f (:east sim) (:north sim) 0.05 out)
        (is (every? #(< (abs (- 3.0 (double %))) 1e-12) (seq out)))))

    (testing "and a zero velocity leaves a field where it was"
      (let [f (field g (fn [p l] (* (Math/sin p) (Math/cos l))))
            out (double-array cells)]
        (sf/advect! g f (double-array cells) (double-array cells) 0.1 out)
        (is (< (rel-error out f) 1e-12))))

    (testing "solid-body rotation returns a field to itself after one turn"
      ;; The strongest accuracy test there is for advection: the exact
      ;; answer after a full revolution is the field you started with, so
      ;; whatever has changed is entirely the scheme's doing.
      (let [omega 1.0
            psi (field g (fn [p _] (* -1.0 omega (Math/sin p))))
            east (double-array cells)
            north (double-array cells)
            start (field g (fn [p l] (* (Math/cos p) (Math/cos (* 2 l)))))
            steps 200
            dt (/ (* 2.0 Math/PI) omega steps)]
        (sf/velocity! g psi east north)
        (let [carry (fn [scheme]
                      (let [cur (double-array cells)
                            tmp (double-array cells)
                            lo (double-array cells)
                            hi (double-array cells)
                            dp (double-array cells) dl (double-array cells)
                            fp (double-array cells) fl (double-array cells)]
                        (System/arraycopy start 0 cur 0 cells)
                        (sf/departures! g east north dt dp dl)
                        (sf/departures! g east north (- dt) fp fl)
                        (dotimes [_ steps]
                          (if (= scheme :plain)
                            (do (sf/sample-at! g cur dp dl tmp)
                                (System/arraycopy tmp 0 cur 0 cells))
                            (let [a (double-array cells)]
                              (sf/sample-at! g cur dp dl a lo hi)
                              (sf/sample-at! g a fp fl tmp)
                              (sf/correct! cur a tmp lo hi cur))))
                        cur))
              plain (carry :plain)
              corrected (carry :corrected)]
          (testing "both stay bounded by what they started with"
            (is (<= (max-abs plain) (+ 1e-9 (max-abs start))))
            (is (<= (max-abs corrected) (+ 1e-9 (max-abs start)))))
          (testing "and the corrected scheme is much the more faithful"
            (is (< (rel-error corrected start) 0.2))
            (is (< (rel-error corrected start) (* 0.5 (rel-error plain start))))))))))

(deftest conservation-test
  (let [build (fn [opts]
                (let [sim (sf/simulation (merge {:grid {:nlat 32 :nlon 64}} opts))
                      g (:grid sim)]
                  (doseq [[c s w] [[[1 0.3 0.2] 8.0 0.22] [[-0.4 0.8 0.1] -7.0 0.25]
                                   [[0.2 -0.9 0.3] 6.0 0.3]]]
                    (sf/add-vortex! g (:vorticity sim) c s w))
                  (sf/update-velocity! sim)
                  sim))]

    (testing "total vorticity stays at zero -- a curl on a closed surface
              has no edge to leak out of"
      (let [final (reduce (fn [s _] (sf/step! s 0.01)) (build {}) (range 100))]
        (is (< (abs (sf/total-vorticity final)) 1e-6))))

    (testing "the corrected scheme keeps the energy and enstrophy that the
              plain one dissipates away"
      (let [run (fn [adv]
                  (let [sim (build {:advection adv})
                        z0 (sf/enstrophy sim)
                        final (reduce (fn [s _] (sf/step! s 0.01)) sim (range 150))]
                    (/ (sf/enstrophy final) z0)))
            plain (run :semi-lagrangian)
            corrected (run :corrected)]
        (is (< plain 0.85))
        (is (> corrected 0.95))))

    (testing "viscosity dissipates enstrophy, and more of it for more of it"
      (let [run (fn [nu]
                  (let [sim (build {:viscosity nu :advection :corrected})
                        z0 (sf/enstrophy sim)
                        final (reduce (fn [s _] (sf/step! s 0.01)) sim (range 100))]
                    (/ (sf/enstrophy final) z0)))
            zs (map run [0.0 1e-4 1e-3])]
        (is (= zs (reverse (sort zs))))
        (is (< (last zs) (first zs)))))))

(defn- zonal-energy-fraction
  "How much of the kinetic energy belongs to the longitude-averaged flow.

  The diagnostic for banding, and far sharper than comparing east against
  north: a field of vortices has almost none of its energy in the zonal
  mean, and a set of jets has almost all of it."
  ^double [sim]
  (let [{:keys [nlat nlon ^doubles area]} (:grid sim)
        nlat (long nlat) nlon (long nlon)
        ^doubles e (:east sim) ^doubles n (:north sim)]
    (loop [j 0 mean-energy 0.0 total 0.0]
      (if (= j nlat)
        (/ mean-energy (max 1e-300 total))
        (let [base (* j nlon)
              a (aget area j)
              ubar (/ (loop [i 0 s 0.0]
                        (if (= i nlon) s (recur (inc i) (+ s (aget e (+ base i))))))
                      nlon)
              tot (loop [i 0 s 0.0]
                    (if (= i nlon)
                      s
                      (recur (inc i) (+ s (* (aget e (+ base i)) (aget e (+ base i)))
                                        (* (aget n (+ base i)) (aget n (+ base i)))))))]
          (recur (inc j) (+ mean-energy (* a nlon ubar ubar)) (+ total (* a tot))))))))

(deftest rotation-test
  (testing "rotation is what holds zonal jets together"
    ;; The signature of the Coriolis term, and the reason the giant
    ;; planets are striped. Seed alternating east-west jets and leave
    ;; them alone: without rotation they are barotropically unstable and
    ;; shred into a field of vortices, and with it the beta effect holds
    ;; them in place -- a parcel pushed off its latitude must give up
    ;; relative vorticity to keep omega + f, which turns it back.
    (let [zonality (fn [rot after]
                     (let [sim (sf/simulation {:grid {:nlat 48 :nlon 128}
                                               :rotation rot :tracer? false})]
                       (sf/zonal-jets! (:grid sim) (:vorticity sim) {})
                       (sf/update-velocity! sim)
                       (zonal-energy-fraction
                        (reduce (fn [s _] (sf/step! s 0.01)) sim (range after)))))]
      (is (> (zonality 0.0 0) 0.99) "jets start perfectly zonal")
      (is (< (zonality 0.0 1200) 0.6) "and without rotation they do not stay that way")
      (is (> (zonality 30.0 1200) 0.9) "with it they do"))))

(deftest coriolis-scheme-test
  (let [jets (fn [rot cor nlat nlon]
               (let [sim (sf/simulation {:grid {:nlat nlat :nlon nlon}
                                         :rotation rot :coriolis cor :tracer? false})]
                 (sf/zonal-jets! (:grid sim) (:vorticity sim) {})
                 (sf/update-velocity! sim)
                 sim))
        energy-ratio (fn [rot cor dt time]
                       (let [sim (jets rot cor 48 128)
                             e0 (sf/kinetic-energy sim)
                             f (reduce (fn [s _] (sf/step! s dt))
                                       sim (range (long (/ time dt))))]
                         (/ (sf/kinetic-energy f) e0)))]

    (testing "energy is conserved while the flow stays smooth"
      ;; An instability moves energy between scales; it never makes any.
      (doseq [rot [10.0 30.0]]
        (is (< 0.85 (energy-ratio rot :semi-implicit 0.02 6.0) 1.2)
            (str "at rotation " rot))))

    (testing "and is lost, honestly, once the flow goes turbulent"
      ;; Without rotation the jets shred, energy cascades to the grid
      ;; scale, and there it is dissipated by the interpolation and the
      ;; limiter -- there being nowhere else for it to go on a grid this
      ;; coarse. Worth pinning as a known property rather than leaving it
      ;; to be discovered as a surprise.
      (is (< 0.3 (energy-ratio 0.0 :semi-implicit 0.02 6.0) 0.9)))

    (testing "and the explicit scheme cannot manage that once the planet spins"
      ;; Rossby waves get faster with rotation, and stepping an
      ;; oscillation from the start of the interval is forward Euler,
      ;; whose amplification exceeds one at every step size. This is the
      ;; whole reason `solve-rossby!` exists, so it is worth pinning: the
      ;; explicit scheme is not merely less accurate here, it diverges.
      (let [explicit (energy-ratio 80.0 :absolute 0.02 6.0)
            implicit (energy-ratio 80.0 :semi-implicit 0.02 6.0)]
        (is (> explicit 100.0) "explicit runs away")
        (is (< 0.7 implicit 1.3) "implicit does not")))

    (testing "the two agree where the explicit one is still valid"
      (let [a (energy-ratio 10.0 :absolute 0.005 4.0)
            b (energy-ratio 10.0 :semi-implicit 0.005 4.0)]
        (is (< (abs (- a b)) 0.1))))

    (testing "with no rotation the two schemes describe the same physics"
      ;; Not bit-for-bit -- they reach it by different arithmetic -- but
      ;; the implicit operator degenerates to the Poisson one when the
      ;; rotation is zero, so they should differ only by round-off.
      (let [field (fn [cor]
                    (let [sim (jets 0.0 cor 32 64)]
                      (vec (:vorticity (reduce (fn [s _] (sf/step! s 0.01))
                                               sim (range 50))))))
            a (field :absolute) b (field :semi-implicit)
            scale (apply max (map abs a))]
        (is (< (apply max (map #(abs (- (double %1) (double %2))) a b))
               (* 1e-10 scale)))))))

(deftest zonal-jets-test
  (let [sim (sf/simulation {:grid {:nlat 48 :nlon 128} :tracer? false})
        g (:grid sim)]
    (sf/zonal-jets! g (:vorticity sim) {:n-jets 9 :amplitude 8.0 :perturbation 0.0})

    (testing "with no perturbation the vorticity depends on latitude alone"
      (let [^doubles w (:vorticity sim)]
        (doseq [j (range 48)]
          (let [row (map #(aget w (+ (* j 128) %)) (range 128))]
            (is (< (- (apply max row) (apply min row)) 1e-12))))))

    (testing "and the flow it makes is purely zonal"
      (sf/update-velocity! sim)
      (is (> (zonal-energy-fraction sim) 0.999))
      (is (< (max-abs (:north sim)) 1e-9)))

    (testing "the jets alternate, and there are as many as asked for"
      (let [^doubles w (:vorticity sim)
            ^doubles phi (:phi g)
            column (map #(aget w (* % 128)) (range 48))
            crossings (count (filter (fn [[a b]] (neg? (* a b)))
                                     (partition 2 1 column)))]
        ;; sin(9 phi) vanishes at phi = k pi/9, and nine of those fall
        ;; strictly inside the grid's latitude range.
        (is (= 9 crossings))
        (is (pos? (aget phi 47)))))

    (testing "a perturbation breaks the symmetry without swamping the jets"
      (let [sim2 (sf/simulation {:grid {:nlat 48 :nlon 128} :tracer? false})]
        (sf/zonal-jets! (:grid sim2) (:vorticity sim2) {:perturbation 0.03})
        (sf/update-velocity! sim2)
        (is (> (zonal-energy-fraction sim2) 0.99))
        (is (pos? (max-abs (:north sim2))))))))
