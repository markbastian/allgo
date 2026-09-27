(ns allgo.astro.sgp4
  "SGP4 and SDP4, the propagators two-line element sets are made for.

  Sources, which are authoritative here wherever anything else differs:

  - Hoots, F. R. and R. L. Roehrich, *Spacetrack Report No. 3: Models for
    Propagation of NORAD Element Sets*, 1980.
  - Vallado, D. A., P. Crawford, R. Hujsak and T. S. Kelso, *Revisiting
    Spacetrack Report #3*, AIAA 2006-6753, and its revisions 1-3, with the
    code that accompanies them:
    https://celestrak.org/publications/AIAA/2006-6753/

  This namespace is a transcription of that code's C++ version (SGP4.cpp,
  the 12 March 2020 revision, which the paper's FAQ notes is the code in
  STK), function for function and expression for expression: the names
  are the C++ names and each expression keeps the C++ order of
  operations, so results agree with it to the last few bits rather than
  merely to the theory. C's `fmod` is an exact
  `fmod` here, Clojure's `rem` on doubles being inexact. The one structural change is
  that `sgp4` is a pure function: where the C++ keeps the resonance
  integrator's last step in the element record to resume from, this
  integrates from the epoch each call, which lands on the same half-day
  grid and gives the same numbers.

  A TLE is not a state vector nor an osculating element set: its mean
  elements are the output of a fit made with exactly this theory and mean
  nothing without it. Output is in TEME, kilometers and kilometers per
  second, at `tsince` minutes from the element set's epoch;
  `allgo.astro.reduction/teme->eci` takes it to GCRF.

  `opsmode` is :a, AFSPC's operation, or :i, the paper's improved mode.
  :a is the default, as the paper's FAQ recommends to match the US Air
  Force, and is the mode the published C++ verification output is in."
  (:require [clojure.math :as math]
            [clojure.string :as str]))

(def ^:private pi math/PI)
(def ^:private twopi (* 2.0 pi))
(def ^:private x2o3 (/ 2.0 3.0))
(def ^:private temp4 1.5e-12)

(defn- fmod
  "C's fmod, which is exact: the remainder of x/y with the sign of x.
  Clojure's `rem` on doubles is x - trunc(x/y) y, which rounds."
  [x y]
  #?(:clj (let [r (math/IEEE-remainder x y)]
            (cond (zero? r) (if (neg? x) -0.0 0.0)
                  (= (neg? r) (neg? x)) r
                  (neg? x) (- r (abs y))
                  :else (+ r (abs y))))
     :cljs (js-mod x y)))

;; ------------------------------------------------------------ dpper

(defn- dpper
  "Lunar-solar periodics: [ep inclp nodep argpp mp]."
  [{:keys [e3 ee2 peo pgho pho pinco plo se2 se3 sgh2 sgh3 sgh4 sh2 sh3 si2
           si3 sl2 sl3 sl4 xgh2 xgh3 xgh4 xh2 xh3 xi2 xi3 xl2 xl3 xl4 zmol zmos]}
   t init ep inclp nodep argpp mp opsmode]
  (let [zns 1.19459e-5
        zes 0.01675
        znl 1.5835218e-4
        zel 0.05490
        zm (if (= init \y) zmos (+ zmos (* zns t)))
        zf (+ zm (* 2.0 zes (math/sin zm)))
        sinzf (math/sin zf)
        f2 (- (* 0.5 sinzf sinzf) 0.25)
        f3 (* -0.5 sinzf (math/cos zf))
        ses (+ (* se2 f2) (* se3 f3))
        sis (+ (* si2 f2) (* si3 f3))
        sls (+ (* sl2 f2) (* sl3 f3) (* sl4 sinzf))
        sghs (+ (* sgh2 f2) (* sgh3 f3) (* sgh4 sinzf))
        shs (+ (* sh2 f2) (* sh3 f3))
        zm (if (= init \y) zmol (+ zmol (* znl t)))
        zf (+ zm (* 2.0 zel (math/sin zm)))
        sinzf (math/sin zf)
        f2 (- (* 0.5 sinzf sinzf) 0.25)
        f3 (* -0.5 sinzf (math/cos zf))
        sel (+ (* ee2 f2) (* e3 f3))
        sil (+ (* xi2 f2) (* xi3 f3))
        sll (+ (* xl2 f2) (* xl3 f3) (* xl4 sinzf))
        sghl (+ (* xgh2 f2) (* xgh3 f3) (* xgh4 sinzf))
        shll (+ (* xh2 f2) (* xh3 f3))
        pe (+ ses sel)
        pinc (+ sis sil)
        pl (+ sls sll)
        pgh (+ sghs sghl)
        ph (+ shs shll)]
    (if (not= init \n)
      [ep inclp nodep argpp mp]
      (let [pe (- pe peo)
            pinc (- pinc pinco)
            pl (- pl plo)
            pgh (- pgh pgho)
            ph (- ph pho)
            inclp (+ inclp pinc)
            ep (+ ep pe)
            sinip (math/sin inclp)
            cosip (math/cos inclp)]
        (if (>= inclp 0.2)
          (let [ph (/ ph sinip)
                pgh (- pgh (* cosip ph))]
            [ep inclp (+ nodep ph) (+ argpp pgh) (+ mp pl)])
          ;; apply periodics with the Lyddane modification
          (let [sinop (math/sin nodep)
                cosop (math/cos nodep)
                alfdp (* sinip sinop)
                betdp (* sinip cosop)
                dalf (+ (* ph cosop) (* pinc cosip sinop))
                dbet (+ (* (- ph) sinop) (* pinc cosip cosop))
                alfdp (+ alfdp dalf)
                betdp (+ betdp dbet)
                nodep (fmod nodep twopi)
                nodep (if (and (< nodep 0.0) (= opsmode \a)) (+ nodep twopi) nodep)
                xls (+ mp argpp (* cosip nodep))
                dls (- (+ pl pgh) (* pinc nodep sinip))
                xls (+ xls dls)
                xnoh nodep
                nodep (math/atan2 alfdp betdp)
                nodep (if (and (< nodep 0.0) (= opsmode \a)) (+ nodep twopi) nodep)
                nodep (if (> (abs (- xnoh nodep)) pi)
                        (if (< nodep xnoh) (+ nodep twopi) (- nodep twopi))
                        nodep)
                mp (+ mp pl)
                argpp (- xls mp (* cosip nodep))]
            [ep inclp nodep argpp mp]))))))

;; ------------------------------------------------------------ dscom

(defn- dscom
  "Deep-space common items, from the Sun's and the Moon's geometry."
  [epoch ep argpp tc inclp nodep np]
  (let [zes 0.01675
        zel 0.05490
        c1ss 2.9864797e-6
        c1l 4.7968065e-7
        zsinis 0.39785416
        zcosis 0.91744867
        zcosgs 0.1945905
        zsings -0.98088458
        nm np
        em ep
        snodm (math/sin nodep)
        cnodm (math/cos nodep)
        sinomm (math/sin argpp)
        cosomm (math/cos argpp)
        sinim (math/sin inclp)
        cosim (math/cos inclp)
        emsq (* em em)
        betasq (- 1.0 emsq)
        rtemsq (math/sqrt betasq)
        day (+ epoch 18261.5 (/ tc 1440.0))
        xnodce (fmod (- 4.5236020 (* 9.2422029e-4 day)) twopi)
        stem (math/sin xnodce)
        ctem (math/cos xnodce)
        zcosil (- 0.91375164 (* 0.03568096 ctem))
        zsinil (math/sqrt (- 1.0 (* zcosil zcosil)))
        zsinhl (/ (* 0.089683511 stem) zsinil)
        zcoshl (math/sqrt (- 1.0 (* zsinhl zsinhl)))
        gam (+ 5.8351514 (* 0.0019443680 day))
        zx (/ (* 0.39785416 stem) zsinil)
        zy (+ (* zcoshl ctem) (* 0.91744867 zsinhl stem))
        zx (math/atan2 zx zy)
        zx (- (+ gam zx) xnodce)
        zcosgl (math/cos zx)
        zsingl (math/sin zx)
        xnoi (/ 1.0 nm)
        ;; one pass of the C++ loop: lsflg 1 is the Sun, 2 the Moon
        pass (fn [zcosg zsing zcosi zsini zcosh zsinh cc]
               (let [a1 (+ (* zcosg zcosh) (* zsing zcosi zsinh))
                     a3 (+ (* (- zsing) zcosh) (* zcosg zcosi zsinh))
                     a7 (+ (* (- zcosg) zsinh) (* zsing zcosi zcosh))
                     a8 (* zsing zsini)
                     a9 (+ (* zsing zsinh) (* zcosg zcosi zcosh))
                     a10 (* zcosg zsini)
                     a2 (+ (* cosim a7) (* sinim a8))
                     a4 (+ (* cosim a9) (* sinim a10))
                     a5 (+ (* (- sinim) a7) (* cosim a8))
                     a6 (+ (* (- sinim) a9) (* cosim a10))
                     x1 (+ (* a1 cosomm) (* a2 sinomm))
                     x2 (+ (* a3 cosomm) (* a4 sinomm))
                     x3 (+ (* (- a1) sinomm) (* a2 cosomm))
                     x4 (+ (* (- a3) sinomm) (* a4 cosomm))
                     x5 (* a5 sinomm)
                     x6 (* a6 sinomm)
                     x7 (* a5 cosomm)
                     x8 (* a6 cosomm)
                     z31 (- (* 12.0 x1 x1) (* 3.0 x3 x3))
                     z32 (- (* 24.0 x1 x2) (* 6.0 x3 x4))
                     z33 (- (* 12.0 x2 x2) (* 3.0 x4 x4))
                     z1 (+ (* 3.0 (+ (* a1 a1) (* a2 a2))) (* z31 emsq))
                     z2 (+ (* 6.0 (+ (* a1 a3) (* a2 a4))) (* z32 emsq))
                     z3 (+ (* 3.0 (+ (* a3 a3) (* a4 a4))) (* z33 emsq))
                     z11 (+ (* -6.0 a1 a5) (* emsq (- (* -24.0 x1 x7) (* 6.0 x3 x5))))
                     z12 (+ (* -6.0 (+ (* a1 a6) (* a3 a5)))
                            (* emsq (- (* -24.0 (+ (* x2 x7) (* x1 x8)))
                                       (* 6.0 (+ (* x3 x6) (* x4 x5))))))
                     z13 (+ (* -6.0 a3 a6) (* emsq (- (* -24.0 x2 x8) (* 6.0 x4 x6))))
                     z21 (+ (* 6.0 a2 a5) (* emsq (- (* 24.0 x1 x5) (* 6.0 x3 x7))))
                     z22 (+ (* 6.0 (+ (* a4 a5) (* a2 a6)))
                            (* emsq (- (* 24.0 (+ (* x2 x5) (* x1 x6)))
                                       (* 6.0 (+ (* x4 x7) (* x3 x8))))))
                     z23 (+ (* 6.0 a4 a6) (* emsq (- (* 24.0 x2 x6) (* 6.0 x4 x8))))
                     z1 (+ z1 z1 (* betasq z31))
                     z2 (+ z2 z2 (* betasq z32))
                     z3 (+ z3 z3 (* betasq z33))
                     s3 (* cc xnoi)
                     s2 (/ (* -0.5 s3) rtemsq)
                     s4 (* s3 rtemsq)
                     s1 (* -15.0 em s4)]
                 {:s1 s1 :s2 s2 :s3 s3 :s4 s4
                  :s5 (+ (* x1 x3) (* x2 x4))
                  :s6 (+ (* x2 x3) (* x1 x4))
                  :s7 (- (* x2 x4) (* x1 x3))
                  :z1 z1 :z2 z2 :z3 z3 :z11 z11 :z12 z12 :z13 z13
                  :z21 z21 :z22 z22 :z23 z23 :z31 z31 :z32 z32 :z33 z33}))
        {ss1 :s1 ss2 :s2 ss3 :s3 ss4 :s4 ss5 :s5 ss6 :s6 ss7 :s7
         sz1 :z1 sz2 :z2 sz3 :z3 sz11 :z11 sz12 :z12 sz13 :z13
         sz21 :z21 sz22 :z22 sz23 :z23 sz31 :z31 sz32 :z32 sz33 :z33}
        (pass zcosgs zsings zcosis zsinis cnodm snodm c1ss)
        {:keys [s1 s2 s3 s4 s5 s6 s7 z1 z2 z3 z11 z12 z13 z21 z22 z23 z31 z32 z33]}
        (pass zcosgl zsingl zcosil zsinil
              (+ (* zcoshl cnodm) (* zsinhl snodm))
              (- (* snodm zcoshl) (* cnodm zsinhl))
              c1l)]
    {:snodm snodm :cnodm cnodm :sinim sinim :cosim cosim :sinomm sinomm :cosomm cosomm
     :day day :em em :emsq emsq :gam gam :rtemsq rtemsq :nm nm
     :peo 0.0 :pinco 0.0 :plo 0.0 :pgho 0.0 :pho 0.0
     :s1 s1 :s2 s2 :s3 s3 :s4 s4 :s5 s5 :s6 s6 :s7 s7
     :ss1 ss1 :ss2 ss2 :ss3 ss3 :ss4 ss4 :ss5 ss5 :ss6 ss6 :ss7 ss7
     :sz1 sz1 :sz2 sz2 :sz3 sz3 :sz11 sz11 :sz12 sz12 :sz13 sz13
     :sz21 sz21 :sz22 sz22 :sz23 sz23 :sz31 sz31 :sz32 sz32 :sz33 sz33
     :z1 z1 :z2 z2 :z3 z3 :z11 z11 :z12 z12 :z13 z13
     :z21 z21 :z22 z22 :z23 z23 :z31 z31 :z32 z32 :z33 z33
     :zmol (fmod (- (+ 4.7199672 (* 0.22997150 day)) gam) twopi)
     :zmos (fmod (+ 6.2565837 (* 0.017201977 day)) twopi)
     :se2 (* 2.0 ss1 ss6)
     :se3 (* 2.0 ss1 ss7)
     :si2 (* 2.0 ss2 sz12)
     :si3 (* 2.0 ss2 (- sz13 sz11))
     :sl2 (* -2.0 ss3 sz2)
     :sl3 (* -2.0 ss3 (- sz3 sz1))
     :sl4 (* -2.0 ss3 (- -21.0 (* 9.0 emsq)) zes)
     :sgh2 (* 2.0 ss4 sz32)
     :sgh3 (* 2.0 ss4 (- sz33 sz31))
     :sgh4 (* -18.0 ss4 zes)
     :sh2 (* -2.0 ss2 sz22)
     :sh3 (* -2.0 ss2 (- sz23 sz21))
     :ee2 (* 2.0 s1 s6)
     :e3 (* 2.0 s1 s7)
     :xi2 (* 2.0 s2 z12)
     :xi3 (* 2.0 s2 (- z13 z11))
     :xl2 (* -2.0 s3 z2)
     :xl3 (* -2.0 s3 (- z3 z1))
     :xl4 (* -2.0 s3 (- -21.0 (* 9.0 emsq)) zel)
     :xgh2 (* 2.0 s4 z32)
     :xgh3 (* 2.0 s4 (- z33 z31))
     :xgh4 (* -18.0 s4 zel)
     :xh2 (* -2.0 s2 z22)
     :xh3 (* -2.0 s2 (- z23 z21))}))

;; ------------------------------------------------------------ dsinit

(defn- dsinit
  "Deep-space initialization: the Sun's and Moon's secular rates and, for
  12- and 24-hour orbits, the resonance terms. Returns the satrec fields
  it sets."
  [xke {:keys [cosim emsq s1 s2 s3 s4 s5 sinim ss1 ss2 ss3 ss4 ss5 sz1 sz3 sz11
               sz13 sz21 sz23 sz31 sz33 z1 z3 z11 z13 z21 z23 z31 z33]}
   argpo tc gsto mo mdot no nodeo nodedot xpidot ecco eccsq em nm inclm]
  (let [q22 1.7891679e-6
        q31 2.1460748e-6
        q33 2.2123015e-7
        root22 1.7891679e-6
        root44 7.3636953e-9
        root54 2.1765803e-9
        rptim 4.37526908801129966e-3
        root32 3.7393792e-7
        root52 1.1428639e-7
        znl 1.5835218e-4
        zns 1.19459e-5
        ;; deep space initialization
        irez (cond (and (>= nm 8.26e-3) (<= nm 9.24e-3) (>= em 0.5)) 2
                   (and (< nm 0.0052359877) (> nm 0.0034906585)) 1
                   :else 0)
        polar? (or (< inclm 5.2359877e-2) (> inclm (- pi 5.2359877e-2)))
        ;; solar terms
        ses (* ss1 zns ss5)
        sis (* ss2 zns (+ sz11 sz13))
        sls (* (- zns) ss3 (- (+ sz1 sz3) 14.0 (* 6.0 emsq)))
        sghs (* ss4 zns (- (+ sz31 sz33) 6.0))
        shs (* (- zns) ss2 (+ sz21 sz23))
        shs (if polar? 0.0 shs)
        shs (if (not= sinim 0.0) (/ shs sinim) shs)
        sgs (- sghs (* cosim shs))
        ;; lunar terms
        dedt (+ ses (* s1 znl s5))
        didt (+ sis (* s2 znl (+ z11 z13)))
        dmdt (- sls (* znl s3 (- (+ z1 z3) 14.0 (* 6.0 emsq))))
        sghl (* s4 znl (- (+ z31 z33) 6.0))
        shll (* (- znl) s2 (+ z21 z23))
        shll (if polar? 0.0 shll)
        domdt (+ sgs sghl)
        dnodt shs
        [domdt dnodt] (if (not= sinim 0.0)
                        [(- domdt (* (/ cosim sinim) shll)) (+ dnodt (/ shll sinim))]
                        [domdt dnodt])
        theta (fmod (+ gsto (* tc rptim)) twopi)
        rates {:irez irez :dedt dedt :didt didt :dmdt dmdt :dnodt dnodt :domdt domdt}]
    (if (zero? irez)
      rates
      (let [aonv (math/pow (/ nm xke) x2o3)]
        (merge
         rates
         (if (= irez 2)
           ;; geopotential resonance for 12-hour orbits
           (let [cosisq (* cosim cosim)
                 em ecco
                 emsq eccsq
                 eoc (* em emsq)
                 g201 (- -0.306 (* (- em 0.64) 0.440))
                 [g211 g310 g322 g410 g422 g520]
                 (if (<= em 0.65)
                   [(+ (- 3.616 (* 13.2470 em)) (* 16.2900 emsq))
                    (+ (- (+ -19.302 (* 117.3900 em)) (* 228.4190 emsq)) (* 156.5910 eoc))
                    (+ (- (+ -18.9068 (* 109.7927 em)) (* 214.6334 emsq)) (* 146.5816 eoc))
                    (+ (- (+ -41.122 (* 242.6940 em)) (* 471.0940 emsq)) (* 313.9530 eoc))
                    (+ (- (+ -146.407 (* 841.8800 em)) (* 1629.014 emsq)) (* 1083.4350 eoc))
                    (+ (- (+ -532.114 (* 3017.977 em)) (* 5740.032 emsq)) (* 3708.2760 eoc))]
                   [(+ (- (+ -72.099 (* 331.819 em)) (* 508.738 emsq)) (* 266.724 eoc))
                    (+ (- (+ -346.844 (* 1582.851 em)) (* 2415.925 emsq)) (* 1246.113 eoc))
                    (+ (- (+ -342.585 (* 1554.908 em)) (* 2366.899 emsq)) (* 1215.972 eoc))
                    (+ (- (+ -1052.797 (* 4758.686 em)) (* 7193.992 emsq)) (* 3651.957 eoc))
                    (+ (- (+ -3581.690 (* 16178.110 em)) (* 24462.770 emsq)) (* 12422.520 eoc))
                    (if (> em 0.715)
                      (+ (- (+ -5149.66 (* 29936.92 em)) (* 54087.36 emsq)) (* 31324.56 eoc))
                      (+ (- 1464.74 (* 4664.75 em)) (* 3763.64 emsq)))])
                 [g533 g521 g532]
                 (if (< em 0.7)
                   [(+ (- (+ -919.22770 (* 4988.6100 em)) (* 9064.7700 emsq)) (* 5542.21 eoc))
                    (+ (- (+ -822.71072 (* 4568.6173 em)) (* 8491.4146 emsq)) (* 5337.524 eoc))
                    (+ (- (+ -853.66600 (* 4690.2500 em)) (* 8624.7700 emsq)) (* 5341.4 eoc))]
                   [(+ (- (+ -37995.780 (* 161616.52 em)) (* 229838.20 emsq)) (* 109377.94 eoc))
                    (+ (- (+ -51752.104 (* 218913.95 em)) (* 309468.16 emsq)) (* 146349.42 eoc))
                    (+ (- (+ -40023.880 (* 170470.89 em)) (* 242699.48 emsq)) (* 115605.82 eoc))])
                 sini2 (* sinim sinim)
                 f220 (* 0.75 (+ 1.0 (* 2.0 cosim) cosisq))
                 f221 (* 1.5 sini2)
                 f321 (* 1.875 sinim (- 1.0 (* 2.0 cosim) (* 3.0 cosisq)))
                 f322 (* -1.875 sinim (- (+ 1.0 (* 2.0 cosim)) (* 3.0 cosisq)))
                 f441 (* 35.0 sini2 f220)
                 f442 (* 39.3750 sini2 sini2)
                 f522 (* 9.84375 sinim
                         (+ (* sini2 (- 1.0 (* 2.0 cosim) (* 5.0 cosisq)))
                            (* 0.33333333 (+ (+ -2.0 (* 4.0 cosim)) (* 6.0 cosisq)))))
                 f523 (* sinim
                         (+ (* 4.92187512 sini2 (+ (- -2.0 (* 4.0 cosim)) (* 10.0 cosisq)))
                            (* 6.56250012 (- (+ 1.0 (* 2.0 cosim)) (* 3.0 cosisq)))))
                 f542 (* 29.53125 sinim
                         (+ (- 2.0 (* 8.0 cosim))
                            (* cosisq (+ (+ -12.0 (* 8.0 cosim)) (* 10.0 cosisq)))))
                 f543 (* 29.53125 sinim
                         (+ (- -2.0 (* 8.0 cosim))
                            (* cosisq (- (+ 12.0 (* 8.0 cosim)) (* 10.0 cosisq)))))
                 xno2 (* nm nm)
                 ainv2 (* aonv aonv)
                 temp1 (* 3.0 xno2 ainv2)
                 temp (* temp1 root22)
                 d2201 (* temp f220 g201)
                 d2211 (* temp f221 g211)
                 temp1 (* temp1 aonv)
                 temp (* temp1 root32)
                 d3210 (* temp f321 g310)
                 d3222 (* temp f322 g322)
                 temp1 (* temp1 aonv)
                 temp (* 2.0 temp1 root44)
                 d4410 (* temp f441 g410)
                 d4422 (* temp f442 g422)
                 temp1 (* temp1 aonv)
                 temp (* temp1 root52)
                 d5220 (* temp f522 g520)
                 d5232 (* temp f523 g532)
                 temp (* 2.0 temp1 root54)
                 d5421 (* temp f542 g521)
                 d5433 (* temp f543 g533)]
             {:d2201 d2201 :d2211 d2211 :d3210 d3210 :d3222 d3222
              :d4410 d4410 :d4422 d4422 :d5220 d5220 :d5232 d5232
              :d5421 d5421 :d5433 d5433
              :xlamo (fmod (- (+ mo nodeo nodeo) theta theta) twopi)
              :xfact (- (+ mdot dmdt (* 2.0 (- (+ nodedot dnodt) rptim))) no)})
           ;; synchronous resonance terms
           (let [g200 (+ 1.0 (* emsq (+ -2.5 (* 0.8125 emsq))))
                 g310 (+ 1.0 (* 2.0 emsq))
                 g300 (+ 1.0 (* emsq (+ -6.0 (* 6.60937 emsq))))
                 f220 (* 0.75 (+ 1.0 cosim) (+ 1.0 cosim))
                 f311 (- (* 0.9375 sinim sinim (+ 1.0 (* 3.0 cosim))) (* 0.75 (+ 1.0 cosim)))
                 f330 (+ 1.0 cosim)
                 f330 (* 1.875 f330 f330 f330)
                 del1 (* 3.0 nm nm aonv aonv)
                 del2 (* 2.0 del1 f220 g200 q22)
                 del3 (* 3.0 del1 f330 g300 q33 aonv)
                 del1 (* del1 f311 g310 q31 aonv)]
             {:del1 del1 :del2 del2 :del3 del3
              :xlamo (fmod (- (+ mo nodeo argpo) theta) twopi)
              :xfact (- (+ (- (+ mdot xpidot) rptim) dmdt domdt dnodt) no)})))))))

;; ------------------------------------------------------------ dspace

(defn- dspace
  "Deep-space contributions to the mean elements at `t` minutes: the
  secular rates and the resonance integration, in half-day steps from the
  epoch. [em argpm inclm mm nodem nm]."
  [{:keys [irez d2201 d2211 d3210 d3222 d4410 d4422 d5220 d5232 d5421 d5433
           dedt del1 del2 del3 didt dmdt dnodt domdt argpo argpdot gsto xfact
           xlamo no-unkozai]}
   t tc em argpm inclm mm nodem nm]
  (let [no no-unkozai
        fasx2 0.13130908
        fasx4 2.8843198
        fasx6 0.37448087
        g22 5.7686396
        g32 0.95240898
        g44 1.8014998
        g52 1.0508330
        g54 4.4108898
        rptim 4.37526908801129966e-3
        stepp 720.0
        stepn -720.0
        step2 259200.0
        ;; calculate deep space resonance effects
        theta (fmod (+ gsto (* tc rptim)) twopi)
        em (+ em (* dedt t))
        inclm (+ inclm (* didt t))
        argpm (+ argpm (* domdt t))
        nodem (+ nodem (* dnodt t))
        mm (+ mm (* dmdt t))]
    (if (zero? irez)
      [em argpm inclm mm nodem nm]
      (let [delt (if (> t 0.0) stepp stepn)
            sn math/sin
            cs math/cos
            ;; dot terms
            derivs
            (fn [xli xni atime]
              (if (not= irez 2)
                ;; near-synchronous resonance terms
                (let [xndt (+ (* del1 (sn (- xli fasx2)))
                              (* del2 (sn (* 2.0 (- xli fasx4))))
                              (* del3 (sn (* 3.0 (- xli fasx6)))))
                      xldot (+ xni xfact)
                      xnddt (+ (* del1 (cs (- xli fasx2)))
                               (* 2.0 del2 (cs (* 2.0 (- xli fasx4))))
                               (* 3.0 del3 (cs (* 3.0 (- xli fasx6)))))]
                  [xndt xldot (* xnddt xldot)])
                ;; near-half-day resonance terms
                (let [xomi (+ argpo (* argpdot atime))
                      x2omi (+ xomi xomi)
                      x2li (+ xli xli)
                      xndt (+ (* d2201 (sn (- (+ x2omi xli) g22)))
                              (* d2211 (sn (- xli g22)))
                              (* d3210 (sn (- (+ xomi xli) g32)))
                              (* d3222 (sn (- (+ (- xomi) xli) g32)))
                              (* d4410 (sn (- (+ x2omi x2li) g44)))
                              (* d4422 (sn (- x2li g44)))
                              (* d5220 (sn (- (+ xomi xli) g52)))
                              (* d5232 (sn (- (+ (- xomi) xli) g52)))
                              (* d5421 (sn (- (+ xomi x2li) g54)))
                              (* d5433 (sn (- (+ (- xomi) x2li) g54))))
                      xldot (+ xni xfact)
                      xnddt (+ (* d2201 (cs (- (+ x2omi xli) g22)))
                               (* d2211 (cs (- xli g22)))
                               (* d3210 (cs (- (+ xomi xli) g32)))
                               (* d3222 (cs (- (+ (- xomi) xli) g32)))
                               (* d5220 (cs (- (+ xomi xli) g52)))
                               (* d5232 (cs (- (+ (- xomi) xli) g52)))
                               (* 2.0 (+ (* d4410 (cs (- (+ x2omi x2li) g44)))
                                         (* d4422 (cs (- x2li g44)))
                                         (* d5421 (cs (- (+ xomi x2li) g54)))
                                         (* d5433 (cs (- (+ (- xomi) x2li) g54))))))]
                  [xndt xldot (* xnddt xldot)])))
            ;; the integrator, from the epoch -- the state the C++ resets
            ;; to whenever it cannot resume from its last call
            [xli xni xndt xldot xnddt ft]
            (loop [xli xlamo xni no atime 0.0]
              (let [[xndt xldot xnddt] (derivs xli xni atime)]
                (if (>= (abs (- t atime)) stepp)
                  (recur (+ xli (* xldot delt) (* xndt step2))
                         (+ xni (* xndt delt) (* xnddt step2))
                         (+ atime delt))
                  [xli xni xndt xldot xnddt (- t atime)])))
            nm (+ xni (* xndt ft) (* xnddt ft ft 0.5))
            xl (+ xli (* xldot ft) (* xndt ft ft 0.5))
            mm (if (not= irez 1)
                 (+ (- xl (* 2.0 nodem)) (* 2.0 theta))
                 (+ (- xl nodem argpm) theta))
            dndt (- nm no)
            nm (+ no dndt)]
        [em argpm inclm mm nodem nm]))))

;; ------------------------------------------------------------- initl

(defn gstime
  "Greenwich mean sidereal time, radians, at the UT1 Julian Date `jdut1`
  (IAU-82), as SGP4 computes it."
  [jdut1]
  (let [deg2rad (/ pi 180.0)
        tut1 (/ (- jdut1 2451545.0) 36525.0)
        temp (+ (* -6.2e-6 tut1 tut1 tut1) (* 0.093104 tut1 tut1)
                (* (+ (* 876600.0 3600) 8640184.812866) tut1) 67310.54841)
        temp (fmod (/ (* temp deg2rad) 240.0) twopi)]
    (if (< temp 0.0) (+ temp twopi) temp)))

(defn- initl [xke j2 ecco epoch inclo no-kozai]
  (let [eccsq (* ecco ecco)
        omeosq (- 1.0 eccsq)
        rteosq (math/sqrt omeosq)
        cosio (math/cos inclo)
        cosio2 (* cosio cosio)
        ;; un-kozai the mean motion
        ak (math/pow (/ xke no-kozai) x2o3)
        d1 (/ (* 0.75 j2 (- (* 3.0 cosio2) 1.0)) (* rteosq omeosq))
        del (/ d1 (* ak ak))
        adel (* ak (- 1.0 (* del del) (* del (+ (/ 1.0 3.0) (/ (* 134.0 del del) 81.0)))))
        del (/ d1 (* adel adel))
        no-unkozai (/ no-kozai (+ 1.0 del))
        ao (math/pow (/ xke no-unkozai) x2o3)
        sinio (math/sin inclo)
        po (* ao omeosq)
        con42 (- 1.0 (* 5.0 cosio2))
        con41 (- (- con42) cosio2 cosio2)]
    ;; the C++ also works out AFSPC's 1970-based sidereal time here but has
    ;; used gstime in both modes since 2008
    {:ao ao :con41 con41 :con42 con42 :cosio cosio :cosio2 cosio2 :eccsq eccsq
     :omeosq omeosq :posq (* po po) :rp (* ao (- 1.0 ecco)) :rteosq rteosq :sinio sinio
     :gsto (gstime (+ epoch 2433281.5)) :no-unkozai no-unkozai}))

;; --------------------------------------------------------- constants

(defn getgravconst
  "The constants SGP4 may be run with: :wgs72old, :wgs72 (what TLEs are
  fitted with, and the default) or :wgs84."
  [whichconst]
  (let [[mus radiusearthkm xke j2 j3 j4]
        (case whichconst
          :wgs72old [398600.79964 6378.135 0.0743669161 0.001082616 -0.00000253881 -0.00000165597]
          :wgs72 (let [mus 398600.8 re 6378.135]
                   [mus re (/ 60.0 (math/sqrt (/ (* re re re) mus)))
                    0.001082616 -0.00000253881 -0.00000165597])
          :wgs84 (let [mus 398600.5 re 6378.137]
                   [mus re (/ 60.0 (math/sqrt (/ (* re re re) mus)))
                    0.00108262998905 -0.00000253215306 -0.00000161098761]))]
    {:tumin (/ 1.0 xke) :mus mus :radiusearthkm radiusearthkm :xke xke
     :j2 j2 :j3 j3 :j4 j4 :j3oj2 (/ j3 j2)}))

;; ---------------------------------------------------------- sgp4init

(defn sgp4init
  "The element record `sgp4` propagates, from mean elements: `:epoch`
  (days from 1949 December 31 0h UT), `:bstar`, `:ndot`, `:nddot`,
  `:ecco`, `:argpo`, `:inclo`, `:mo`, `:no-kozai` (radians and radians per
  minute) and `:nodeo`. Options `:whichconst` (default :wgs72) and
  `:opsmode` (default :a)."
  [{:keys [epoch bstar ecco argpo inclo mo no-kozai nodeo whichconst opsmode]
    :or {whichconst :wgs72 opsmode :a} :as elements}]
  (let [{:keys [radiusearthkm xke j2 j3oj2 j4 tumin] :as grav} (getgravconst whichconst)
        opsmode (case opsmode (:a \a) \a (:i \i) \i)
        ss (+ (/ 78.0 radiusearthkm) 1.0)
        qzms2ttemp (/ (- 120.0 78.0) radiusearthkm)
        qzms2t (* qzms2ttemp qzms2ttemp qzms2ttemp qzms2ttemp)
        {:keys [ao con41 con42 cosio cosio2 eccsq omeosq posq rp rteosq sinio gsto no-unkozai]}
        (initl xke j2 ecco epoch inclo no-kozai)
        a (math/pow (* no-unkozai tumin) (/ -2.0 3.0))
        satrec (merge
                (dissoc elements :whichconst)
                grav
                {:opsmode opsmode :method \n :isimp 0 :irez 0
                 :aycof 0.0 :con41 con41 :cc1 0.0 :cc4 0.0 :cc5 0.0 :d2 0.0 :d3 0.0 :d4 0.0
                 :delmo 0.0 :eta 0.0 :argpdot 0.0 :omgcof 0.0 :sinmao 0.0 :t2cof 0.0
                 :t3cof 0.0 :t4cof 0.0 :t5cof 0.0 :x1mth2 0.0 :x7thm1 0.0 :mdot 0.0
                 :nodedot 0.0 :xlcof 0.0 :xmcof 0.0 :nodecf 0.0
                 :gsto gsto :no-unkozai no-unkozai :a a
                 :alta (- (* a (+ 1.0 ecco)) 1.0)
                 :altp (- (* a (- 1.0 ecco)) 1.0)})]
    (if-not (or (>= omeosq 0.0) (>= no-unkozai 0.0))
      satrec
      (let [isimp (if (< rp (+ (/ 220.0 radiusearthkm) 1.0)) 1 0)
            perige (* (- rp 1.0) radiusearthkm)
            ;; for perigees below 156 km, s and qoms2t are altered
            [sfour qzms24]
            (if (< perige 156.0)
              (let [sfour (if (< perige 98.0) 20.0 (- perige 78.0))
                    qzms24temp (/ (- 120.0 sfour) radiusearthkm)]
                [(+ (/ sfour radiusearthkm) 1.0)
                 (* qzms24temp qzms24temp qzms24temp qzms24temp)])
              [ss qzms2t])
            pinvsq (/ 1.0 posq)
            tsi (/ 1.0 (- ao sfour))
            eta (* ao ecco tsi)
            etasq (* eta eta)
            eeta (* ecco eta)
            psisq (abs (- 1.0 etasq))
            coef (* qzms24 (math/pow tsi 4.0))
            coef1 (/ coef (math/pow psisq 3.5))
            cc2 (* coef1 no-unkozai
                   (+ (* ao (+ 1.0 (* 1.5 etasq) (* eeta (+ 4.0 etasq))))
                      (* (/ (* 0.375 j2 tsi) psisq) con41
                         (+ 8.0 (* 3.0 etasq (+ 8.0 etasq))))))
            cc1 (* bstar cc2)
            cc3 (if (> ecco 1.0e-4)
                  (/ (* -2.0 coef tsi j3oj2 no-unkozai sinio) ecco)
                  0.0)
            x1mth2 (- 1.0 cosio2)
            cc4 (* 2.0 no-unkozai coef1 ao omeosq
                   (- (+ (* eta (+ 2.0 (* 0.5 etasq)))
                         (* ecco (+ 0.5 (* 2.0 etasq))))
                      (* (/ (* j2 tsi) (* ao psisq))
                         (+ (* -3.0 con41 (+ (- 1.0 (* 2.0 eeta))
                                             (* etasq (- 1.5 (* 0.5 eeta)))))
                            (* 0.75 x1mth2 (- (* 2.0 etasq) (* eeta (+ 1.0 etasq)))
                               (math/cos (* 2.0 argpo)))))))
            cc5 (* 2.0 coef1 ao omeosq (+ 1.0 (* 2.75 (+ etasq eeta)) (* eeta etasq)))
            cosio4 (* cosio2 cosio2)
            temp1 (* 1.5 j2 pinvsq no-unkozai)
            temp2 (* 0.5 temp1 j2 pinvsq)
            temp3 (* -0.46875 j4 pinvsq pinvsq no-unkozai)
            mdot (+ no-unkozai (* 0.5 temp1 rteosq con41)
                    (* 0.0625 temp2 rteosq (+ (- 13.0 (* 78.0 cosio2)) (* 137.0 cosio4))))
            argpdot (+ (* -0.5 temp1 con42)
                       (* 0.0625 temp2 (+ (- 7.0 (* 114.0 cosio2)) (* 395.0 cosio4)))
                       (* temp3 (+ (- 3.0 (* 36.0 cosio2)) (* 49.0 cosio4))))
            xhdot1 (* (- temp1) cosio)
            nodedot (+ xhdot1 (* (+ (* 0.5 temp2 (- 4.0 (* 19.0 cosio2)))
                                    (* 2.0 temp3 (- 3.0 (* 7.0 cosio2))))
                                 cosio))
            xpidot (+ argpdot nodedot)
            delmotemp (+ 1.0 (* eta (math/cos mo)))
            satrec (merge satrec
                          {:isimp isimp :eta eta :cc1 cc1 :cc4 cc4 :cc5 cc5 :x1mth2 x1mth2
                           :mdot mdot :argpdot argpdot :nodedot nodedot
                           :omgcof (* bstar cc3 (math/cos argpo))
                           :xmcof (if (> ecco 1.0e-4) (/ (* (- x2o3) coef bstar) eeta) 0.0)
                           :nodecf (* 3.5 omeosq xhdot1 cc1)
                           :t2cof (* 1.5 cc1)
                           :xlcof (/ (* -0.25 j3oj2 sinio (+ 3.0 (* 5.0 cosio)))
                                     (if (> (abs (+ cosio 1.0)) 1.5e-12) (+ 1.0 cosio) temp4))
                           :aycof (* -0.5 j3oj2 sinio)
                           :delmo (* delmotemp delmotemp delmotemp)
                           :sinmao (math/sin mo)
                           :x7thm1 (- (* 7.0 cosio2) 1.0)})
            ;; deep space, for periods of 225 minutes and more. The C++
            ;; also calls dpper here, with init 'y', which changes nothing.
            satrec (if (>= (/ (* 2.0 pi) no-unkozai) 225.0)
                     (let [tc 0.0
                           ds (dscom epoch ecco argpo tc inclo nodeo no-unkozai)]
                       (merge satrec
                              (select-keys ds [:e3 :ee2 :peo :pgho :pho :pinco :plo :se2 :se3
                                               :sgh2 :sgh3 :sgh4 :sh2 :sh3 :si2 :si3 :sl2 :sl3
                                               :sl4 :xgh2 :xgh3 :xgh4 :xh2 :xh3 :xi2 :xi3 :xl2
                                               :xl3 :xl4 :zmol :zmos])
                              {:method \d :isimp 1}
                              (dsinit xke ds argpo tc gsto mo mdot no-unkozai nodeo nodedot
                                      xpidot ecco eccsq (:em ds) (:nm ds) inclo)))
                     satrec)]
        (if (not= (:isimp satrec) 1)
          (let [cc1sq (* cc1 cc1)
                d2 (* 4.0 ao tsi cc1sq)
                temp (/ (* d2 tsi cc1) 3.0)
                d3 (* (+ (* 17.0 ao) sfour) temp)
                d4 (* 0.5 temp ao tsi (+ (* 221.0 ao) (* 31.0 sfour)) cc1)]
            (assoc satrec :d2 d2 :d3 d3 :d4 d4
                   :t3cof (+ d2 (* 2.0 cc1sq))
                   :t4cof (* 0.25 (+ (* 3.0 d3) (* cc1 (+ (* 12.0 d2) (* 10.0 cc1sq)))))
                   :t5cof (* 0.2 (+ (* 3.0 d4) (* 12.0 cc1 d3) (* 6.0 d2 d2)
                                    (* 15.0 cc1sq (+ (* 2.0 d2) cc1sq))))))
          satrec)))))

;; -------------------------------------------------------------- sgp4

(defn sgp4
  "Propagate `satrec` to `tsince` minutes from its epoch: `{:r :v :error 0}`,
  TEME km and km/s, with the singly averaged mean elements `:am :em :im
  :Om :om :mm :nm` the C++ also records. Where the theory breaks down,
  `:error` is the C++ code and there is no state: 1 mean eccentricity out
  of range, 2 mean motion not positive, 3 perturbed eccentricity out of
  range, 4 semilatus rectum negative. Code 6, orbit decayed, still carries
  the state computed, as the C++ does."
  [{:keys [mo mdot argpo argpdot nodeo nodedot nodecf cc1 bstar cc4 t2cof isimp
           omgcof eta xmcof delmo d2 d3 d4 cc5 sinmao t3cof t4cof t5cof no-unkozai
           ecco inclo method radiusearthkm xke j2 j3oj2 opsmode]
    :as satrec}
   tsince]
  (let [vkmpersec (/ (* radiusearthkm xke) 60.0)
        t tsince
        ;; update for secular gravity and atmospheric drag
        xmdf (+ mo (* mdot t))
        argpdf (+ argpo (* argpdot t))
        nodedf (+ nodeo (* nodedot t))
        t2 (* t t)
        nodem (+ nodedf (* nodecf t2))
        tempa (- 1.0 (* cc1 t))
        tempe (* bstar cc4 t)
        templ (* t2cof t2)
        [mm argpm tempa tempe templ]
        (if (not= isimp 1)
          (let [delomg (* omgcof t)
                delmtemp (+ 1.0 (* eta (math/cos xmdf)))
                delm (* xmcof (- (* delmtemp delmtemp delmtemp) delmo))
                temp (+ delomg delm)
                mm (+ xmdf temp)
                t3 (* t2 t)
                t4 (* t3 t)]
            [mm
             (- argpdf temp)
             (- tempa (* d2 t2) (* d3 t3) (* d4 t4))
             (+ tempe (* bstar cc5 (- (math/sin mm) sinmao)))
             (+ templ (* t3cof t3) (* t4 (+ t4cof (* t t5cof))))])
          [xmdf argpdf tempa tempe templ])
        [em argpm inclm mm nodem nm]
        (if (= method \d)
          (dspace satrec t t ecco argpm inclo mm nodem no-unkozai)
          [ecco argpm inclo mm nodem no-unkozai])]
    (if (<= nm 0.0)
      {:error 2}
      (let [am (* (math/pow (/ xke nm) x2o3) tempa tempa)
            nm (/ xke (math/pow am 1.5))
            em (- em tempe)]
        (if (or (>= em 1.0) (< em -0.001))
          {:error 1}
          (let [em (if (< em 1.0e-6) 1.0e-6 em)
                mm (+ mm (* no-unkozai templ))
                xlm (+ mm argpm nodem)
                nodem (fmod nodem twopi)
                argpm (fmod argpm twopi)
                xlm (fmod xlm twopi)
                mm (fmod (- xlm argpm nodem) twopi)
                mean {:am am :em em :im inclm :Om nodem :om argpm :mm mm :nm nm}
                sinim (math/sin inclm)
                cosim (math/cos inclm)
                ;; add lunar-solar periodics
                [ep xincp nodep argpp mp]
                (if (= method \d)
                  (let [[ep xincp nodep argpp mp]
                        (dpper satrec t \n em inclm nodem argpm mm opsmode)]
                    (if (< xincp 0.0)
                      [ep (- xincp) (+ nodep pi) (- argpp pi) mp]
                      [ep xincp nodep argpp mp]))
                  [em inclm nodem argpm mm])]
            (if (and (= method \d) (or (< ep 0.0) (> ep 1.0)))
              (assoc mean :error 3)
              (let [[sinip cosip aycof xlcof]
                    (if (= method \d)
                      (let [sinip (math/sin xincp)
                            cosip (math/cos xincp)]
                        [sinip cosip (* -0.5 j3oj2 sinip)
                         (/ (* -0.25 j3oj2 sinip (+ 3.0 (* 5.0 cosip)))
                            (if (> (abs (+ cosip 1.0)) 1.5e-12) (+ 1.0 cosip) temp4))])
                      [sinim cosim (:aycof satrec) (:xlcof satrec)])
                    ;; long period periodics
                    axnl (* ep (math/cos argpp))
                    temp (/ 1.0 (* am (- 1.0 (* ep ep))))
                    aynl (+ (* ep (math/sin argpp)) (* temp aycof))
                    xl (+ mp argpp nodep (* temp xlcof axnl))
                    ;; solve Kepler's equation. As in the C++, sineo1 and
                    ;; coseo1 are those of eo1 before its last correction.
                    u (fmod (- xl nodep) twopi)
                    [sineo1 coseo1]
                    (loop [eo1 u tem5 9999.9 ktr 1 sc nil]
                      (if (and (>= (abs tem5) 1.0e-12) (<= ktr 10))
                        (let [sineo1 (math/sin eo1)
                              coseo1 (math/cos eo1)
                              tem5 (- 1.0 (* coseo1 axnl) (* sineo1 aynl))
                              tem5 (/ (- (+ (- u (* aynl coseo1)) (* axnl sineo1)) eo1) tem5)
                              tem5 (if (>= (abs tem5) 0.95) (if (> tem5 0.0) 0.95 -0.95) tem5)]
                          (recur (+ eo1 tem5) tem5 (inc ktr) [sineo1 coseo1]))
                        sc))
                    ;; short period preliminary quantities
                    ecose (+ (* axnl coseo1) (* aynl sineo1))
                    esine (- (* axnl sineo1) (* aynl coseo1))
                    el2 (+ (* axnl axnl) (* aynl aynl))
                    pl (* am (- 1.0 el2))]
                (if (< pl 0.0)
                  (assoc mean :error 4)
                  (let [rl (* am (- 1.0 ecose))
                        rdotl (/ (* (math/sqrt am) esine) rl)
                        rvdotl (/ (math/sqrt pl) rl)
                        betal (math/sqrt (- 1.0 el2))
                        temp (/ esine (+ 1.0 betal))
                        sinu (* (/ am rl) (- sineo1 aynl (* axnl temp)))
                        cosu (* (/ am rl) (+ (- coseo1 axnl) (* aynl temp)))
                        su (math/atan2 sinu cosu)
                        sin2u (* (+ cosu cosu) sinu)
                        cos2u (- 1.0 (* 2.0 sinu sinu))
                        temp (/ 1.0 pl)
                        temp1 (* 0.5 j2 temp)
                        temp2 (* temp1 temp)
                        [con41 x1mth2 x7thm1]
                        (if (= method \d)
                          (let [cosisq (* cosip cosip)]
                            [(- (* 3.0 cosisq) 1.0) (- 1.0 cosisq) (- (* 7.0 cosisq) 1.0)])
                          [(:con41 satrec) (:x1mth2 satrec) (:x7thm1 satrec)])
                        ;; update for short period periodics
                        mrt (+ (* rl (- 1.0 (* 1.5 temp2 betal con41)))
                               (* 0.5 temp1 x1mth2 cos2u))
                        su (- su (* 0.25 temp2 x7thm1 sin2u))
                        xnode (+ nodep (* 1.5 temp2 cosip sin2u))
                        xinc (+ xincp (* 1.5 temp2 cosip sinip cos2u))
                        mvt (- rdotl (/ (* nm temp1 x1mth2 sin2u) xke))
                        rvdot (+ rvdotl (/ (* nm temp1 (+ (* x1mth2 cos2u) (* 1.5 con41))) xke))
                        ;; orientation vectors
                        sinsu (math/sin su)
                        cossu (math/cos su)
                        snod (math/sin xnode)
                        cnod (math/cos xnode)
                        sini (math/sin xinc)
                        cosi (math/cos xinc)
                        xmx (* (- snod) cosi)
                        xmy (* cnod cosi)
                        ux (+ (* xmx sinsu) (* cnod cossu))
                        uy (+ (* xmy sinsu) (* snod cossu))
                        uz (* sini sinsu)
                        vx (- (* xmx cossu) (* cnod sinsu))
                        vy (- (* xmy cossu) (* snod sinsu))
                        vz (* sini cossu)]
                    (assoc mean
                           :r [(* (* mrt ux) radiusearthkm)
                               (* (* mrt uy) radiusearthkm)
                               (* (* mrt uz) radiusearthkm)]
                           :v [(* (+ (* mvt ux) (* rvdot vx)) vkmpersec)
                               (* (+ (* mvt uy) (* rvdot vy)) vkmpersec)
                               (* (+ (* mvt uz) (* rvdot vz)) vkmpersec)]
                           ;; sgp4fix: the satellite has decayed
                           :error (if (< mrt 1.0) 6 0))))))))))))

;; ------------------------------------------------------ TLE and epoch

(defn days2mdhms
  "[mon day hr minute sec] of day-of-year `days` in `year`, by the C++'s
  arithmetic, in which a leap year is one divisible by four."
  [year days]
  (let [dayofyr (long (math/floor days))
        lmonth [0 31 (if (zero? (mod year 4)) 29 28) 31 30 31 30 31 31 30 31 30 31]
        [i inttemp] (loop [i 1 inttemp 0]
                      (if (and (> dayofyr (+ inttemp (lmonth i))) (< i 12))
                        (recur (inc i) (+ inttemp (lmonth i)))
                        [i inttemp]))
        temp (* (- days dayofyr) 24.0)
        hr (long (math/floor temp))
        temp (* (- temp hr) 60.0)
        minute (long (math/floor temp))
        sec (* (- temp minute) 60.0)]
    [i (- dayofyr inttemp) hr minute sec]))

(defn jday
  "[jd jdFrac], a Julian Date split into its day and its fraction, as the
  C++ computes it."
  [year mon day hr minute sec]
  (let [jd (+ (- (* 367.0 year)
                 (math/floor (* (* 7 (+ year (math/floor (/ (+ mon 9) 12.0)))) 0.25)))
              (math/floor (/ (* 275 mon) 9.0))
              day 1721013.5)
        jdfrac (/ (+ sec (* minute 60.0) (* hr 3600.0)) 86400.0)]
    (if (> (abs jdfrac) 1.0)
      (let [dtt (math/floor jdfrac)] [(+ jd dtt) (- jdfrac dtt)])
      [jd jdfrac])))

(defn- implied-decimal
  "The C++'s reading of a TLE field like '-12345-3': its sign and a
  decimal point before the digits, times ten to the exponent."
  [sign digits exponent]
  (* (parse-double (str (if (= (str sign) "-") "-" "") "." (str/replace digits " " "0")))
     (math/pow 10.0 (parse-long (str/replace (str/trim exponent) "+" "")))))

(defn twoline2rv
  "The initialized element record of a two-line element set, read by
  column as the C++ `twoline2rv` reads it. Options as for `sgp4init`. Also
  carries `:satnum`, `:epochyr`, `:epochdays`, `:jdsatepoch` and
  `:jdsatepochF`."
  ([longstr1 longstr2] (twoline2rv longstr1 longstr2 {}))
  ([longstr1 longstr2 opts]
   (let [deg2rad (/ pi 180.0)
         xpdotp (/ 1440.0 (* 2.0 pi))
         col (fn [s a b] (str/trim (subs s a b)))
         epochyr (parse-long (col longstr1 18 20))
         epochdays (parse-double (col longstr1 20 32))
         ;; temporary fix for years from 1957-2056, as in the C++
         year (if (< epochyr 57) (+ epochyr 2000) (+ epochyr 1900))
         [mon day hr minute sec] (days2mdhms year epochdays)
         [jdsatepoch jdsatepochF] (jday year mon day hr minute sec)]
     (sgp4init
      (merge {:satnum (col longstr1 2 7)
              :epochyr epochyr
              :epochdays epochdays
              :jdsatepoch jdsatepoch
              :jdsatepochF jdsatepochF
              :epoch (- (+ jdsatepoch jdsatepochF) 2433281.5)
              :ndot (/ (parse-double (col longstr1 33 43)) (* xpdotp 1440.0))
              :nddot (/ (implied-decimal (nth longstr1 44) (subs longstr1 45 50) (subs longstr1 50 52))
                        (* xpdotp 1440.0 1440))
              :bstar (implied-decimal (nth longstr1 53) (subs longstr1 54 59) (subs longstr1 59 61))
              :inclo (* (parse-double (col longstr2 8 16)) deg2rad)
              :nodeo (* (parse-double (col longstr2 17 25)) deg2rad)
              :ecco (parse-double (str "." (str/replace (subs longstr2 26 33) " " "0")))
              :argpo (* (parse-double (col longstr2 34 42)) deg2rad)
              :mo (* (parse-double (col longstr2 43 51)) deg2rad)
              :no-kozai (/ (parse-double (col longstr2 52 63)) xpdotp)}
             opts)))))
