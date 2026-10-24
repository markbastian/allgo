(ns allgo.physics.convex
  "Contact as one convex problem, solved by Newton's method: the shape of
  MuJoCo's solver (Todorov, \"Convex and analytically-invertible dynamics
  with contacts and constraints\", ICRA 2014) and of Drake's SAP (Castro,
  Han and Permenter, \"An unconstrained convex formulation of compliant
  contact\", IEEE T-RO 2022), whose formulation this follows.

  ## Why

  `allgo.physics.solver`'s first three solvers visit contacts one at a
  time, each correcting itself against its neighbors' current answers.
  What the floor knows reaches the top of a column one contact per sweep,
  so a column taller than its iterations can carry never learns its own
  weight, and leans. Here every constraint is solved at once: each Newton
  step is one linear system coupling every body to every other through
  what joins them, and it converges in a handful of steps however tall
  the stack.

  ## The problem

  The velocities `v` at the end of the step minimize

      l(v) = 1/2 (v - v*)' M (v - v*)  +  sum over constraints of l_k

  where `v*` is where everything would be going with nothing joining it
  and `M` the masses. A constraint's own term is `1/2 g' R g`, with its
  impulse `g` a projection of `y = -R^-1 (vc - vhat)`: `vc` its rows'
  velocity, `vhat` the velocity it asks for, and `R` a small
  regularization that makes it very stiff rather than rigid. Three
  kinds, by what the impulse may be:

    contact      three rows -- normal, then two tangents -- and the
                 impulse inside the friction cone, projected in `R`'s
                 norm
    unilateral   one row, the impulse a push and never a pull: a joint
                 limit
    bilateral    one row, the impulse either way: a pin

  The gradient is `M (v - v*) - J' g` and the Hessian `M + J' G J`, `G`
  the derivative of the impulse, so a Newton step is one solve with that
  matrix, and a backtracking line search keeps every step a descent.
  `l` is convex whatever the friction, which is what the cone in `R`'s
  norm buys over the usual complementarity problem: there is one answer
  and Newton finds it.

  ## Blocks

  What moves is a set of blocks, each with its own dense mass matrix: six
  numbers for a loose body -- linear velocity and angular, in the world
  -- or a whole articulated model's generalized velocity, which is what
  MuJoCo's form is built for. A constraint joins at most two blocks, and
  each of its rows is a stretch of numbers against each.

  ## What this namespace knows

  Only the arithmetic, on flat arrays, as `solve` describes. Where the
  rows come from -- collision detection, joints, warm starting, moving
  things afterward -- is `allgo.physics.solver`'s and
  `allgo.physics.world`'s, and this is a `:solver` in each."
  (:require [allgo.array :as a]
            [clojure.math :as math]))

(def contact 0)
(def unilateral 1)
(def bilateral 2)

;; ---------------------------------------------------------------------------
;; One contact: its cone, its impulse and the impulse's derivative

(defn- project!
  "The impulse of one contact and its derivative, from its `y`.

  Projects `y` (normal, tangent, tangent) onto the friction cone in the
  norm of `R` (`rn` along the normal, `rt` along both tangents) into
  `g`, and writes `G` = -dg/dvc, a symmetric 3x3, into `gm` row-major.
  Returns the contact's cost, `1/2 g' R g`.

  In the stretched coordinates `z = R^1/2 g` the cone is a round one of
  slope `mu sqrt(rt/rn)` and the projection is the Euclidean one: inside
  it is the identity, in its polar cone zero, and between the two it
  lands on the cone's surface. `G` is that projection's derivative
  pulled back, `R^-1/2 P' R^-1/2`."
  [mu rn rt ^doubles y ^doubles g ^doubles gm]
  (let [mu (double mu) rn (double rn) rt (double rt)
        sn (math/sqrt rn) st (math/sqrt rt)
        yn (* sn (aget y 0)) y1 (* st (aget y 1)) y2 (* st (aget y 2))
        m (* mu (/ st sn))
        r (math/sqrt (+ (* y1 y1) (* y2 y2)))]
    (cond
      ;; sticking: the impulse is whatever the contact asks for
      (<= r (* m yn))
      (do (aset g 0 (/ yn sn)) (aset g 1 (/ y1 st)) (aset g 2 (/ y2 st))
          (dotimes [i 9] (aset gm i 0.0))
          (aset gm 0 (/ 1.0 rn)) (aset gm 4 (/ 1.0 rt)) (aset gm 8 (/ 1.0 rt))
          (* 0.5 (+ (* yn yn) (* y1 y1) (* y2 y2))))
      ;; separating: nothing
      (<= (* m r) (- yn))
      (do (aset g 0 0.0) (aset g 1 0.0) (aset g 2 0.0)
          (dotimes [i 9] (aset gm i 0.0))
          0.0)
      ;; sliding: on the cone's surface, the friction opposing the slip
      :else
      (let [k (/ 1.0 (+ 1.0 (* m m)))
            zn (* k (+ yn (* m r)))
            u1 (/ y1 r) u2 (/ y2 r)
            z1 (* m zn u1) z2 (* m zn u2)
            ;; P' in the stretched coordinates, normal first
            pnn k
            pn1 (* k m u1) pn2 (* k m u2)
            c (/ (* m zn) r)
            p11 (+ (* k m m u1 u1) (* c (- 1.0 (* u1 u1))))
            p22 (+ (* k m m u2 u2) (* c (- 1.0 (* u2 u2))))
            p12 (- (* k m m u1 u2) (* c u1 u2))
            ;; pulled back by R^-1/2 on both sides
            in (/ 1.0 sn) it (/ 1.0 st)]
        (aset g 0 (* zn in)) (aset g 1 (* z1 it)) (aset g 2 (* z2 it))
        (aset gm 0 (* pnn in in)) (aset gm 1 (* pn1 in it)) (aset gm 2 (* pn2 in it))
        (aset gm 3 (* pn1 in it)) (aset gm 4 (* p11 it it)) (aset gm 5 (* p12 it it))
        (aset gm 6 (* pn2 in it)) (aset gm 7 (* p12 it it)) (aset gm 8 (* p22 it it))
        (* 0.5 (+ (* zn zn) (* z1 z1) (* z2 z2)))))))

(defn- one-row!
  "The impulse of a one-row constraint and its derivative: a push, never
  a pull, for a unilateral one; either way for a bilateral one. Into
  `g` and `gm` as `project!` writes them, and its cost returned."
  [unilateral? rn ^doubles y ^doubles g ^doubles gm]
  (let [rn (double rn)
        yn (aget y 0)
        on? (or (not unilateral?) (pos? yn))]
    (dotimes [i 9] (aset gm i 0.0))
    (aset g 1 0.0) (aset g 2 0.0)
    (if on?
      (do (aset g 0 yn) (aset gm 0 (/ 1.0 rn)) (* 0.5 rn yn yn))
      (do (aset g 0 0.0) 0.0))))

;; ---------------------------------------------------------------------------
;; The problem's arithmetic
;;
;; Blocks: `:nb` of them, block `b` its `:bsize` numbers from `:boff` in
;; the velocity vector and its dense mass matrix from `:moff` in `:mass`
;; (and inverse in `:minv`), row-major. Constraints: `:nk` of them, each
;; its `:kind`, its rows from `:krow` (three for a contact, one otherwise),
;; joining blocks `:ka` and `:kb` (-1 for nothing that moves), and each row
;; `r` its numbers against them from `:va` and `:vb` in `:vals`.

(defn- rows-of ^long [^ints kind k]
  (if (zero? (aget kind (long k))) 3 1))

(defn- row-velocities!
  "Every row's velocity at `v`, into `out`: one pass over the
  constraints, the arrays read once rather than per row."
  [{:keys [^long nk ^ints kind ^ints krow ^ints ka ^ints kb ^ints boff ^ints bsize
           ^ints va ^ints vb ^doubles vals]} ^doubles v ^doubles out]
  (dotimes [k nk]
    (let [a (aget ka k) b (aget kb k) r0 (aget krow k)
          nr (if (zero? (aget kind k)) 3 1)
          oa (long (if (neg? a) 0 (aget boff a))) na (long (if (neg? a) 0 (aget bsize a)))
          ob (long (if (neg? b) 0 (aget boff b))) nb (long (if (neg? b) 0 (aget bsize b)))]
      (dotimes [i nr]
        (let [r (+ r0 i) sa (aget va r) sb (aget vb r)]
          (aset out r
                (+ (double (loop [j 0 acc 0.0]
                             (if (== j na) acc (recur (inc j) (+ acc (* (aget vals (+ sa j)) (aget v (+ oa j))))))))
                   (double (loop [j 0 acc 0.0]
                             (if (== j nb) acc (recur (inc j) (+ acc (* (aget vals (+ sb j)) (aget v (+ ob j))))))))))))))
  out)

(defn- block-mul!
  "`out` = (`mats`, `M` or `M^-1`, block by block) times `x`."
  [{:keys [^long nb ^ints boff ^ints bsize ^ints moff]} ^doubles mats ^doubles x ^doubles out]
  (dotimes [b nb]
    (let [o (aget boff b) n (aget bsize b) m (aget moff b)]
      (dotimes [r n]
        (aset out (+ o r)
              (double (loop [c 0 acc 0.0]
                        (if (== c n) acc (recur (inc c) (+ acc (* (aget mats (+ m (* r n) c)) (aget x (+ o c))))))))))))
  out)

(defn- impulses!
  "Every constraint's impulse at row velocities `vc`, into `gam` (three
  a constraint) and, when `gms` is given, its `G` into `gms` (nine).
  Returns the constraints' share of the cost."
  ^double [{:keys [^long nk ^ints kind ^ints krow ^doubles rn ^doubles rt ^doubles vhat ^doubles mu]}
           ^doubles vc ^doubles gam gms]
  (let [y (a/f64 3) g (a/f64 3) gm (a/f64 9)]
    (loop [k 0 total 0.0]
      (if (== k nk)
        total
        (let [kd (aget kind k) r0 (aget krow k)
              cost (if (zero? kd)
                     (let [vn (aget vc r0) v1 (aget vc (+ r0 1)) v2 (aget vc (+ r0 2))]
                       (aset y 0 (/ (- (aget vhat k) vn) (aget rn k)))
                       (aset y 1 (/ (- v1) (aget rt k)))
                       (aset y 2 (/ (- v2) (aget rt k)))
                       (project! (aget mu k) (aget rn k) (aget rt k) y g gm))
                     (let [vn (aget vc r0)]
                       (aset y 0 (/ (- (aget vhat k) vn) (aget rn k)))
                       (one-row! (== kd unilateral) (aget rn k) y g gm)))
              k3 (* 3 k)]
          (aset gam k3 (aget g 0)) (aset gam (+ k3 1) (aget g 1)) (aset gam (+ k3 2) (aget g 2))
          (when gms
            (let [^doubles gms gms k9 (* 9 k)]
              (dotimes [i 9] (aset gms (+ k9 i) (aget gm i)))))
          (recur (inc k) (+ total (double cost))))))))

(defn- evaluate!
  "The cost at `v`, each constraint's impulse into `gam` (three a
  constraint) and, when `gms` is given, its `G` into `gms` (nine)."
  [{:keys [^long nd ^doubles mass ^doubles vstar ^long nr] :as problem}
   ^doubles v ^doubles gam gms ^doubles scratch]
  (let [dv (a/f64 nd)
        ^doubles vc (row-velocities! problem v (a/f64 nr))]
    (dotimes [i nd] (aset dv i (- (aget v i) (aget vstar i))))
    (block-mul! problem mass dv scratch)
    (let [kinetic (loop [i 0 s 0.0] (if (== i nd) s (recur (inc i) (+ s (* (aget dv i) (aget scratch i))))))]
      (+ (* 0.5 kinetic) (impulses! problem vc gam gms)))))

(defn- add-jt-gamma!
  "`out` += J' `gam`."
  [{:keys [^long nk ^ints kind ^ints krow ^ints ka ^ints kb ^ints boff ^ints bsize
           ^ints va ^ints vb ^doubles vals]} ^doubles gam ^doubles out]
  (dotimes [k nk]
    (let [a (aget ka k) b (aget kb k) r0 (aget krow k)
          nr (if (zero? (aget kind k)) 3 1)]
      (dotimes [i nr]
        (let [x (aget gam (+ (* 3 k) i)) r (+ r0 i)]
          (when-not (zero? x)
            (when-not (neg? a)
              (let [o (aget boff a) s (aget va r)]
                (dotimes [j (aget bsize a)] (aset out (+ o j) (+ (aget out (+ o j)) (* x (aget vals (+ s j))))))))
            (when-not (neg? b)
              (let [o (aget boff b) s (aget vb r)]
                (dotimes [j (aget bsize b)] (aset out (+ o j) (+ (aget out (+ o j)) (* x (aget vals (+ s j))))))))))))))

;; ---------------------------------------------------------------------------
;; The Newton matrix, stored by its envelope
;;
;; Two blocks are coupled only if something joins them, so the matrix is
;; almost all zeros, and factoring it densely -- 1,536 rows for a 16x16
;; wall -- took 600 ms a step. Numbered so that blocks that touch are
;; near each other (reverse Cuthill-McKee), the nonzeros and everything
;; the factor fills in lie close to the diagonal: each row is kept from
;; its first nonzero to the diagonal and no further, and a column of
;; bricks is a band six or twelve wide however tall it is.

(defn- ordering
  "Reverse Cuthill-McKee over the blocks, touching where a constraint
  joins two: each block's new position."
  [nb nk ^ints ka ^ints kb]
  (let [nb (long nb)
        adj (let [adj (object-array nb)]
              (dotimes [b nb] (aset adj b (transient #{})))
              (dotimes [k (long nk)]
                (let [a (aget ka k) b (aget kb k)]
                  (when (and (>= a 0) (>= b 0) (not (== a b)))
                    (aset adj a (conj! (aget adj a) b))
                    (aset adj b (conj! (aget adj b) a)))))
              (mapv #(persistent! %) adj))
        degree (fn [b] (count (nth adj b)))
        ^ints visited (a/i32 nb)
        order (loop [order (transient []) starts (sort-by degree (range nb))]
                (if-let [s (first starts)]
                  (if (pos? (aget visited (long s)))
                    (recur order (rest starts))
                    (do (aset visited (long s) 1)
                        (recur (loop [order order
                                      queue (conj #?(:clj clojure.lang.PersistentQueue/EMPTY
                                                     :cljs cljs.core/PersistentQueue.EMPTY)
                                                  s)]
                                 (if-let [b (peek queue)]
                                   (let [fresh (sort-by degree (remove #(pos? (aget visited (long %))) (nth adj b)))]
                                     (doseq [n fresh] (aset visited (long n) 1))
                                     (recur (conj! order b) (into (pop queue) fresh)))
                                   order))
                               (rest starts))))
                  (persistent! order)))
        ^ints perm (a/i32 nb)]
    (doseq [[k b] (map-indexed vector (rseq order))]
      (aset perm (long b) (int k)))
    {:perm perm :adj adj}))

(defn- envelope
  "The envelope for blocks placed by `perm` and joined as `adj`: each
  block's first row in the new numbering, and for every row where it
  starts -- the first column of the lowest-numbered block it touches, or
  its own -- and where it is stored."
  [{:keys [^long nb ^long nd ^ints bsize]} ^ints perm adj]
  (let [^ints by-place (a/i32 nb)
        _ (dotimes [b nb] (aset by-place (aget perm b) (int b)))
        ^ints noff (a/i32 nb)
        _ (loop [p 0 off 0]
            (when (< p nb)
              (let [b (aget by-place p)]
                (aset noff b (int off))
                (recur (inc p) (+ off (aget bsize b))))))
        ^ints first-col (a/i32 nd)
        ^ints start (a/i32 (inc nd))]
    (dotimes [b nb]
      (let [lowest (reduce (fn [m q] (if (< (aget perm (long q)) (aget perm (long m))) q m)) b (nth adj b))
            col (aget noff (long lowest))]
        (dotimes [r (aget bsize b)] (aset first-col (+ (aget noff b) r) col))))
    (dotimes [i nd]
      (aset start (inc i) (int (+ (aget start i) (- (inc i) (aget first-col i))))))
    {:first first-col :start start :size (aget start nd) :noff noff}))

(defn- at
  "Where row `i`, column `j` (`j` <= `i`, inside the envelope) is stored."
  ^long [^ints first-col ^ints start ^long i ^long j]
  (+ (aget start i) (- j (aget first-col i))))

(defn- factor!
  "Cholesky, in place, of the envelope-stored matrix in `h`: the lower
  factor overwrites it. Nothing outside the envelope fills in. False if
  the matrix is not positive definite."
  [^doubles h {:keys [^ints first ^ints start]} nd]
  ;; Each row is a run of the array, so `row i` + k is its column k: the
  ;; inner product of two rows is a plain walk along both runs.
  (let [nd (long nd)
        row (fn ^long [^long i] (- (aget start i) (aget first i)))]
    (loop [i 0]
      (if (== i nd)
        true
        (let [fi (aget first i)
              ri (long (row i))
              ok? (loop [j fi]
                    (if (> j i)
                      true
                      (let [rj (long (row j))
                            k0 (max fi (aget first j))
                            s (loop [k k0 s (aget h (+ ri j))]
                                (if (== k j) s (recur (inc k) (- s (* (aget h (+ ri k)) (aget h (+ rj k)))))))]
                        (if (< j i)
                          (do (aset h (+ ri j) (/ s (aget h (+ rj j))))
                              (recur (inc j)))
                          (if (pos? s)
                            (do (aset h (+ ri i) (math/sqrt s)) true)
                            false)))))]
          (if ok? (recur (inc i)) false))))))

(defn- substitute!
  "Solves `L L' x = b` with `factor!`'s `L`, overwriting `b`."
  [^doubles h {:keys [^ints first ^ints start]} nd ^doubles b]
  ;; as in `factor!`, row `i` + k is its column k
  (let [nd (long nd)]
    (dotimes [i nd]
      (let [ri (- (aget start i) (aget first i))
            s (loop [k (aget first i) s (aget b i)]
                (if (== k i) s (recur (inc k) (- s (* (aget h (+ ri k)) (aget b k))))))]
        (aset b i (/ s (aget h (+ ri i))))))
    (loop [i (dec nd)]
      (when (>= i 0)
        (let [ri (- (aget start i) (aget first i))
              x (/ (aget b i) (aget h (+ ri i)))]
          (aset b i x)
          (loop [k (aget first i)]
            (when (< k i)
              (aset b k (- (aget b k) (* (aget h (+ ri k)) x)))
              (recur (inc k))))
          (recur (dec i)))))
    b))

(defn- hessian!
  "`h` = M + J' G J, its lower triangle into the envelope, rows and
  columns renumbered as `env` places the blocks."
  [{:keys [^long nb ^long nk ^ints kind ^ints krow ^ints ka ^ints kb ^ints bsize ^ints moff
           ^ints va ^ints vb ^doubles vals ^doubles mass]}
   ^doubles gms {:keys [^ints first ^ints start size ^ints noff]} ^doubles h]
  (dotimes [i (long size)] (aset h i 0.0))
  (dotimes [b nb]
    (let [o (aget noff b) n (aget bsize b) m (aget moff b)]
      (dotimes [r n]
        (dotimes [c (inc r)]
          (let [idx (at first start (+ o r) (+ o c))]
            (aset h idx (+ (aget h idx) (aget mass (+ m (* r n) c)))))))))
  ;; Per constraint, for each block it joins, `w_i` = sum over j of
  ;; G_ij row_j against that block; then the block pair's entries are
  ;; sum over i of row_i against the first block times w_i.
  (let [stride (long (loop [b 0 m 1] (if (== b nb) m (recur (inc b) (max m (aget bsize b))))))
        w (a/f64 (* 3 stride))]
    (dotimes [k nk]
      (let [nr (rows-of kind k) r0 (aget krow k) k9 (* 9 k)
            ka-k (aget ka k) kb-k (aget kb k)]
        (dotimes [sq 2]
          (let [q (if (zero? sq) ka-k kb-k)]
            (when-not (neg? q)
              (let [nq (aget bsize q)
                    voff (if (zero? sq) va vb)]
                ;; w_i over block q
                (dotimes [i nr]
                  (dotimes [s nq]
                    (aset w (+ (* i stride) s)
                          (double (loop [j 0 acc 0.0]
                                    (if (== j nr) acc
                                        (recur (inc j) (+ acc (* (aget gms (+ k9 (* 3 i) j))
                                                                 (aget vals (+ (aget ^ints voff (+ r0 j)) s)))))))))))
                (dotimes [sp 2]
                  (let [p (if (zero? sp) ka-k kb-k)]
                    (when-not (neg? p)
                      (let [np (aget bsize p)
                            poff (if (zero? sp) va vb)
                            pn (aget noff p) qn (aget noff q)]
                        ;; the lower triangle only: block p at or below q
                        (when (>= pn qn)
                          (dotimes [r np]
                            ;; a row's columns from block q on are one run
                            ;; of the envelope, so its place is found once
                            (let [row (+ pn r)
                                  base (+ (- (aget start row) (aget first row)) qn)
                                  ns (if (> pn qn) nq (inc r))]
                              (dotimes [s ns]
                                (let [x (loop [i 0 acc 0.0]
                                          (if (== i nr) acc
                                              (recur (inc i) (+ acc (* (aget vals (+ (aget ^ints poff (+ r0 i)) r))
                                                                       (aget w (+ (* i stride) s)))))))
                                      idx (+ base s)]
                                  (aset h idx (+ (aget h idx) x))))))))))))))))))
  h)

;; Last step's factor, kept
;;
;; Factoring is most of a step on a big island: on a 16x16 wall, 4 ms a
;; factor of 1,536 rows, averaging 75 stored a row, and about two a step.
;; But the matrix is indexed by blocks, not constraints, and from one
;; step to the next it hardly changes: a sticking contact's G is the
;; constant R^-1 inside the cone, and the bodies move by microns. So the
;; last factor taken, kept across steps, is a preconditioner for the
;; next matrix, and conjugate gradients with it reaches 1e-8 in three to
;; eight iterations -- each a product with M + J' G J, never formed, and
;; a substitution -- on a factor eleven steps old. Any positive definite
;; factor gives the right direction; a stale one only takes more
;; iterations, and past `max-cg` of them the matrix is factored afresh.

(defn- hessian-mul!
  "`out` = (M + J' G J) `x`, without forming it: `jx` (a row each) and
  `t` (three a constraint) are scratch."
  [{:keys [^long nk ^ints kind ^ints krow ^doubles mass] :as problem}
   ^doubles gms ^doubles x ^doubles jx ^doubles t ^doubles out]
  (row-velocities! problem x jx)
  (dotimes [k nk]
    (let [nr (rows-of kind k) r0 (aget krow k) k3 (* 3 k) k9 (* 9 k)]
      (dotimes [i nr]
        (aset t (+ k3 i)
              (double (loop [j 0 s 0.0]
                        (if (== j nr)
                          s
                          (recur (inc j) (+ s (* (aget gms (+ k9 (* 3 i) j)) (aget jx (+ r0 j))))))))))))
  (block-mul! problem mass x out)
  (add-jt-gamma! problem t out)
  out)

(defn- precondition!
  "`out` = (L L')^-1 `x` with `factor`'s L, renumbered into its envelope's
  order and back through `y`."
  [{:keys [^long nb ^ints boff ^ints bsize]} {:keys [env ^doubles h]} ^doubles x ^doubles y ^doubles out]
  (let [^ints noff (:noff env) nd (alength x)]
    (dotimes [b nb]
      (let [o (aget boff b) p (aget noff b)]
        (dotimes [j (aget bsize b)] (aset y (+ p j) (aget x (+ o j))))))
    (substitute! h env nd y)
    (dotimes [b nb]
      (let [o (aget boff b) p (aget noff b)]
        (dotimes [j (aget bsize b)] (aset out (+ o j) (aget y (+ p j))))))
    out))

(defn- dot ^double [^doubles a ^doubles b]
  (loop [i 0 s 0.0] (if (== i (alength a)) s (recur (inc i) (+ s (* (aget a i) (aget b i)))))))

(defn- conjugate-gradients!
  "Solves (M + J' G J) `x` = `b` into `x`, preconditioned by `factor`, to
  a residual `tolerance` of `b`'s: the iterations it took, or nil if
  `max-cg` were not enough. `ws` is scratch, from `cg-scratch`."
  [problem gms factor ^doubles b ^doubles x tolerance max-cg
   {:keys [^doubles r ^doubles z ^doubles d ^doubles hd ^doubles y ^doubles jx ^doubles t]}]
  (let [nd (alength b)
        goal (* (double tolerance) (math/sqrt (dot b b)))]
    (dotimes [i nd] (aset x i 0.0) (aset r i (aget b i)))
    (precondition! problem factor r y z)
    (dotimes [i nd] (aset d i (aget z i)))
    (loop [it 0 rz (dot r z)]
      (cond
        (<= (math/sqrt (dot r r)) goal) it
        (>= it (long max-cg)) nil
        :else (let [_ (hessian-mul! problem gms d jx t hd)
                    a (/ rz (dot d hd))]
                (dotimes [i nd]
                  (aset x i (+ (aget x i) (* a (aget d i))))
                  (aset r i (- (aget r i) (* a (aget hd i)))))
                (precondition! problem factor r y z)
                (let [rz' (dot r z) beta (/ rz' rz)]
                  (dotimes [i nd] (aset d i (+ (aget z i) (* beta (aget d i)))))
                  (recur (inc it) rz')))))))

(defn- cg-scratch [{:keys [^long nd ^long nr ^long nk]}]
  {:r (a/f64 nd) :z (a/f64 nd) :d (a/f64 nd) :hd (a/f64 nd) :y (a/f64 nd)
   :jx (a/f64 nr) :t (a/f64 (* 3 nk))})

;; ---------------------------------------------------------------------------
;; Solving

(defn pushed
  "Where everything would be going if `impulses` -- three a constraint,
  as `solve` returns them -- were applied to `vstar`: `vstar + M^-1 J' g`.
  Last step's impulses through this are the warm start, the guess the
  Newton steps begin from."
  [{:keys [^long nd ^doubles vstar ^doubles minv] :as problem} ^doubles impulses]
  (let [f (a/f64 nd) v (a/f64 nd)]
    (add-jt-gamma! problem impulses f)
    (block-mul! problem minv f v)
    (dotimes [i nd] (aset v i (+ (aget v i) (aget vstar i))))
    v))

(defn cost
  "The step's cost at `v`."
  [{:keys [^long nd ^long nk] :as problem} ^doubles v]
  (evaluate! problem v (a/f64 (* 3 nk)) nil (a/f64 nd)))

(defn- line-search
  "How far to go along `dir` from `v`: where the cost stops falling.

  The cost is convex along any line, so its slope there is increasing
  and the step is that slope's root -- bracketed by 0, where it is
  `slope`, and by the full Newton step, and found by false position.
  SAP does the same. Halving a step until it is good enough, as this did
  at first, takes short steps whenever contacts change between sticking
  and sliding along the way, and a ragdoll on a pile of bricks needed
  thirty Newton steps a frame for it. Leaves `trial` at the point
  chosen.

  Every row's velocity is linear along the line, so each trial is the
  constraints' projections alone: the rows are multiplied out once, not
  once a trial."
  [{:keys [^long nd ^long nr ^long nk ^ints kind ^doubles mass ^doubles vstar] :as problem}
   ^doubles v ^doubles dir slope ^doubles trial ^doubles trial-gam ^doubles scratch]
  (let [^doubles md (block-mul! problem mass dir scratch)
        ^doubles jd (row-velocities! problem dir (a/f64 nr))
        ^doubles vc0 (row-velocities! problem v (a/f64 nr))
        vc (a/f64 nr)
        ^ints krow (:krow problem)
        ;; the slope at `v + a dir`: dir . M (v + a dir - v*) - (J dir) . g
        dv0 (loop [i 0 s 0.0] (if (== i nd) s (recur (inc i) (+ s (* (- (aget v i) (aget vstar i)) (aget md i))))))
        dd (loop [i 0 s 0.0] (if (== i nd) s (recur (inc i) (+ s (* (aget dir i) (aget md i))))))
        slope-at (fn [alpha]
                   (let [alpha (double alpha)]
                     (dotimes [i nr] (aset vc i (+ (aget vc0 i) (* alpha (aget jd i)))))
                     (impulses! problem vc trial-gam nil)
                     (- (+ dv0 (* alpha dd))
                        (loop [k 0 s 0.0]
                          (if (== k nk)
                            s
                            (let [r0 (aget krow k) k3 (* 3 k)]
                              (recur (inc k)
                                     (+ s (loop [i 0 s2 0.0]
                                            (if (== i (rows-of kind k))
                                              s2
                                              (recur (inc i) (+ s2 (* (aget jd (+ r0 i)) (aget trial-gam (+ k3 i)))))))))))))))
        slope (double slope)
        ;; past the full step while the cost is still falling there
        [lo slo hi shi] (loop [lo 0.0 slo slope hi 1.0]
                          (let [s (double (slope-at hi))]
                            (if (and (neg? s) (< hi 64.0))
                              (recur hi s (* 2.0 hi))
                              [lo slo hi s])))
        alpha (if (<= (double shi) 0.0)
                hi
                (loop [lo (double lo) slo (double slo) hi (double hi) shi (double shi) n 0]
                  (let [a (/ (- (* lo shi) (* hi slo)) (- shi slo))
                        a (if (and (> a lo) (< a hi)) a (* 0.5 (+ lo hi)))
                        s (double (slope-at a))]
                    (if (or (<= (abs s) (* 1e-6 (abs slope))) (>= n 40) (< (- hi lo) 1e-12))
                      a
                      (if (neg? s)
                        (recur a s hi (* 0.5 shi) (inc n))
                        (recur lo (* 0.5 slo) a s (inc n)))))))
        alpha (double alpha)]
    (dotimes [i nd] (aset trial i (+ (aget v i) (* alpha (aget dir i)))))
    alpha))

(defn- newton
  "`solve` on one island: Newton steps from `v0` (or `vstar`) until the
  gradient vanishes. Each direction is by conjugate gradients on
  `factor` (last step's, kept) while that converges in `max-cg`, else
  by a factor taken afresh, which is kept in its place and returned as
  `:factor`. One that took more than `refactor-after` is worn, and is
  replaced at the next direction: past about a dozen iterations a
  factor is the cheaper of the two. An island smaller than `cg-from`
  rows is factored afresh every direction, which there costs less than
  the iterations would.

  At most `max-iterations`, a hundred as MuJoCo has it. Thirty was the
  limit once, and twenty ragdolls landing in a pile take up to sixty:
  stopped short, a step's answer is lower in cost than where it began
  and nothing like the minimum -- joints turning at hundreds of radians
  a second -- and the pile flew apart in a second and a half."
  [{:keys [^long nb ^long nd ^long nk ^ints ka ^ints kb ^doubles mass ^doubles vstar] :as problem}
   v0 {:keys [max-iterations tolerance factor max-cg refactor-after cg-tolerance cg-from]
       :or {max-iterations 100 tolerance 1e-10 max-cg 15 refactor-after 8 cg-tolerance 1e-8 cg-from 100}}]
  (let [v (a/f64 nd)
        ^doubles v-init (or v0 vstar)
        _ (dotimes [i nd] (aset v i (aget v-init i)))
        gam (a/f64 (* 3 nk))
        gms (a/f64 (* 9 nk))
        ;; the ordering only when a factor is taken here
        env (delay (let [{:keys [perm adj]} (ordering nb nk ka kb)] (envelope problem perm adj)))
        fresh! (fn []
                 (let [env @env h (a/f64 (long (:size env)))]
                   (hessian! problem gms env h)
                   (when (factor! h env nd) {:env env :h h :nd nd})))
        grad (a/f64 nd)
        scratch (a/f64 nd)
        trial (a/f64 nd)
        trial-gam (a/f64 (* 3 nk))
        dir (a/f64 nd)
        ws (cg-scratch problem)]
    (loop [it 0 cost (evaluate! problem v gam gms scratch)
           fac (when (= nd (:nd factor)) factor)]
       ;; the gradient, M (v - v*) - J' g
      (dotimes [i nd] (aset trial i (- (aget v i) (aget vstar i))))
      (block-mul! problem mass trial grad)
      (dotimes [i nd] (aset grad i (- (aget grad i))))
      (add-jt-gamma! problem gam grad)
      (let [gnorm (loop [i 0 s 0.0] (if (== i nd) s (recur (inc i) (max s (abs (aget grad i))))))]
        ;; `grad` holds minus the gradient from here
        (if (or (>= it (long max-iterations)) (<= gnorm (double tolerance)) (zero? nd))
          {:v v :impulses gam :iterations it :cost cost :factor fac}
          (let [cg (when (and fac (not (:worn? fac)) (>= nd (long cg-from)))
                     (conjugate-gradients! problem gms fac grad dir cg-tolerance max-cg ws))
                fac (cond
                      (nil? cg) (when-let [f (fresh!)]
                                  (precondition! problem f grad (:y ws) dir)
                                  f)
                      (> (long cg) (long refactor-after)) (assoc fac :worn? true)
                      :else fac)]
            (if-not fac
              {:v v :impulses gam :iterations it :cost cost :failed :not-positive-definite}
              (let [slope (- (dot grad dir))
                    alpha (line-search problem v dir slope trial trial-gam scratch)
                    step (loop [i 0 s 0.0] (if (== i nd) s (recur (inc i) (max s (abs (* (double alpha) (aget dir i)))))))]
                (dotimes [i nd] (aset v i (aget trial i)))
                (let [cost' (evaluate! problem v gam gms scratch)]
                  (if (< step 1e-12)
                    {:v v :impulses gam :iterations (inc it) :cost cost' :factor fac}
                    (recur (inc it) cost' fac)))))))))))

(defn- islands
  "The blocks constraints join, directly or through others: each
  island's blocks and constraints, in their own order. A block nothing
  touches is in none."
  [nb nk ^ints ka ^ints kb]
  (let [nb (long nb)
        ^ints root (a/i32 (range nb))
        find (fn [b] (loop [b (long b)]
                       (let [p (aget root b)]
                         (if (== p b) b (do (aset root b (aget root p)) (recur p))))))]
    (dotimes [k (long nk)]
      (let [a (aget ka k) b (aget kb k)]
        (when (and (>= a 0) (>= b 0))
          (aset root (long (find a)) (int (find b))))))
    (let [key (fn [k] (find (if (neg? (aget ka k)) (aget kb k) (aget ka k))))
          ks (group-by key (range nk))]
      (vec (for [[r cs] ks]
             [(vec (filter #(== r (find %)) (range nb))) cs])))))

(defn- part
  "The island of `blocks` and constraints `ks` as a problem of its own.
  On arrays throughout: a wall of 255 bricks beside one alone is copied
  here every step."
  [{:keys [^long nb ^ints bsize ^ints boff ^ints moff ^doubles mass ^doubles minv ^doubles vstar
           ^ints kind ^ints krow ^ints ka ^ints kb ^ints va ^ints vb ^doubles vals
           ^doubles rn ^doubles rt ^doubles vhat ^doubles mu]} blocks ks]
  (let [^ints bl (a/i32 blocks) ^ints kl (a/i32 ks)
        nb' (alength bl) nk (alength kl)
        copy! (fn [^doubles from f ^doubles to t n]
                (let [f (long f) t (long t)] (dotimes [i (long n)] (aset to (+ t i) (aget from (+ f i))))))
        ;; old block -> new, -1 for none
        ^ints at (a/i32 nb -1)
        _ (dotimes [i nb'] (aset at (aget bl i) (int i)))
        ^ints bsize' (a/i32 nb') ^ints boff' (a/i32 nb') ^ints moff' (a/i32 nb')
        [nd nm] (loop [i 0 o 0 m 0]
                  (if (== i nb')
                    [o m]
                    (let [n (aget bsize (aget bl i))]
                      (aset bsize' i n) (aset boff' i (int o)) (aset moff' i (int m))
                      (recur (inc i) (+ o n) (+ m (* n n))))))
        nd (long nd)
        mass' (a/f64 (long nm)) minv' (a/f64 (long nm)) vstar' (a/f64 nd)
        _ (dotimes [i nb']
            (let [q (aget bl i) n (aget bsize q)]
              (copy! mass (aget moff q) mass' (aget moff' i) (* n n))
              (copy! minv (aget moff q) minv' (aget moff' i) (* n n))
              (copy! vstar (aget boff q) vstar' (aget boff' i) n)))
        size (fn ^long [^long q] (if (neg? q) 0 (aget bsize q)))
        ^ints kind' (a/i32 nk) ^ints krow' (a/i32 nk) ^ints ka' (a/i32 nk) ^ints kb' (a/i32 nk)
        rn' (a/f64 nk) rt' (a/f64 nk) vhat' (a/f64 nk) mu' (a/f64 nk)
        [nr nv] (loop [i 0 r 0 v 0]
                  (if (== i nk)
                    [r v]
                    (let [k (aget kl i) n (rows-of kind k) a (aget ka k) b (aget kb k)]
                      (aset kind' i (aget kind k)) (aset krow' i (int r))
                      (aset ka' i (int (if (neg? a) -1 (aget at a))))
                      (aset kb' i (int (if (neg? b) -1 (aget at b))))
                      (aset rn' i (aget rn k)) (aset rt' i (aget rt k))
                      (aset vhat' i (aget vhat k)) (aset mu' i (aget mu k))
                      (recur (inc i) (+ r n) (+ v (* n (+ (long (size a)) (long (size b)))))))))
        nr (long nr)
        ^ints va' (a/i32 nr) ^ints vb' (a/i32 nr)
        vals' (a/f64 (long nv))]
    (loop [i 0 s 0]
      (when (< i nk)
        (let [k (aget kl i) r0 (aget krow k) na (long (size (aget ka k))) nb (long (size (aget kb k)))
              s (long (loop [j 0 s s]
                        (if (== j (rows-of kind k))
                          s
                          (let [r (+ r0 j) r' (+ (aget krow' i) j)]
                            (aset va' r' (int s))
                            (copy! vals (aget va r) vals' s na)
                            (aset vb' r' (int (+ s na)))
                            (copy! vals (aget vb r) vals' (+ s na) nb)
                            (recur (inc j) (+ s na nb))))))]
          (recur (inc i) s))))
    {:nb nb' :nd nd :nk nk :nr nr
     :bsize bsize' :boff boff' :moff moff' :mass mass' :minv minv' :vstar vstar'
     :kind kind' :krow krow' :ka ka' :kb kb' :va va' :vb vb' :vals vals'
     :rn rn' :rt rt' :vhat vhat' :mu mu'}))

(defn solve
  "The velocities that minimize the step's cost, and the impulses that
  produce them.

  `problem` is flat arrays (see the comment above `rows-of`):

    :nb :bsize :boff   the blocks, and :nd their total size
    :mass :minv :moff  each block's mass matrix and its inverse, dense
    :vstar             where everything is going with nothing joining it
    :nk :kind :krow    the constraints: `contact`, `unilateral` or
                       `bilateral`, and where their rows start
    :ka :kb            the two blocks each joins, -1 for neither
    :va :vb :vals      each row's numbers against the first block and the
                       second
    :rn :rt            each constraint's regularization, normal and
                       tangent
    :vhat              the velocity each asks for along its first row
    :mu                each contact's friction coefficient

  and `v0` a first guess, or nil to start from `vstar`. Returns `{:v
  :impulses :iterations :cost :factors}`, the impulses three a constraint
  and the iterations those of the island that took most.

  `:factors` is each island's Newton factor, by its blocks. Given back
  as `opts`' `:factors` next step it preconditions the Newton steps
  there (see the comment above `hessian-mul!`), so a caller keeps it only
  while its blocks still mean the same bodies; a stale one costs time,
  never the answer.

  Each island -- blocks joined, however indirectly, by constraints -- is
  solved by itself. The cost is a sum over islands, so this is the same
  minimum, but a Newton step's line search is one length for the whole
  direction: solved together, a ragdoll struggling through its contacts
  shortened and lengthened the steps of bricks it was not touching, and
  they lifted off the floor and came back for twenty steps before the
  ragdoll was done. Apart, the bricks take one or two."
  ([problem] (solve problem nil {}))
  ([problem v0] (solve problem v0 {}))
  ([{:keys [^long nb ^long nd ^long nk ^ints ka ^ints kb ^ints boff ^ints bsize ^doubles vstar] :as problem} v0 opts]
   (let [parts (islands nb nk ka kb)]
     (if (and (== 1 (count parts)) (== nb (count (first (first parts)))))
       (let [blocks (first (first parts))
             r (newton problem v0 (assoc opts :factor (get (:factors opts) blocks)))]
         (-> r (dissoc :factor) (assoc :factors (if-let [f (:factor r)] {blocks f} {}))))
       (let [v (a/f64 nd) gam (a/f64 (* 3 nk))
             ^doubles v-init (or v0 vstar)]
         ;; anything nothing touches goes where it was going
         (dotimes [i nd] (aset v i (aget vstar i)))
         (reduce (fn [acc [blocks ks]]
                   (let [p (part problem blocks ks)
                         ^ints poff (:boff p)
                         ^ints bl (a/i32 blocks) ^ints kl (a/i32 ks)
                         guess (let [g (a/f64 (long (:nd p)))]
                                 (dotimes [i (alength bl)]
                                   (let [q (aget bl i) o (aget boff q) po (aget poff i)]
                                     (dotimes [j (aget bsize q)] (aset g (+ po j) (aget v-init (+ o j))))))
                                 g)
                         {:keys [^doubles impulses iterations cost failed factor] ^doubles pv :v}
                         (newton p guess (assoc opts :factor (get (:factors opts) blocks)))]
                     (dotimes [i (alength bl)]
                       (let [q (aget bl i) o (aget boff q) po (aget poff i)]
                         (dotimes [j (aget bsize q)] (aset v (+ o j) (aget pv (+ po j))))))
                     (dotimes [i (alength kl)]
                       (let [k3 (* 3 (aget kl i)) i3 (* 3 i)]
                         (dotimes [t 3] (aset gam (+ k3 t) (aget impulses (+ i3 t))))))
                     (cond-> (-> acc
                                 (update :iterations max iterations)
                                 (update :cost + cost))
                       factor (assoc-in [:factors blocks] factor)
                       failed (assoc :failed failed))))
                 {:v v :impulses gam :iterations 0 :cost 0.0 :factors {}}
                 parts))))))

(defn spd-inverse
  "The inverse of the symmetric positive definite `n`x`n` matrix `m`
  (flat, row-major), by Cholesky on arrays: what a model's mass matrix
  is recovered from its inverse with, every step. Nil if `m` is not
  positive definite."
  [m n]
  (let [n (long n)
        ^doubles l (a/f64 m)
        env {:first (a/i32 (repeat n 0))
             :start (a/i32 (map #(quot (* (long %) (inc (long %))) 2) (range (inc n))))}
        ;; the lower triangle, row by row, as `factor!` keeps an envelope
        ^doubles packed (a/f64 (quot (* n (inc n)) 2))]
    (dotimes [i n]
      (dotimes [j (inc i)]
        (aset packed (+ (quot (* i (inc i)) 2) j) (aget l (+ (* i n) j)))))
    (when (factor! packed env n)
      (let [out (a/f64 (* n n)) col (a/f64 n)]
        (dotimes [j n]
          (dotimes [i n] (aset col i (if (== i j) 1.0 0.0)))
          (substitute! packed env n col)
          (dotimes [i n] (aset out (+ (* i n) j) (aget col i))))
        out))))

;; ---------------------------------------------------------------------------
;; Building a problem

(defn builder
  "An empty problem to add blocks and constraints to, as a mutable map
  of growing vectors; `build` turns it into `solve`'s flat arrays."
  []
  {:bsize (transient []) :mass (transient []) :minv (transient []) :vstar (transient [])
   :kind (transient []) :ka (transient []) :kb (transient [])
   :rn (transient []) :rt (transient []) :vhat (transient []) :mu (transient [])
   :rows (transient [])})

(defn add-block
  "Adds a block of `n` numbers, its dense `mass` and `minv` (row-major,
  n*n numbers each) and `vstar`. Returns `[builder index]`."
  [b n mass minv vstar]
  (let [i (count (:bsize b))]
    [(-> b
         (update :bsize conj! (long n))
         (update :mass #(reduce conj! % mass))
         (update :minv #(reduce conj! % minv))
         (update :vstar #(reduce conj! % vstar)))
     i]))

(defn add-constraint
  "Adds a constraint of `kind` joining blocks `ka` and `kb` (-1 for none),
  with `rows` -- each `[against-a against-b]`, numbers as long as the
  block or nil -- and its `rn`, `rt`, `vhat` and `mu`."
  [b kind ka kb rows rn rt vhat mu]
  (-> b
      (update :kind conj! (long kind))
      (update :ka conj! (long ka))
      (update :kb conj! (long kb))
      (update :rn conj! (double rn))
      (update :rt conj! (double rt))
      (update :vhat conj! (double vhat))
      (update :mu conj! (double mu))
      (update :rows #(reduce conj! % rows))))

(defn build
  "`builder`'s growing vectors as `solve`'s flat arrays."
  [b]
  (let [bsize (persistent! (:bsize b))
        kind (persistent! (:kind b))
        rows (persistent! (:rows b))
        ka (persistent! (:ka b)) kb (persistent! (:kb b))
        nb (count bsize) nk (count kind)
        boff (vec (reductions + 0 bsize))
        moff (vec (reductions + 0 (map #(* % %) bsize)))
        krow (vec (reductions + 0 (map #(if (zero? (long %)) 3 1) kind)))
        nr (count rows)
        ;; each row's two stretches, laid end to end
        [va vb vals] (loop [r 0 va (transient []) vb (transient []) vals (transient [])]
                       (if (== r nr)
                         [(persistent! va) (persistent! vb) (persistent! vals)]
                         (let [[ra rb] (nth rows r)
                               s (count vals)
                               vals (reduce conj! vals (or ra []))
                               s2 (count vals)
                               vals (reduce conj! vals (or rb []))]
                           (recur (inc r) (conj! va s) (conj! vb s2) vals))))]
    {:nb nb :nd (long (nth boff nb)) :nk nk :nr nr
     :bsize (a/i32 bsize) :boff (a/i32 (butlast boff)) :moff (a/i32 (butlast moff))
     :mass (a/f64 (persistent! (:mass b))) :minv (a/f64 (persistent! (:minv b)))
     :vstar (a/f64 (persistent! (:vstar b)))
     :kind (a/i32 kind) :krow (a/i32 (butlast krow)) :ka (a/i32 ka) :kb (a/i32 kb)
     :va (a/i32 va) :vb (a/i32 vb) :vals (a/f64 vals)
     :rn (a/f64 (persistent! (:rn b))) :rt (a/f64 (persistent! (:rt b)))
     :vhat (a/f64 (persistent! (:vhat b))) :mu (a/f64 (persistent! (:mu b)))}))
