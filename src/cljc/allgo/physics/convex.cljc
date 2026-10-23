(ns allgo.physics.convex
  "Contact as one convex problem, solved by Newton's method: the shape of
  MuJoCo's solver (Todorov, \"Convex and analytically-invertible dynamics
  with contacts and constraints\", ICRA 2014) and of Drake's SAP (Castro,
  Han and Permenter, \"An unconstrained convex formulation of compliant
  contact\", IEEE T-RO 2022), whose formulation this follows.

  ## Why

  `allgo.physics.solver`'s three solvers visit contacts one at a time,
  each correcting itself against its neighbors' current answers. What the
  floor knows reaches the top of a column one contact per sweep, so a
  column taller than its iterations can carry never learns its own
  weight, and leans. Here every contact is solved at once: each Newton
  step is one linear system coupling every body to every other through
  the contacts, and it converges in a handful of steps however tall the
  stack.

  ## The problem

  The velocities `v` at the end of the step minimize

      l(v) = 1/2 (v - v*)' M (v - v*)  +  sum over contacts of l_c

  where `v*` is where the bodies would be going with nothing touching
  and `M` their masses. A contact's own term is `1/2 g' R g`, with its
  impulse `g` the projection of `y = -R^-1 (vc - vhat)` onto its friction
  cone in `R`'s norm: `vc` its relative velocity -- separating along the
  normal, sliding along the two tangents -- `vhat` the velocity it asks
  for, which pushes an overlap out, and `R` a small regularization that
  makes the contact very stiff rather than rigid. The gradient is
  `M (v - v*) - J' g` and the Hessian `M + J' G J`, `G` the derivative
  of the impulse, so a Newton step is one solve with that matrix, and a
  backtracking line search keeps every step a descent.

  `l` is convex whatever the friction, which is what the cone in `R`'s
  norm buys over the usual complementarity problem: there is one answer
  and Newton finds it.

  ## What this namespace knows

  Only the arithmetic, on flat arrays: a problem of bodies' degrees of
  freedom and contacts' Jacobian rows, as `solve` describes. Where those
  come from -- collision detection, warm starting, sleeping, moving the
  bodies afterward -- is `allgo.physics.solver`'s, and this is one more
  of its `:solver`s."
  (:require [allgo.array :as a]
            [clojure.math :as math]))

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

;; ---------------------------------------------------------------------------
;; The Newton matrix, stored by its envelope
;;
;; Two bodies are coupled only if they touch, so the matrix is almost all
;; zeros, and factoring it densely -- 1,536 rows for a 16x16 wall -- took
;; 600 ms a step. Numbered so that bodies that touch are near each other
;; (reverse Cuthill-McKee), the nonzeros and everything the factor fills
;; in lie close to the diagonal: each row is kept from its first nonzero
;; to the diagonal and no further, and a column of bricks is a band six
;; or twelve wide however tall it is.

(defn- ordering
  "Reverse Cuthill-McKee over the `nb` bodies, touching where a contact
  joins two that move: each body's new position."
  [nb nc ^ints oa ^ints ob]
  (let [nb (long nb)
        adj (let [adj (object-array nb)]
              (dotimes [b nb] (aset adj b (transient #{})))
              (dotimes [c (long nc)]
                (let [a (aget oa c) b (aget ob c)]
                  (when (and (>= a 0) (>= b 0) (not= a b))
                    (let [a (quot a 6) b (quot b 6)]
                      (aset adj a (conj! (aget adj a) b))
                      (aset adj b (conj! (aget adj b) a))))))
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
  "Where each row of the matrix starts and is stored, for bodies placed
  by `perm` and joined as `adj`: a row of body `p`'s block begins at the
  first column of the lowest-numbered body it touches, or its own."
  [nb ^ints perm adj]
  (let [nb (long nb) nd (* 6 nb)
        ^ints first-col (a/i32 nd)
        ^ints start (a/i32 (inc nd))]
    (dotimes [b nb]
      (let [p (aget perm b)
            lowest (reduce (fn [m q] (min (long m) (long (aget perm (long q))))) p (nth adj b))]
        (dotimes [r 6] (aset first-col (+ (* 6 p) r) (int (* 6 (long lowest)))))))
    (dotimes [i nd]
      (aset start (inc i) (int (+ (aget start i) (- (inc i) (aget first-col i))))))
    {:first first-col :start start :size (aget start nd)}))

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
      (if (= i nd)
        true
        (let [fi (aget first i)
              ri (long (row i))
              ok? (loop [j fi]
                    (if (> j i)
                      true
                      (let [rj (long (row j))
                            k0 (max fi (aget first j))
                            s (loop [k k0 s (aget h (+ ri j))]
                                (if (= k j) s (recur (inc k) (- s (* (aget h (+ ri k)) (aget h (+ rj k)))))))]
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
  (let [nd (long nd)]
    (dotimes [i nd]
      (let [s (loop [k (aget first i) s (aget b i)]
                (if (= k i) s (recur (inc k) (- s (* (aget h (at first start i k)) (aget b k))))))]
        (aset b i (/ s (aget h (at first start i i))))))
    (loop [i (dec nd)]
      (when (>= i 0)
        (let [x (/ (aget b i) (aget h (at first start i i)))]
          (aset b i x)
          (loop [k (aget first i)]
            (when (< k i)
              (aset b k (- (aget b k) (* (aget h (at first start i k)) x)))
              (recur (inc k))))
          (recur (dec i)))))
    b))

;; ---------------------------------------------------------------------------
;; The problem

(defn- mass-mul!
  "`out` = M `x`: each body's six degrees of freedom, three linear by its
  mass and three angular by its world inertia."
  [^doubles mass ^doubles inertia nb ^doubles x ^doubles out]
  (dotimes [b (long nb)]
    (let [o (* 6 b) m (aget mass b) i9 (* 9 b)]
      (aset out o (* m (aget x o)))
      (aset out (+ o 1) (* m (aget x (+ o 1))))
      (aset out (+ o 2) (* m (aget x (+ o 2))))
      (dotimes [r 3]
        (aset out (+ o 3 r) (+ (* (aget inertia (+ i9 (* 3 r))) (aget x (+ o 3)))
                               (* (aget inertia (+ i9 (* 3 r) 1)) (aget x (+ o 4)))
                               (* (aget inertia (+ i9 (* 3 r) 2)) (aget x (+ o 5))))))))
  out)

(defn- row-dot
  "One Jacobian row of contact `c`, direction `k`, against `v`."
  [^doubles jac ^ints oa ^ints ob c k ^doubles v]
  (let [c (long c) k (long k)
        base (* 36 c) r (+ base (* 12 k))
        a (aget oa c) b (aget ob c)]
    (+ (if (neg? a) 0.0
           (loop [j 0 s 0.0] (if (= j 6) s (recur (inc j) (+ s (* (aget jac (+ r j)) (aget v (+ a j))))))))
       (if (neg? b) 0.0
           (loop [j 0 s 0.0] (if (= j 6) s (recur (inc j) (+ s (* (aget jac (+ r 6 j)) (aget v (+ b j)))))))))))

(defn- evaluate!
  "The cost at `v`, with each contact's impulse into `gam` and, when
  `gms` is given, its `G` into `gms`."
  [{:keys [^long nb ^long nc ^doubles mass ^doubles inertia ^doubles vstar
           ^doubles jac ^ints oa ^ints ob ^doubles rn ^doubles rt ^doubles vhat ^doubles mu]}
   ^doubles v ^doubles gam gms ^doubles scratch]
  (let [nd (* 6 nb)
        dv (a/f64 nd)
        y (a/f64 3) g (a/f64 3) gm (a/f64 9)]
    (dotimes [i nd] (aset dv i (- (aget v i) (aget vstar i))))
    (mass-mul! mass inertia nb dv scratch)
    (let [kinetic (loop [i 0 s 0.0] (if (= i nd) s (recur (inc i) (+ s (* (aget dv i) (aget scratch i))))))]
      (loop [c 0 total (* 0.5 kinetic)]
        (if (= c nc)
          total
          (let [vn (row-dot jac oa ob c 0 v)
                v1 (row-dot jac oa ob c 1 v)
                v2 (row-dot jac oa ob c 2 v)]
            (aset y 0 (/ (- (aget vhat c) vn) (aget rn c)))
            (aset y 1 (/ (- v1) (aget rt c)))
            (aset y 2 (/ (- v2) (aget rt c)))
            (let [cost (project! (aget mu c) (aget rn c) (aget rt c) y g gm)
                  c3 (* 3 c)]
              (aset gam c3 (aget g 0)) (aset gam (+ c3 1) (aget g 1)) (aset gam (+ c3 2) (aget g 2))
              (when gms
                (let [^doubles gms gms c9 (* 9 c)]
                  (dotimes [i 9] (aset gms (+ c9 i) (aget gm i)))))
              (recur (inc c) (+ total cost)))))))))

(defn- add-jt-gamma!
  "`out` += J' `gam`."
  [{:keys [^long nc ^doubles jac ^ints oa ^ints ob]} ^doubles gam ^doubles out]
  (dotimes [c nc]
    (let [a (aget oa c) b (aget ob c) c3 (* 3 c)]
      (dotimes [k 3]
        (let [gk (aget gam (+ c3 k)) r (+ (* 36 c) (* 12 k))]
          (when-not (neg? a)
            (dotimes [j 6] (aset out (+ a j) (+ (aget out (+ a j)) (* gk (aget jac (+ r j)))))))
          (when-not (neg? b)
            (dotimes [j 6] (aset out (+ b j) (+ (aget out (+ b j)) (* gk (aget jac (+ r 6 j))))))))))))

(defn- hessian!
  "`h` = M + J' G J, its lower triangle into the envelope, rows and
  columns renumbered by `perm`."
  [{:keys [^long nb ^long nc ^doubles mass ^doubles inertia ^doubles jac ^ints oa ^ints ob]}
   ^doubles gms ^ints perm {:keys [^ints first ^ints start size] :as _env} ^doubles h]
  (dotimes [i (long size)] (aset h i 0.0))
  (dotimes [b nb]
    (let [o (* 6 (aget perm b)) i9 (* 9 b)]
      (dotimes [r 3]
        (aset h (at first start (+ o r) (+ o r)) (aget mass b))
        (dotimes [s (inc r)]
          (let [idx (at first start (+ o 3 r) (+ o 3 s))]
            (aset h idx (+ (aget h idx) (aget inertia (+ i9 (* 3 r) s)))))))))
  (let [gj (a/f64 36)]
    (dotimes [c nc]
      (let [c9 (* 9 c) base (* 36 c)
            a (aget oa c) b (aget ob c)]
        ;; G J for this contact: three rows of twelve
        (dotimes [k 3]
          (dotimes [j 12]
            (aset gj (+ (* 12 k) j)
                  (+ (* (aget gms (+ c9 (* 3 k))) (aget jac (+ base j)))
                     (* (aget gms (+ c9 (* 3 k) 1)) (aget jac (+ base 12 j)))
                     (* (aget gms (+ c9 (* 3 k) 2)) (aget jac (+ base 24 j)))))))
        ;; J' (G J), into the blocks of the bodies it touches, below the
        ;; diagonal only
        (dotimes [side-p 2]
          (dotimes [side-q 2]
            (let [p (if (zero? side-p) a b) q (if (zero? side-q) a b)]
              (when (and (>= p 0) (>= q 0))
                (let [pn (* 6 (aget perm (quot p 6))) qn (* 6 (aget perm (quot q 6)))
                      pj (* 6 side-p) qj (* 6 side-q)]
                  (when (>= pn qn)
                    (dotimes [r 6]
                      (dotimes [s 6]
                        (when (or (> pn qn) (>= r s))
                          (let [idx (at first start (+ pn r) (+ qn s))]
                            (aset h idx (+ (aget h idx)
                                           (* (aget jac (+ base pj r)) (aget gj (+ qj s)))
                                           (* (aget jac (+ base 12 pj r)) (aget gj (+ 12 qj s)))
                                           (* (aget jac (+ base 24 pj r)) (aget gj (+ 24 qj s))))))))))))))))))
  h)

(defn pushed
  "Where the bodies would be going if `impulses` -- three a contact, as
  `solve` returns them -- were applied to `vstar`: `vstar + M^-1 J' g`.
  Last step's impulses through this are the warm start, the guess the
  Newton steps begin from. `problem` needs `:inv-mass` and
  `:inv-inertia`, the inverses of `:mass` and `:inertia`."
  [{:keys [^long nb ^doubles vstar ^doubles inv-mass ^doubles inv-inertia] :as problem} ^doubles impulses]
  (let [nd (* 6 nb)
        f (a/f64 nd)
        v (a/f64 nd)]
    (add-jt-gamma! problem impulses f)
    (dotimes [b nb]
      (let [o (* 6 b) im (aget inv-mass b) i9 (* 9 b)]
        (dotimes [r 3]
          (aset v (+ o r) (+ (aget vstar (+ o r)) (* im (aget f (+ o r)))))
          (aset v (+ o 3 r) (+ (aget vstar (+ o 3 r))
                               (* (aget inv-inertia (+ i9 (* 3 r))) (aget f (+ o 3)))
                               (* (aget inv-inertia (+ i9 (* 3 r) 1)) (aget f (+ o 4)))
                               (* (aget inv-inertia (+ i9 (* 3 r) 2)) (aget f (+ o 5))))))))
    v))

(defn cost
  "The step's cost at `v`."
  ^double [{:keys [^long nb ^long nc] :as problem} ^doubles v]
  (evaluate! problem v (a/f64 (* 3 nc)) nil (a/f64 (* 6 nb))))

(defn solve
  "The velocities that minimize the step's cost, and the impulses that
  produce them.

  `problem` is flat arrays:

    :nb       bodies that move, each six degrees of freedom: linear
              velocity, then angular, in the world
    :mass     each one's mass, and :inertia its world inertia tensor,
              nine numbers row-major
    :vstar    where each is going with nothing touching it, 6 nb
    :nc       contacts
    :oa :ob   where each contact's two bodies start in the velocity
              vector, or -1 for a body that does not move
    :jac      each contact's three rows -- normal, then two tangents --
              of twelve: six for its first body, six for its second
    :rn :rt   the regularization along the normal and the tangents
    :vhat     the separating velocity each asks for along its normal
    :mu       each one's friction coefficient

  and `v0` a first guess, or nil to start from `vstar`. Returns `{:v
  :impulses :iterations :cost}`, the impulses three a contact in the
  order of `jac`'s rows."
  ([problem] (solve problem nil {}))
  ([problem v0] (solve problem v0 {}))
  ([{:keys [^long nb ^long nc ^doubles mass ^doubles inertia ^doubles vstar ^ints oa ^ints ob] :as problem}
    v0 {:keys [max-iterations tolerance] :or {max-iterations 30 tolerance 1e-10}}]
   (let [nd (* 6 nb)
         v (a/f64 nd)
         ^doubles v-init (or v0 vstar)
         _ (dotimes [i nd] (aset v i (aget v-init i)))
         gam (a/f64 (* 3 nc))
         gms (a/f64 (* 9 nc))
         {:keys [^ints perm adj]} (ordering nb nc oa ob)
         env (envelope nb perm adj)
         h (a/f64 (long (:size env)))
         pdir (a/f64 nd)
         grad (a/f64 nd)
         scratch (a/f64 nd)
         trial (a/f64 nd)
         trial-gam (a/f64 (* 3 nc))]
     (loop [it 0 cost (evaluate! problem v gam gms scratch)]
       ;; the gradient, M (v - v*) - J' g
       (dotimes [i nd] (aset trial i (- (aget v i) (aget vstar i))))
       (mass-mul! mass inertia nb trial grad)
       (dotimes [i nd] (aset grad i (- (aget grad i))))
       (add-jt-gamma! problem gam grad)
       (dotimes [i nd] (aset grad i (- (aget grad i))))
       (let [gnorm (loop [i 0 s 0.0] (if (= i nd) s (recur (inc i) (max s (abs (aget grad i))))))]
         (if (or (>= it (long max-iterations)) (<= gnorm (double tolerance)) (zero? nd))
           {:v v :impulses gam :iterations it :cost cost}
           (do
             (hessian! problem gms perm env h)
             (if-not (factor! h env nd)
               {:v v :impulses gam :iterations it :cost cost :failed :not-positive-definite}
               (let [dir (a/f64 nd)
                     ;; renumbered into the envelope's order and back
                     _ (dotimes [b nb]
                         (let [o (* 6 b) p (* 6 (aget perm b))]
                           (dotimes [j 6] (aset pdir (+ p j) (- (aget grad (+ o j)))))))
                     _ (substitute! h env nd pdir)
                     _ (dotimes [b nb]
                         (let [o (* 6 b) p (* 6 (aget perm b))]
                           (dotimes [j 6] (aset dir (+ o j) (aget pdir (+ p j))))))
                     slope (loop [i 0 s 0.0] (if (= i nd) s (recur (inc i) (+ s (* (aget grad i) (aget dir i))))))
                     ;; Armijo backtracking: the cost is convex, so a short
                     ;; enough step along a descent direction always does
                     [alpha _] (loop [alpha 1.0 tries 0]
                                 (dotimes [i nd] (aset trial i (+ (aget v i) (* alpha (aget dir i)))))
                                 (let [c' (evaluate! problem trial trial-gam nil scratch)]
                                   (if (or (<= c' (+ cost (* 1e-4 alpha slope))) (>= tries 30))
                                     [alpha c']
                                     (recur (* 0.5 alpha) (inc tries)))))
                     step (loop [i 0 s 0.0] (if (= i nd) s (recur (inc i) (max s (abs (* (double alpha) (aget dir i)))))))]
                 (dotimes [i nd] (aset v i (aget trial i)))
                 (let [cost' (evaluate! problem v gam gms scratch)]
                   (if (< step 1e-12)
                     {:v v :impulses gam :iterations (inc it) :cost cost'}
                     (recur (inc it) cost'))))))))))))
