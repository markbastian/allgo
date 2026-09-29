(ns allgo.numerics.fem
  "The finite-element method (Chapra and Canale, *Numerical Methods for
  Engineers*, chapter 31), in its five steps: the domain cut into
  elements; on each, the solution approximated by linear shape functions
  and the element equations derived from the differential equation by
  Galerkin's weighting; the element equations assembled into the system's
  by the topology; the boundary conditions imposed; the system solved; and
  the postprocessing that recovers the derived quantities.

  In one dimension, k u'' + f(x) = 0 -- the book's heated rod -- on line
  elements; in two, the Poisson equation -k (u_xx + u_yy) = f(x, y) on
  triangles with linear shape functions."
  (:require [allgo.numerics.linear-systems :as ls]
            [allgo.numerics.quadrature :as quad]))

(defn- impose
  "The assembled system `K u = F` with the nodes of `fixed`, a map of node
  to value, held there: their columns moved to the right side and their
  rows replaced by u_i = value."
  [K F fixed]
  (let [n (count F)
        F (reduce (fn [F [j v]] (mapv (fn [i fi] (if (contains? fixed i) fi (- fi (* (get-in K [i j]) v)))) (range n) F))
                  F fixed)
        K (reduce (fn [K [j _]]
                    (-> (mapv #(assoc % j 0.0) K)
                        (assoc j (assoc (vec (repeat n 0.0)) j 1.0))))
                  K fixed)
        F (reduce (fn [F [j v]] (assoc F j (double v))) F fixed)]
    [K F]))

(defn line-elements
  "k u'' + f(x) = 0 on the nodes `xs` (ascending; the elements between
  them) by linear finite elements: each element's equations (k/L) [[1 -1]
  [-1 1]] {u} = {int f N} + the end fluxes, the source integrated by
  Gauss-Legendre quadrature, assembled and solved. The ends are held at
  `:left` and `:right` (values), or given `:left-flux` / `:right-flux` --
  k du/dx there. Returns `{:u :left-flux :right-flux}`, the nodal values
  and, from the first and last element's equations, the end fluxes k
  du/dx (the book's postprocessing, eq. 31.26)."
  [xs {:keys [k f left right left-flux right-flux] :or {k 1.0 f (constantly 0.0)}}]
  (let [xs (mapv double xs) n (count xs)
        elements (mapv (fn [e] [e (inc e)]) (range (dec n)))
        zero-K (vec (repeat n (vec (repeat n 0.0))))
        [K F] (reduce (fn [[K F] [a b]]
                        (let [x1 (xs a) x2 (xs b) L (- x2 x1) c (/ k L)
                              fa (quad/gauss-legendre #(* (f %) (/ (- x2 %) L)) x1 x2 5)
                              fb (quad/gauss-legendre #(* (f %) (/ (- % x1) L)) x1 x2 5)]
                          [(-> K (update-in [a a] + c) (update-in [a b] - c) (update-in [b a] - c) (update-in [b b] + c))
                           (-> F (update a + fa) (update b + fb))]))
                      [zero-K (vec (repeat n 0.0))] elements)
        ;; a given flux enters the end's equation; k du/dx at the left end
        ;; is -(the outward flux term), at the right +
        F (cond-> F left-flux (update 0 - left-flux) right-flux (update (dec n) + right-flux))
        fixed (cond-> {} left (assoc 0 left) right (assoc (dec n) right))
        [K' F'] (impose K F fixed)
        u (ls/gauss K' F')
        ;; postprocessing: the end rows of the unconstrained system give the
        ;; fluxes -- (K u - F) at the first node is -k u'(x1), at the last
        ;; +k u'(xn)
        residual (fn [i] (- (reduce + (map * (K i) u)) (F i)))]
    {:u u
     :left-flux (or left-flux (- (residual 0)))
     :right-flux (or right-flux (residual (dec n)))}))

(defn triangles
  "-k (u_xx + u_yy) = f(x, y) by linear triangular elements: `nodes` the
  points `[[x y] ...]`, `elements` the triangles as triples of node
  indices (counterclockwise or not), `fixed` a map of node index to its
  value. Each triangle's stiffness k A (grad N_i . grad N_j), A its area,
  and load f at its centroid times A/3 per node; assembled, constrained
  and solved. Returns the nodal values."
  [nodes elements {:keys [k f fixed] :or {k 1.0 f (constantly 0.0) fixed {}}}]
  (let [n (count nodes)
        zero-K (vec (repeat n (vec (repeat n 0.0))))
        [K F] (reduce (fn [[K F] [i j m]]
                        (let [[xi yi] (nodes i) [xj yj] (nodes j) [xm ym] (nodes m)
                              area2 (- (* (- xj xi) (- ym yi)) (* (- xm xi) (- yj yi)))
                              A (* 0.5 (abs area2))
                              ;; the gradients of the shape functions: b/(2A), c/(2A)
                              bs [(- yj ym) (- ym yi) (- yi yj)]
                              cs [(- xm xj) (- xi xm) (- xj xi)]
                              idx [i j m]
                              load (* (f (/ (+ xi xj xm) 3.0) (/ (+ yi yj ym) 3.0)) (/ A 3.0))]
                          [(reduce (fn [K [p q]]
                                     (update-in K [(idx p) (idx q)] +
                                                (/ (* k (+ (* (bs p) (bs q)) (* (cs p) (cs q)))) (* 4.0 A))))
                                   K (for [p (range 3) q (range 3)] [p q]))
                           (reduce (fn [F p] (update F (idx p) + load)) F (range 3))]))
                      [zero-K (vec (repeat n 0.0))] elements)
        [K' F'] (impose K F fixed)]
    (ls/gauss K' F')))
