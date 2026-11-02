(ns allgo.vision.structure
  "Structure from motion: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 25 -- the shape of a scene and
  how the camera moved, both from a sequence of views of it, by Tomasi
  and Kanade's factorization (\"Shape and Motion from Image Streams under
  Orthography: a Factorization Method\", International Journal of
  Computer Vision 9, 1992). It solves together the chapter's relative
  orientation problem (how the views are turned to one another) and its
  depth estimation problem (where the points are), given the
  correspondences: which image point is which scene point in each view.

  Under orthography the image of point p in frame f is (i_f . p + a_f,
  j_f . p + b_f), i_f and j_f the camera's unit axes. So the measurements
  -- the u's of every frame stacked over the v's, a row a frame, a column
  a point -- less each row's mean, which takes out the translation, are
  the product of a 2F x 3 matrix of camera axes and a 3 x P matrix of
  points: rank three, whatever the number of frames and points. A step a
  stage:

    measured  the 2F x P matrix W of image coordinates
    centered  each row less its mean
    factored  the best rank-three factorization, M S, from the three
              largest singular values (here the eigenvectors of W W^T,
              by `allgo.numerics.eigen`'s Jacobi); `:residual` the size
              of what is left over, which noise makes nonzero
    metric    the factorization is only up to an invertible Q, M Q and
              Q^-1 S: Q is fixed by the camera axes being unit length and
              square to each other -- linear in L = Q Q^T, six unknowns,
              by least squares -- then the axes turned so the first
              frame's are x and y.

  What remains unknown is a mirror image: under orthography a shape and
  its reflection in the image plane, turning the other way, look the
  same. Frames are `[[u v] ...]`, a point each, every point in every
  frame."
  (:require [allgo.geometry.vec3 :as v]
            [allgo.numerics.eigen :as eigen]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(defn- columns-times
  "A^T B for matrices given as rows."
  [A B]
  (lin/mat-mul (lin/transpose A) B))

(defn- center [W]
  (mapv (fn [row] (let [m (/ (reduce + row) (count row))] (mapv #(- % m) row))) W))

(defn- factor
  "The rank-three factorization of `W`: `[M S residual]`, M = U3 sqrt(D3)
  and S = sqrt(D3)^-1 U3^T W, U3 the eigenvectors of W W^T for its three
  largest eigenvalues D3, the singular values squared; `residual` the
  root of the sum of the rest, the size of W - M S."
  [W]
  (let [{:keys [values vectors]} (eigen/jacobi (lin/mat-mul W (lin/transpose W)))
        n (count values)
        top (take-last 3 (range n))
        U3 (mapv (fn [row] (mapv row top)) vectors)
        roots (mapv #(math/sqrt (max 0.0 (double (values %)))) top)
        M (mapv (fn [row] (mapv #(* %1 (math/sqrt %2)) row roots)) U3)
        UtW (columns-times U3 W)
        S (mapv (fn [k] (mapv #(/ % (math/sqrt (roots k))) (UtW k))) (range 3))]
    [M S (math/sqrt (reduce + (map #(max 0.0 (double %)) (drop-last 3 values))))]))

(defn- quadratic-row
  "The coefficients of a^T L b in L's six unknowns, `[L11 L12 L13 L22 L23
  L33]`."
  [[a1 a2 a3] [b1 b2 b3]]
  [(* a1 b1) (+ (* a1 b2) (* a2 b1)) (+ (* a1 b3) (* a3 b1)) (* a2 b2) (+ (* a2 b3) (* a3 b2)) (* a3 b3)])

(defn- metric
  "M and S made Euclidean: L = Q Q^T from every frame's axes being unit and
  square, i^T L i = j^T L j = 1 and i^T L j = 0, in the least-squares
  sense; Q from L's eigenvectors, the eigenvalues clamped positive; then
  everything turned so the first frame's axes are x and y."
  [M S]
  (let [F (quot (count M) 2)
        rows (for [f (range F) :let [i (M f) j (M (+ f F))]
                   [a b rhs] [[i i 1.0] [j j 1.0] [i j 0.0]]]
               [(quadratic-row a b) rhs])
        A (mapv first rows) y (mapv second rows)
        [l11 l12 l13 l22 l23 l33] (lin/cholesky-solve (columns-times A A) (lin/mat-vec (lin/transpose A) y))
        L [[l11 l12 l13] [l12 l22 l23] [l13 l23 l33]]
        {:keys [values vectors]} (eigen/jacobi L)
        roots (mapv #(math/sqrt (max 1e-12 (double %))) values)
        ;; Q = V sqrt(D), V orthogonal, so Q^-1 = sqrt(D)^-1 V^T
        Q (mapv (fn [row] (mapv * row roots)) vectors)
        Qinv (mapv (fn [k] (mapv #(/ (double %) (double (roots k))) (mapv #(% k) vectors))) (range 3))
        M' (lin/mat-mul M Q)
        S' (lin/mat-mul Qinv S)
        ;; the first frame's axes, made exactly orthonormal, as the world's x and y
        i1 (v/normalize (M' 0))
        j1 (v/normalize (v/sub (M' F) (v/scale i1 (v/dot i1 (M' F)))))
        R [i1 j1 (v/cross i1 j1)]]
    [(lin/mat-mul M' (lin/transpose R)) (lin/mat-mul R S')]))

(defn factorization-seq
  "Tomasi and Kanade's factorization of `frames`, a stage a step: states
  `{:stage :W}`, then with `:M` (2F x 3, frame f's camera axes rows f and
  F + f) and `:S` (3 x P, a point a column), and `:residual` once
  factored."
  [frames]
  (let [F (count frames)
        W (into (mapv (fn [fr] (mapv first fr)) frames) (mapv (fn [fr] (mapv second fr)) frames))]
    (take-while some?
                (iterate (fn [{:keys [stage W M S] :as st}]
                           (case stage
                             :measured (assoc st :stage :centered :W (center W) :translation (mapv #(/ (reduce + %) (count %)) W))
                             :centered (let [[M S r] (factor W)] (assoc st :stage :factored :M M :S S :residual r))
                             :factored (let [[M S] (metric M S)] (assoc st :stage :metric :M M :S S))
                             :metric nil))
                         {:stage :measured :W W :frames F}))))

(defn factorize
  "The shape and motion `frames` show: `{:points [[x y z] ...] :axes [[i
  j] ...] :residual}`, the points in the first frame's camera
  coordinates (its image x and y, depth z) up to a reflection in depth,
  and each frame's camera axes."
  [frames]
  (let [{:keys [M S frames residual]} (last (factorization-seq frames))]
    {:points (lin/transpose S)
     :axes (mapv (fn [f] [(M f) (M (+ f frames))]) (range frames))
     :residual residual}))
