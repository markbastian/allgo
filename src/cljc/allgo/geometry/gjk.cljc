(ns allgo.geometry.gjk
  "Gilbert-Johnson-Keerthi distance between convex sets, and EPA for
  penetration depth once they overlap. After Gino van den Bergen,
  \"Collision Detection in Interactive 3D Environments\".

  A convex set is given only by a *support mapping*: a function from a
  direction to the point of the set farthest along it. Nothing else about
  the shape is needed, so the same code serves point clouds, spheres, boxes
  and their Minkowski sums.

  GJK works on the Minkowski difference A - B, where the distance between
  the sets equals the distance from that body to the origin. It walks a
  simplex of at most four points toward the origin, at each step keeping the
  sub-simplex nearest to it, until no support point makes further progress.

  The sub-simplex step is Johnson's algorithm in van den Bergen's account.
  This uses the equivalent Voronoi-region formulation instead: the two agree
  on the answer, but Johnson's recursion over subset determinants is the
  part of GJK he singles out as numerically delicate, since it can affirm a
  sub-simplex it has already rejected and stall the loop. The region tests
  need no such consistency between separately computed determinants."
  (:require [allgo.geometry.vec3 :as v]
            [allgo.math :as am]
            [clojure.math :as math]))

(def ^:private origin [0.0 0.0 0.0])

;; ---------------------------------------------------------------- supports

(defn point-cloud
  "Support mapping for the convex hull of `points`."
  [points]
  (let [pts (vec points)]
    (fn [dir] (apply max-key #(v/dot % dir) pts))))

(defn sphere
  "Support mapping for a ball."
  [center radius]
  (fn [dir] (v/add center (v/scale (v/normalize dir) radius))))

(defn box
  "Support mapping for an axis-aligned box."
  [[x0 y0 z0] [x1 y1 z1]]
  (fn [[dx dy dz]]
    [(if (pos? dx) x1 x0) (if (pos? dy) y1 y0) (if (pos? dz) z1 z0)]))

(defn translate
  "Shift a support mapping by `offset`."
  [support offset]
  (fn [dir] (v/add (support dir) offset)))

;; --------------------------------------------- closest point on a simplex

(defn- weighted [pts weights]
  (reduce (fn [acc [i w]] (v/add acc (v/scale (nth pts i) w))) origin weights))

(defn- closest-seg [pts i j]
  (let [a  (nth pts i) b (nth pts j)
        ab (v/sub b a)
        dd (v/length-squared ab)]
    (if (< dd 1e-20)
      {:weights {i 1.0}}
      (let [t (am/clamp (/ (- (v/dot a ab)) dd) 0.0 1.0)]
        {:weights {i (- 1.0 t) j t}}))))

(defn- best [cands pts]
  (apply min-key #(v/length-squared (weighted pts (:weights %))) cands))

(defn- closest-tri
  "Nearest point of triangle `i j k` to the origin, as barycentric weights."
  [pts i j k]
  (let [a (nth pts i) b (nth pts j) c (nth pts k)
        n (v/cross (v/sub b a) (v/sub c a))
        nn (v/length-squared n)]
    (if (< nn 1e-20)
      (best [(closest-seg pts i j) (closest-seg pts i k) (closest-seg pts j k)] pts)
      (let [u (/ (v/dot (v/cross b c) n) nn)
            v (/ (v/dot (v/cross c a) n) nn)
            w (/ (v/dot (v/cross a b) n) nn)]
        (if (and (>= u -1e-12) (>= v -1e-12) (>= w -1e-12))
          {:weights {i u j v k w}}
          (best [(closest-seg pts i j) (closest-seg pts j k) (closest-seg pts i k)] pts))))))

(defn- outside-face?
  "Is the origin on the far side of plane `a b c` from `d`?"
  [a b c d]
  (let [n  (v/cross (v/sub b a) (v/sub c a))
        so (v/dot n (v/sub origin a))
        sd (v/dot n (v/sub d a))]
    (neg? (* so sd))))

(defn- closest-simplex
  "The point of `pts` (1 to 4 points) nearest the origin, given as weights
  over the sub-simplex that supports it. `:inside` marks the origin enclosed
  by a full tetrahedron -- the sets overlap and GJK is done."
  [pts]
  (case (count pts)
    1 {:weights {0 1.0}}
    2 (closest-seg pts 0 1)
    3 (closest-tri pts 0 1 2)
    4 (let [[a b c d] pts
            vol   (abs (v/dot (v/cross (v/sub b a) (v/sub c a)) (v/sub d a)))
            scale (max (v/length-squared (v/sub b a)) (v/length-squared (v/sub c a)) (v/length-squared (v/sub d a)) 1e-30)
            faces [(closest-tri pts 0 1 2) (closest-tri pts 0 1 3)
                   (closest-tri pts 0 2 3) (closest-tri pts 1 2 3)]]
        ;; A near-flat tetrahedron makes every same-side test ambiguous, and
        ;; concluding "enclosed" from four ambiguous signs reports an overlap
        ;; that is not there. It happens whenever the Minkowski difference is
        ;; a thin slab -- two broad faces almost touching -- which is exactly
        ;; the case where the true distance is small but nonzero. Only trust
        ;; enclosure from a tetrahedron with real volume.
        (if (and (> vol (* 1e-10 (math/pow scale 1.5)))
                 (not (outside-face? a b c d))
                 (not (outside-face? a b d c))
                 (not (outside-face? a c d b))
                 (not (outside-face? b c d a)))
          {:weights {} :inside true}
          (best faces pts)))))

;; ------------------------------------------------------------------- GJK

(def ^:private max-iterations 64)
(def ^:private epa-iterations 96)
(def ^:private epa-tolerance 1e-7)

(defn distance
  "Distance between convex sets `sa` and `sb`, each a support mapping.

  Returns `{:distance d :direction v :on-a p :on-b q}`, where `v` runs from
  the nearest point of A to that of B. When the sets overlap the distance is
  0 and `:simplex` holds the enclosing tetrahedron, which `penetration` can
  expand. Witness points are recovered by carrying, for each simplex vertex,
  the supports on A and B that produced it."
  [sa sb]
  (let [sup (fn [dir] (let [a (sa dir) b (sb (v/scale dir -1.0))]
                        {:w (v/sub a b) :a a :b b}))]
    (loop [simplex [(sup [1.0 0.0 0.0])] i 0]
      (let [pts    (mapv :w simplex)
            {:keys [weights inside]} (closest-simplex pts)]
        (if inside
          {:distance 0.0 :direction origin :on-a origin :on-b origin
           :simplex simplex :overlap? true}
          (let [v     (weighted pts weights)
                kept  (mapv #(nth simplex %) (sort (keys weights)))
                remap (into {} (map-indexed (fn [n idx] [n (weights idx)])
                                            (sort (keys weights))))
                vv    (v/length-squared v)]
            (if (< vv 1e-20)
              {:distance 0.0 :direction origin :overlap? true :simplex simplex
               :on-a origin :on-b origin}
              (let [w (sup (v/scale v -1.0))]
                ;; No progress toward the origin means v is already the
                ;; minimum-norm point of the whole Minkowski difference.
                ;; Exhausting the iteration budget means v has stopped
                ;; improving measurably, so it is the answer to report --
                ;; not grounds to claim an overlap that was never found.
                (if (or (>= i max-iterations)
                        (<= (- vv (v/dot v (:w w))) (* 1e-12 vv))
                        (some #(< (v/length-squared (v/sub (:w %) (:w w))) 1e-20) kept))
                  {:distance  (v/length v)
                   :direction (v/scale (v/normalize v) -1.0)
                   :on-a      (reduce (fn [acc [n wt]] (v/add acc (v/scale (:a (nth kept n)) wt))) origin remap)
                   :on-b      (reduce (fn [acc [n wt]] (v/add acc (v/scale (:b (nth kept n)) wt))) origin remap)
                   :overlap?  false
                   :simplex   kept}
                  (recur (conj kept w) (inc i)))))))))))

(defn intersects?
  "Whether two convex sets share a point."
  [sa sb]
  (boolean (:overlap? (distance sa sb))))

;; ------------------------------------------------------------------- EPA

(defn- face-of
  "One face of the expanding polytope, with an outward normal and the signed
  distance from the origin to its plane.

  Outwardness is judged against `interior`, a point known to be inside the
  polytope, rather than against the origin. The origin frequently lies
  exactly *on* a face -- shapes that touch flush, which is most of them --
  and the sign of its offset then says nothing about which way the face
  points, leaving EPA to expand inward and report a depth of zero."
  [pts [i j k] interior]
  (let [a (nth pts i) b (nth pts j) c (nth pts k)
        n (v/cross (v/sub b a) (v/sub c a))
        m (v/length n)]
    (when (> m 1e-12)
      (let [u        (v/scale n (/ m))
            outward? (pos? (v/dot u (v/sub a interior)))
            u        (if outward? u (v/scale u -1.0))]
        {:idx (if outward? [i j k] [i k j]) :normal u :dist (v/dot u a)}))))

(defn- expand-to-tetra
  "Grow a degenerate simplex into a tetrahedron enclosing the origin.

  GJK stops as soon as it has confirmed overlap, which often leaves a point,
  segment or triangle -- the origin lying exactly on a face or edge of the
  Minkowski difference. EPA needs a closed volume to expand, so missing
  dimensions are recovered by taking supports along directions orthogonal to
  what is already there."
  [sup pts]
  (loop [pts (vec pts)]
    (case (count pts)
      4 pts
      1 (let [p (first pts)]
          (when-let [q (some (fn [d] (let [s (sup d)]
                                       (when (> (v/length-squared (v/sub s p)) 1e-12) s)))
                             [[1.0 0 0] [-1.0 0 0] [0 1.0 0] [0 -1.0 0] [0 0 1.0] [0 0 -1.0]])]
            (recur [p q])))
      2 (let [[a b] pts
              ab   (v/sub b a)
              axis (if (< (abs (nth (v/normalize ab) 0)) 0.9) [1.0 0 0] [0 1.0 0])
              d    (v/cross ab axis)]
          (when-let [q (some (fn [dir] (let [s (sup dir)]
                                         (when (> (v/length-squared (v/cross (v/sub s a) ab)) 1e-12) s)))
                             [d (v/scale d -1.0) (v/cross ab d) (v/scale (v/cross ab d) -1.0)])]
            (recur [a b q])))
      3 (let [[a b c] pts
              n (v/normalize (v/cross (v/sub b a) (v/sub c a)))]
          (when-let [q (some (fn [dir] (let [s (sup dir)]
                                         (when (> (abs (v/dot (v/sub s a) n)) 1e-9) s)))
                             [n (v/scale n -1.0)])]
            (recur [a b c q]))))))

(defn penetration
  "Penetration depth and direction for overlapping convex sets, by the
  Expanding Polytope Algorithm: grow the tetrahedron GJK finished with
  toward the boundary of the Minkowski difference, always through the face
  nearest the origin, until no support point lies measurably beyond it.

  Returns `{:depth d :normal n}` -- translating B by `n * d` separates them.
  nil when the sets do not actually overlap."
  [sa sb]
  (let [{:keys [simplex overlap?]} (distance sa sb)
        sup (fn [dir] (v/sub (sa dir) (sb (v/scale dir -1.0))))]
    (when overlap?
      (when-let [tetra (expand-to-tetra sup (mapv :w simplex))]
        (let [interior (v/scale (reduce v/add tetra) 0.25)]
          (loop [pts   tetra
                 faces (keep #(face-of tetra % interior) [[0 1 2] [0 1 3] [0 2 3] [1 2 3]])
                 i     0]
            (if (empty? faces)
              {:depth 0.0 :normal [0.0 1.0 0.0]}
              (let [{:keys [normal dist]} (apply min-key :dist faces)
                    p    (sup normal)
                    d    (v/dot p normal)]
              ;; The nearest face is always a lower bound on the depth and
              ;; rises toward it, so running out of iterations still yields
              ;; the best answer so far. Smooth shapes never terminate on
              ;; tolerance -- every step just refines the polytope
              ;; approximating the surface -- and reporting zero there would
              ;; discard a converging result.
                (if (or (<= (- d dist) epa-tolerance) (>= i epa-iterations))
                  {:depth dist :normal normal}
                ;; Remove every face the new point can see, then re-close the
                ;; hull over the horizon those faces left behind.
                  (let [n     (count pts)
                        pts'  (conj pts p)
                        seen  (filter #(> (- (v/dot p (:normal %)) (:dist %)) -1e-12) faces)
                        kept  (remove (set seen) faces)
                        edges (frequencies (mapcat (fn [{[a b c] :idx}]
                                                     [#{a b} #{b c} #{c a}]) seen))
                        horizon (keep (fn [[e c]] (when (= 1 c) (vec e))) edges)
                        new-f (keep (fn [[a b]] (face-of pts' [a b n] interior)) horizon)]
                    (recur pts' (into (vec kept) new-f) (inc i))))))))))))
