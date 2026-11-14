(ns allgo.robotics.cylindrical
  "Cylindrical algebraic decomposition of the plane, Collins's: LaValle,
  *Planning Algorithms*, Cambridge University Press, 2006, section 6.4.2.

  A set of polynomials in x and y cuts the plane into cells on each of
  which every polynomial keeps its sign, so that any set the polynomials
  describe -- a semi-algebraic set, the free space of a planning problem
  -- is a union of cells, and a question about it is answered by one
  sample point per cell.

    1. Project. The places in x where the cells must change are where a
       curve folds over in y (its discriminant in y, the resultant of f
       and df/dy, vanishes), where two curves cross (their resultant),
       and where a curve runs off to infinity (its leading coefficient in
       y vanishes): Figure 6.33. These are polynomials in x alone; their
       real roots, found exactly by `allgo.numerics.algebraic`, cut the
       line into 0-cells and the open 1-cells between them (6.20). A root
       where the shared root or double root is complex changes nothing
       in the real plane and is skipped.
    2. Lift. Over each cell of the line stands a cylinder, which the
       curves' real roots in y over it -- sections -- cut into sectors
       (6.24-6.27). Between two critical values the sections neither
       cross nor fold, so they are found once, at a sample x, and over a
       critical value they are found there.
    3. Sample. On the line, each 1-cell's midpoint, or one past its end
       when unbounded, and each 0-cell's root; in each cylinder, the same
       along y.

  For planning, each free 2-cell is a vertex of a roadmap, joined to the
  free vertical 1-cells over the critical values beside it -- the
  technical detail the book notes, which cells of neighboring cylinders
  meet, settled by following each section to its limit -- and a path
  between two cells keeps to them, sliding across each from one side to
  the other.

  The polynomials must have integer coefficients, as the book's are in
  Q, and be square-free. The projection is exact. The lifting over a
  critical value is numerical: its roots in y are found in floating
  point, and two closer than a millionth are one -- the precision
  parameter the book says a practical implementation chooses."
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.graph :as graph]
            [allgo.numerics.algebraic :as alg]
            [allgo.numerics.multivariate :as mv]
            [allgo.numerics.polynomial :as poly]
            [allgo.search :as search]
            [clojure.math :as math]))

(def ^:private tol 1e-6)

(defn- close?
  ([a b] (close? tol a b))
  ([t a b] (<= (abs (- a b)) (* t (max 1.0 (abs a) (abs b))))))

(defn- sign [x] (cond (pos? x) 1 (neg? x) -1 :else 0))

;; ---------------------------------------------------------------- the line

(defn samples
  "Sample points for the cells the ascending `roots` cut the line into,
  in order: one less than the first, each root, the midpoint between two,
  one more than the last (6.20); 0 when there are no roots."
  [roots]
  (if (empty? roots)
    [0.0]
    (vec (concat [(dec (first roots))]
                 (mapcat (fn [[a b]] [a (* 0.5 (+ a b))]) (partition 2 1 roots))
                 [(peek roots) (inc (peek roots))]))))

(defn- distinct-sorted [xs]
  (reduce (fn [acc x] (if (and (seq acc) (<= (- x (peek acc)) (* 1e-12 (max 1.0 (abs x))))) acc (conj acc x)))
          [] (sort xs)))

(defn line-cells
  "The decomposition of the line by polynomials in one variable,
  `fs` (integer coefficients, highest power first): alternating 1-cells
  and 0-cells, each with its dimension, sample point and the sign of each
  polynomial there."
  [fs]
  (let [roots-of (mapv alg/real-roots fs)
        roots (distinct-sorted (apply concat roots-of))]
    (vec (map-indexed (fn [k x]
                        {:dim (if (odd? k) 0 1)
                         :sample x
                         :signs (mapv (fn [f rs] (if (and (odd? k) (some #(= x %) rs)) 0 (sign (poly/horner f x))))
                                      fs roots-of)})
                      (samples roots)))))

;; --------------------------------------------------- polynomials in y

(defn- in-y
  "`f`, a polynomial in x and y, as a vector of polynomials in x, the
  coefficients of y's powers highest first."
  [f]
  (when-not (every? #(== % (math/floor %)) (vals f))
    (throw (ex-info "Cylindrical decomposition takes integer coefficients" {:f f})))
  (let [dy (reduce max 0 (map #(% 1) (keys f)))
        dx (reduce max 0 (map #(% 0) (keys f)))]
    (vec (for [j (range dy -1 -1)]
           (alg/trim (vec (for [i (range dx -1 -1)] (double (get f [i j] 0.0)))))))))

(defn- y-derivative [fy]
  (let [d (dec (count fy))]
    (vec (map-indexed (fn [i c] (alg/scale c (double (- d i)))) (butlast fy)))))

(defn- at-x [fy x] (mapv #(poly/horner % x) fy))

(defn- size-at
  "How large f's coefficients in y could be at x, for telling one that
  vanishes from one that is merely small."
  [fy x]
  (reduce max 0.0 (map #(poly/horner (mapv abs %) (abs x)) fy)))

(defn- y-roots
  "The real roots of the polynomial in y with coefficients `cs`,
  ascending; a double root shows as two close ones."
  [cs]
  (let [big (reduce max 0.0 (map abs cs))
        cs (vec (drop-while #(<= (abs %) (* 1e-12 big)) cs))
        n (dec (count cs))]
    (if (< n 1)
      []
      (->> (case n
             1 [[(/ (- (cs 1)) (cs 0)) 0.0]]
             2 (apply poly/quadratic cs)
             3 (apply poly/cubic cs)
             4 (apply poly/quartic cs)
             (poly/bairstow cs))
           (filter (fn [[re im]] (<= (abs im) (* tol (max 1.0 (abs re))))))
           (map first)
           sort
           vec))))

(defn- cluster
  "Sections from `[y i]` pairs sorted by y: roots closer than `t`
  relatively merged, each with the polynomials vanishing there."
  [t pts]
  (mapv (fn [group] {:y (/ (reduce + (map first group)) (count group)) :zero (set (map second group))})
        (reduce (fn [groups [y :as pt]]
                  (if (and (seq groups) (close? t (first (peek (peek groups))) y))
                    (conj (pop groups) (conj (peek groups) pt))
                    (conj groups [pt])))
                [] pts)))

(defn- column-at
  "The sections over `x`: the real roots in y of each polynomial not
  vanishing on the whole vertical line there, merged when closer than
  `t` relatively; and those that do, a vertical line of their zeros.
  Over a critical value roots meet and `t` is the tolerance; between
  critical values they never do, and `t` is 0 -- near a crossing two
  sections are as close as you like."
  [fys x t]
  (let [evals (map-indexed (fn [i fy] [i (at-x fy x) (size-at fy x)]) fys)
        lines (set (keep (fn [[i cs s]] (when (every? #(<= (abs %) (* 1e-9 s)) cs) i)) evals))]
    {:roots (cluster t (sort-by first (for [[i cs] evals :when (not (lines i)) y (y-roots cs)] [y i])))
     :lines lines}))

;; ------------------------------------------------------------ projection

(defn- nonzero [p what]
  (if (empty? p)
    (throw (ex-info (str "The polynomials are not square-free and coprime: a " (name (:kind what)) " vanishes identically") what))
    p))

(defn- project [fys]
  (let [n (count fys) deg #(dec (count (fys %)))]
    (vec (concat
          (for [i (range n)]
            {:kind :leading :of [i] :poly (alg/square-free (first (fys i)))})
          (for [i (range n) :when (>= (deg i) 2)]
            {:kind :fold :of [i]
             :poly (alg/square-free (nonzero (alg/resultant (fys i) (y-derivative (fys i))) {:kind :fold :of [i]}))})
          (for [i (range n) j (range (inc i) n) :when (and (pos? (deg i)) (pos? (deg j)))]
            {:kind :crossing :of [i j]
             :poly (alg/square-free (nonzero (alg/resultant (fys i) (fys j)) {:kind :crossing :of [i j]}))})))))

(defn projection
  "The projection of polynomials `fs` in x and y onto x, as polynomials
  in x, highest power first, each with its kind -- :leading coefficient
  in y, :fold (the discriminant in y) or :crossing (the resultant of a
  pair) -- and the polynomials it is `:of`."
  [fs]
  (project (mapv in-y fs)))

(defn- real-event?
  "Whether the projection polynomial `p`'s root `x` is a fold or a
  crossing in the real plane, not only among complex y."
  [fys {:keys [kind of]} x]
  (case kind
    :leading true
    :fold (boolean (some (fn [[a b]] (close? a b)) (partition 2 1 (y-roots (at-x (fys (first of)) x)))))
    :crossing (let [[i j] of
                    bs (y-roots (at-x (fys j) x))]
                (boolean (some (fn [a] (some #(close? a %) bs)) (y-roots (at-x (fys i) x)))))))

(defn- critical-values
  "The critical values on the x axis, ascending, each with why: the
  real roots of the projection that are real events."
  [fys proj]
  (->> (for [{:keys [poly] :as p} proj
             :when (pos? (alg/degree poly))
             x (alg/real-roots poly)
             :when (real-event? fys p x)]
         {:x x :why [(dissoc p :poly)]})
       (sort-by :x)
       (reduce (fn [acc {:keys [x why] :as c}]
                 (if (and (seq acc) (<= (- x (:x (peek acc))) (* 1e-9 (max 1.0 (abs x)))))
                   (conj (pop acc) (update (peek acc) :why into why))
                   (conj acc c)))
               [])))

;; ------------------------------------------------------------- lifting

(defn- cells-of
  "The cells of a cylinder, bottom to top: sectors at even positions,
  sections at odd ones, each with its sample and signs."
  [fs {:keys [index kind x roots lines]}]
  (vec (map-indexed (fn [p y]
                      (let [section? (odd? p)
                            zero (into lines (when section? (:zero (roots (quot p 2)))))]
                        {:id [index p]
                         :dim (- (if (= kind :interval) 2 1) (if section? 1 0))
                         :sample [x y]
                         :signs (mapv (fn [i f] (if (zero i) 0 (sign (mv/evaluate f [x y])))) (range) fs)}))
                    (samples (mapv :y roots)))))

(defn decompose
  "The cylindrical algebraic decomposition of the plane by `fs`,
  polynomials in x and y as `allgo.numerics.multivariate` builds them,
  with integer coefficients: `:projection`, `:critical` values in x with
  why, and `:columns`, the cylinders left to right -- over open intervals
  at even positions, over critical values at odd ones -- each with its
  `:x` sample, its `:lo` and `:hi` ends (nil when unbounded), its
  `:roots` -- sections, each with the polynomials vanishing there -- and
  its `:cells`."
  [fs]
  (let [fys (mapv in-y fs)
        proj (project fys)
        crit (critical-values fys proj)
        xs (mapv :x crit)
        k (count xs)]
    {:polys fs
     :y-polys fys
     :projection proj
     :critical crit
     :columns (vec (map-indexed
                    (fn [c x]
                      (let [i (quot c 2)
                            column (merge {:index c :x x}
                                          (if (even? c)
                                            {:kind :interval
                                             :lo (when (pos? i) (xs (dec i)))
                                             :hi (when (< i k) (xs i))}
                                            {:kind :point :lo x :hi x})
                                          (column-at fys x (if (even? c) 0.0 tol)))]
                        (assoc column :cells (cells-of fs column))))
                    (samples xs)))}))

(defn cells
  "Every cell of the decomposition."
  [cad]
  (mapcat :cells (:columns cad)))

(defn cell
  "The cell with id `[column position]`."
  [cad [c p]]
  (get-in cad [:columns c :cells p]))

(defn decide
  "The first cell whose signs satisfy `pred`, or nil: the decision
  problem, answered by trying one sample per cell."
  [cad pred]
  (first (filter (comp pred :signs) (cells cad))))

;; ---------------------------------------------------- neighboring cylinders

(defn- limits
  "Where the sections over interval column `c` end at its end `b`: for
  each, the index of its limit among the sections over b, or -1 or the
  count of them when it runs off to minus or plus infinity. Found by
  following each to just short of b and taking the nearest."
  [{:keys [columns y-polys]} c b]
  (let [{:keys [lo hi roots]} (columns c)
        width (if (and lo hi) (- hi lo) 1.0)
        near (+ b (* (if (= b hi) -1.0 1.0) 1e-9 (max 1.0 width)))
        ys (mapv :y (:roots (column-at y-polys near 0.0)))
        bs (mapv :y (:roots (columns (if (= b hi) (inc c) (dec c)))))
        far (* 1e6 (reduce max 1.0 (map abs bs)))]
    (when-not (= (count ys) (count roots))
      (throw (ex-info "Lost a section near a critical value" {:column c :at near :found ys :expected (count roots)})))
    (mapv (fn [y]
            (cond
              (< y (- far)) -1
              (> y far) (count bs)
              :else (apply min-key #(abs (- y (bs %))) (range (count bs)))))
          ys)))

(defn- touching
  "The positions of the vertical sectors over end `b` of interval column
  `c` that sector `s` of it meets: those between the limits of the
  sections below and above it. None when the two meet in a point."
  [cad c s b]
  (let [lim (limits cad c b)
        n (count (get-in cad [:columns (if (= b (get-in cad [:columns c :hi])) (inc c) (dec c)) :roots]))
        lower (if (zero? s) -1 (lim (dec s)))
        upper (if (= s (count lim)) n (lim s))]
    (range (inc lower) (inc upper))))

;; --------------------------------------------------------------- roadmap

(defn roadmap
  "The roadmap of the cells whose signs satisfy `free?`: each free 2-cell
  joined to the free vertical 1-cells beside it in the neighboring
  cylinders, and two sectors of one cylinder joined through a free
  section between them. `:graph` maps each cell id to its neighbors and
  the distances between their samples."
  [cad free?]
  (let [free (set (map :id (filter (comp free? :signs) (cells cad))))
        sample #(:sample (cell cad %))
        edge (fn [graph a b]
               (let [d (polygon/distance (sample a) (sample b))]
                 (-> graph (assoc-in [a b] d) (assoc-in [b a] d))))
        columns (:columns cad)
        across (for [{:keys [index kind lo hi roots]} columns
                     :when (= kind :interval)
                     s (range (inc (count roots)))
                     :when (free [index (* 2 s)])
                     b (remove nil? [lo hi])
                     t (touching cad index s b)
                     :let [other [(if (= b hi) (inc index) (dec index)) (* 2 t)]]
                     :when (free other)]
                 [[index (* 2 s)] other])
        up (for [{:keys [index roots]} columns
                 s (range (count roots))
                 :let [below [index (* 2 s)] through [index (inc (* 2 s))] above [index (+ 2 (* 2 s))]]
                 :when (and (free below) (free through) (free above))]
             [below above])]
    (assoc cad
           :free free
           :graph (reduce (fn [g [a b]] (edge g a b))
                          (into {} (map (fn [id] [id {}])) (filter (fn [[c p]] (and (even? c) (even? p))) free))
                          (concat across up)))))

(defn components
  "The connected components of the roadmap, as sets of cell ids."
  [{:keys [graph]}]
  (graph/components-by (sort (keys graph)) #(keys (graph %))))

(defn sections
  "The heights of the sections over `x`, ascending, an x between critical
  values: there the cylinder's curves, for drawing them."
  [{:keys [y-polys]} x]
  (mapv :y (:roots (column-at y-polys x 0.0))))

(defn locate
  "The id of the cell holding the point `[x y]`."
  [{:keys [columns y-polys]} [x y]]
  (let [{:keys [index kind roots]}
        (first (filter (fn [{:keys [kind lo hi]}]
                         (if (= kind :point)
                           (<= (abs (- x lo)) (* 1e-12 (max 1.0 (abs x))))
                           (and (or (nil? lo) (> x lo)) (or (nil? hi) (< x hi)))))
                       columns))
        ys (mapv :y (if (= kind :point) roots (:roots (column-at y-polys x 0.0))))
        on (first (filter #(close? y (ys %)) (range (count ys))))]
    [index (if on (inc (* 2 on)) (* 2 (count (filter #(< % y) ys))))]))

(defn- sector-bounds
  "Sector `s` of interval column `c` at `x`, as [below above], nil for an
  unbounded side; at one of the column's ends, the limits of its
  sections there."
  [{:keys [columns y-polys] :as cad} c s x]
  (let [{:keys [lo hi roots]} (columns c)
        n (count roots)]
    (if (or (= x lo) (= x hi))
      (let [lim (limits cad c x)
            bs (mapv :y (:roots (columns (if (= x hi) (inc c) (dec c)))))
            y-of (fn [k] (when (< -1 k (count bs)) (bs k)))]
        [(when (pos? s) (y-of (lim (dec s)))) (when (< s n) (y-of (lim s)))])
      (let [ys (mapv :y (:roots (column-at y-polys x 0.0)))]
        (when (= n (count ys))
          [(when (pos? s) (ys (dec s))) (when (< s n) (ys s))])))))

(defn- to-u
  "Height `y` in a sector with bounds [lo hi] as a coordinate that stays
  inside it as the bounds move: a fraction between two bounds, an offset
  from the one bound of an unbounded sector."
  [[lo hi] y]
  (cond (and lo hi) (/ (- y lo) (- hi lo)) lo (- y lo) hi (- y hi) :else y))

(defn- from-u [[lo hi] u]
  (cond (and lo hi) (+ lo (* u (- hi lo))) lo (+ lo u) hi (+ hi u) :else u))

(defn- across-sector
  "A path from `a` to `b` inside sector `s` of interval column `c`,
  either of them possibly on one of its ends: the sector coordinate
  carried linearly in x while the bounds follow the curves."
  [cad c s [xa ya :as a] [xb yb :as b]]
  (if (== xa xb)
    [a b]
    (let [ua (to-u (sector-bounds cad c s xa) ya)
          ub (to-u (sector-bounds cad c s xb) yb)
          n 32]
      (vec (concat [a]
                   (keep (fn [k]
                           (let [t (/ k n) x (+ xa (* t (- xb xa)))]
                             (when-let [bounds (sector-bounds cad c s x)]
                               [x (from-u bounds (+ ua (* t (- ub ua))))])))
                         (range 1 n))
                   [b])))))

(defn- leg
  "A path from `a` to `b` through the cell with id `[c p]`: across a
  2-cell, straight within anything else."
  [cad [c p] a b]
  (if (and (even? c) (even? p))
    (across-sector cad c (quot p 2) a b)
    [a b]))

(defn link-path
  "The way between the samples of neighboring cells `m` and `n`: up or
  down within a cylinder, or across the 2-cell of the two."
  [rm m n]
  (let [sample #(:sample (cell rm %))]
    (cond
      (= (first m) (first n)) [(sample m) (sample n)]
      (even? (first m)) (leg rm m (sample m) (sample n))
      :else (leg rm n (sample m) (sample n)))))

(defn route
  "A path from `a` to `b` through the free cells of roadmap `rm`, as
  points, or nil if they are not joined."
  [{:keys [graph] :as rm} a b]
  (let [from (locate rm a) to (locate rm b)
        sample #(:sample (cell rm %))]
    (when (and (graph from) (graph to))
      (when-let [nodes (if (= from to)
                         [from]
                         (some->> (search/dijkstra-seq {:start from :goal to
                                                        :neighbors #(keys (graph %))
                                                        :cost (fn [x y] (get-in graph [x y]))})
                                  search/goal-state search/path))]
        (vec (dedupe (concat (leg rm from a (sample from))
                             (mapcat (fn [[m n]] (link-path rm m n)) (partition 2 1 nodes))
                             (leg rm to (sample to) b))))))))
