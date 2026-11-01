(ns allgo.logic.events
  "Events and time: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, chapter 10 (chapter 12 of the 3rd edition, whose
  text this follows).

    event calculus   fluents and events reified: `Happens(e, (t1, t2))`,
                     `Initiates(e, f, t)`, `Terminates(e, f, t)`, and
                     `T(f, t)` -- whether fluent f holds at time t -- by the
                     chapter's two axioms: f holds after an event initiates
                     it unless an event in between clipped it, and fails to
                     hold after one terminates it unless one restored it.
                     A fluent nothing has initiated or terminated is
                     neither: the axioms say nothing about it, and neither
                     does `holds` (nil). `timeline-seq` runs the axioms
                     forward, a time point a step
    processes        liquid events: a process that happens over an interval
                     happens over every subinterval
    intervals        the chapter's interval relations (after Allen):
                     Meet, Before, After, During, Overlap, Begins, Finishes
                     and Equals, with its definitions -- Begins and
                     Finishes asking only that the beginnings, or the ends,
                     agree -- and Overlap one way round, i before j
    time scale       seconds from midnight at the start of 1900, and the
                     chapter's six-argument `Date`

  A narrative is `{:events [{:id :category :interval [t1 t2]} ...]
  :initiates f :terminates f}`, the two functions taking an event and a
  time and returning the fluents it initiates or terminates then; times
  are integers. Fluents are anything comparable by `=`, `[:At :Shankar
  :SF]` say.")

;; ---------------------------------------------------------------------------
;; The time scale

(defn leap-year? [y] (and (zero? (mod y 4)) (or (pos? (mod y 100)) (zero? (mod y 400)))))

(defn- days-before-year
  "Days from 1 January of year 1 to 1 January of year `y`, Gregorian."
  [y]
  (let [y (dec (long y))] (+ (* 365 y) (quot y 4) (- (quot y 100)) (quot y 400))))

(def ^:private days-before-month [0 31 59 90 120 151 181 212 243 273 304 334])

(defn- day-number [day month year]
  (+ (days-before-year year) (long (days-before-month (dec (long month))))
     (if (and (> (long month) 2) (leap-year? year)) 1 0) (dec (long day))))

(defn date
  "The chapter's `Date`: the time, in seconds from midnight at the start of
  1900, of `seconds`, `minutes` and `hours` past midnight on `day` of
  `month` of `year`, read on a clock `utc-offset` hours from Greenwich
  (default 0).

  The text gives the arguments as hours first, but its example --
  Date(0, 20, 21, 24, 1, 1995) = Seconds(3000000000) -- holds only with
  seconds first, and on Pacific Standard Time, eight hours behind
  Greenwich, as do its other numbers; see NOTES."
  [seconds minutes hours day month year & {:keys [utc-offset] :or {utc-offset 0}}]
  (- (+ (* 86400 (- (long (day-number day month year)) (long (day-number 1 1 1900))))
        (* 3600 (long hours)) (* 60 (long minutes)) (long seconds))
     (long (* 3600 utc-offset))))

(defn year-interval
  "The interval of a year, from its first moment to the next year's:
  `AD2001`."
  [year & opts]
  [(apply date 0 0 0 1 1 year opts) (apply date 0 0 0 1 1 (inc (long year)) opts)])

(defn duration [[b e]] (- (long e) (long b)))

;; ---------------------------------------------------------------------------
;; Intervals

(defn meet? [[_ ei] [bj _]] (= ei bj))
(defn before? [[_ ei] [bj _]] (< (long ei) (long bj)))
(defn after? [j i] (before? i j))
(defn during? [[bi ei] [bj ej]] (< (long bj) (long bi) (long ei) (long ej)))
(defn overlap? [[bi ei] [bj ej]] (< (long bi) (long bj) (long ei) (long ej)))
(defn begins? [[bi _] [bj _]] (= bi bj))
(defn finishes? [[_ ei] [_ ej]] (= ei ej))
(defn equals? [i j] (and (begins? i j) (finishes? i j)))

(def relations
  "The chapter's interval relations, by name."
  {:meet meet? :before before? :after after? :during during? :overlap overlap?
   :begins begins? :finishes finishes? :equals equals?})

(defn relations-of
  "The names of the relations that hold between intervals `i` and `j`."
  [i j]
  (into #{} (keep (fn [[k r]] (when (r i j) k))) relations))

;; ---------------------------------------------------------------------------
;; Event calculus

(defn happens?
  "Happens(e, i): whether event `e` takes place over interval `i` --
  exactly, or, for an event of one of the `processes` categories, over
  any part of where it takes place."
  ([e i] (happens? e i #{}))
  ([{:keys [interval category]} [t2 t3] processes]
   (let [[t1 t4] interval]
     (or (= interval [t2 t3])
         (and (contains? processes category) (<= (long t1) (long t2)) (< (long t2) (long t3)) (<= (long t3) (long t4)))))))

(defn- touching
  "The events of `narrative` that start in [t1, t2) and `relation`
  (`:initiates` or `:terminates`) fluent `f` when they do."
  [narrative relation f t1 t2]
  (filter (fn [{[t] :interval :as e}]
            (and (<= (double t1) (double t)) (< (double t) (double t2))
                 (contains? (set ((relation narrative) e t)) f)))
          (:events narrative)))

(defn clipped?
  "Clipped(f, (t1, t2)): whether some event starting in [t1, t2)
  terminates `f`."
  [narrative f [t1 t2]]
  (boolean (seq (touching narrative :terminates f t1 t2))))

(defn restored?
  "Restored(f, (t1, t2)): whether some event starting in [t1, t2)
  initiates `f`."
  [narrative f [t1 t2]]
  (boolean (seq (touching narrative :initiates f t1 t2))))

(defn holds
  "T(f, t) by the chapter's axioms, read straight off them: true if an
  event initiated `f` at some t1 < t and it was not clipped in [t1, t);
  false if one terminated it and it was not restored; nil if neither
  follows -- nothing has said, or an initiation and a termination at the
  same moment each undo the other."
  [narrative f t]
  (let [starts (fn [rel] (map (comp first :interval) (touching narrative rel f ##-Inf t)))
        yes (some #(not (clipped? narrative f [% t])) (starts :initiates))
        no (some #(not (restored? narrative f [% t])) (starts :terminates))]
    (cond yes true no false :else nil)))

(defn holds-over
  "T(f, (t1, t2)): whether `f` holds at every time point of [t1, t2)."
  [narrative f [t1 t2]]
  (every? #(true? (holds narrative f %)) (range t1 t2)))

(defn- timeline-step
  [{:keys [narrative t last-on last-off]}]
  ;; The events starting at t decide the fluents from t + 1: each
  ;; fluent's latest initiation and latest termination so far.
  (let [starting (filter #(= t (first (:interval %))) (:events narrative))
        mark (fn [m rel] (reduce (fn [m e] (reduce #(assoc %1 %2 t) m ((rel narrative) e t))) m starting))
        last-on (mark last-on :initiates)
        last-off (mark last-off :terminates)
        t (inc (long t))
        status (into {} (for [f (distinct (concat (keys last-on) (keys last-off)))
                              :let [on (get last-on f ##-Inf) off (get last-off f ##-Inf)]]
                          [f (cond (> (double on) (double off)) true (> (double off) (double on)) false :else nil)]))]
    {:narrative narrative :t t :last-on last-on :last-off last-off :events starting
     :holds (into #{} (keep (fn [[f v]] (when (true? v) f))) status)
     :not-holds (into #{} (keep (fn [[f v]] (when (false? v) f))) status)
     :undecided (into #{} (keep (fn [[f v]] (when (nil? v) f))) status)}))

(defn timeline-seq
  "T at every time point from `t0` on, a step a time point: each state's
  `:holds`, `:not-holds` and `:undecided` the fluents the axioms make
  true, false and neither at `:t`, and `:events` the events that started
  at the time before. A fluent's status is that of its latest initiation
  or termination before t -- the same answer as `holds`, made a step at
  a time rather than by searching the whole past."
  [narrative t0]
  (iterate timeline-step {:narrative narrative :t t0 :last-on {} :last-off {} :events []
                          :holds #{} :not-holds #{} :undecided #{}}))
