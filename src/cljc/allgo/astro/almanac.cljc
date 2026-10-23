(ns allgo.astro.almanac
  "The next or previous of each recurring event Meeus's chapters predict:
  lunar phases and eclipses, the seasons, the Moon's perigee and apogee,
  the Earth's perihelion, the planets' oppositions, conjunctions,
  greatest elongations and stations, and the Galilean satellites' shadows crossing
  Jupiter.

  Each of those chapters answers a different question -- the event
  nearest a given year -- so this turns them into the one a calendar
  asks: what comes after now? The event nearest the present is either
  the next one or the last, and stepping by its mean period finds the
  other. Eclipses are the exception, since most lunations have none: they
  are sought lunation by lunation until one turns up.

  Times are MJD (TT)."
  (:require [allgo.astro.eclipse :as eclipse]
            [allgo.astro.jupiter-moons :as jm]
            [allgo.astro.lunar-events :as lunar]
            [allgo.astro.phenomena :as phenomena]
            [allgo.astro.planet-orbits :as orbits]
            [allgo.astro.solar :as solar]
            [allgo.astro.time :as time]
            [clojure.math :as math]))

(def ^:private synodic-month 29.530588861)

(defn- searched [planet kind]
  (fn [mjd direction] (phenomena/search planet kind mjd direction)))

(def events
  "Each event: `:name`, `:period` in days (roughly -- only used to step
  from one to the next), and `:at`, which gives the MJD of the event
  nearest a decimal year, or nil if there is none that lunation. An
  event found by searching from a moment instead has `:find`, which takes
  the MJD and a direction, 1 or -1, and gives the MJD of the next or
  previous one."
  (array-map
   :new-moon      {:name "New Moon" :period synodic-month :at #(lunar/moon-phase % 0.0)}
   :first-quarter {:name "First quarter" :period synodic-month :at #(lunar/moon-phase % 0.25)}
   :full-moon     {:name "Full Moon" :period synodic-month :at #(lunar/moon-phase % 0.5)}
   :last-quarter  {:name "Last quarter" :period synodic-month :at #(lunar/moon-phase % 0.75)}
   :solar-eclipse {:name "Solar eclipse" :period synodic-month :sparse true
                   :at #(:mjd (eclipse/solar %))}
   :lunar-eclipse {:name "Lunar eclipse" :period synodic-month :sparse true
                   :at #(:mjd (eclipse/lunar %))}
   :perigee       {:name "Lunar perigee" :period 27.55455 :at #(first (lunar/perigee %))}
   :apogee        {:name "Lunar apogee" :period 27.55455 :at #(first (lunar/apogee %))}
   :march-equinox {:name "March equinox" :period 365.2422 :at #(solar/season-exact (math/round (double %)) :march)}
   :june-solstice {:name "June solstice" :period 365.2422 :at #(solar/season-exact (math/round (double %)) :june)}
   :september-equinox {:name "September equinox" :period 365.2422
                       :at #(solar/season-exact (math/round (double %)) :september)}
   :december-solstice {:name "December solstice" :period 365.2422
                       :at #(solar/season-exact (math/round (double %)) :december)}
   :perihelion    {:name "Earth at perihelion" :period 365.2596 :at #(orbits/perihelion :earth %)}
   :aphelion      {:name "Earth at aphelion" :period 365.2596 :at #(orbits/aphelion :earth %)}
   ;; The planets: from Meeus's series for the rows of his table 36.B that
   ;; were verified, and by searching VSOP87 for the rest (see
   ;; `phenomena/search`).
   :mercury-inferior {:name "Mercury at inferior conjunction" :period 115.8775
                      :at #(phenomena/phenomenon :mercury-inferior-conjunction %)}
   :mercury-superior {:name "Mercury at superior conjunction" :period 115.8775
                      :at #(phenomena/phenomenon :mercury-superior-conjunction %)}
   :mercury-east  {:name "Mercury at greatest eastern elongation" :period 115.8775
                   :at #(first (phenomena/phenomenon :mercury-greatest-east-elongation %))}
   :mercury-west  {:name "Mercury at greatest western elongation" :period 115.8775
                   :at #(first (phenomena/phenomenon :mercury-greatest-west-elongation %))}
   :mercury-first-station {:name "Mercury stationary, turning retrograde" :period 115.8775
                           :find (searched :mercury :first-station)}
   :mercury-second-station {:name "Mercury stationary, turning direct" :period 115.8775
                            :find (searched :mercury :second-station)}
   :venus-inferior {:name "Venus at inferior conjunction" :period 583.9214
                    :at #(phenomena/phenomenon :venus-inferior-conjunction %)}
   :venus-superior {:name "Venus at superior conjunction" :period 583.9214
                    :find (searched :venus :superior-conjunction)}
   :venus-east {:name "Venus at greatest eastern elongation" :period 583.9214
                :find (searched :venus :greatest-east-elongation)}
   :venus-west {:name "Venus at greatest western elongation" :period 583.9214
                :find (searched :venus :greatest-west-elongation)}
   :venus-first-station {:name "Venus stationary, turning retrograde" :period 583.9214
                         :find (searched :venus :first-station)}
   :venus-second-station {:name "Venus stationary, turning direct" :period 583.9214
                          :find (searched :venus :second-station)}
   :mars-opposition {:name "Mars at opposition" :period 779.9361
                     :at #(phenomena/phenomenon :mars-opposition %)}
   :mars-conjunction {:name "Mars in conjunction with the Sun" :period 779.9361
                      :find (searched :mars :conjunction)}
   :mars-first-station {:name "Mars stationary, turning retrograde" :period 779.9361
                        :find (searched :mars :first-station)}
   :mars-second-station {:name "Mars stationary, turning direct" :period 779.9361
                         :at #(phenomena/phenomenon :mars-second-station %)}
   :jupiter-opposition {:name "Jupiter at opposition" :period 398.8840
                        :at #(phenomena/phenomenon :jupiter-opposition %)}
   :jupiter-conjunction {:name "Jupiter in conjunction with the Sun" :period 398.8840
                         :find (searched :jupiter :conjunction)}
   :saturn-opposition {:name "Saturn at opposition" :period 378.0919
                       :at #(phenomena/phenomenon :saturn-opposition %)}
   :saturn-conjunction {:name "Saturn in conjunction with the Sun" :period 378.0919
                        :at #(phenomena/phenomenon :saturn-conjunction %)}
   :uranus-opposition {:name "Uranus at opposition" :period 369.6560
                       :at #(phenomena/phenomenon :uranus-opposition %)}
   :uranus-conjunction {:name "Uranus in conjunction with the Sun" :period 369.6560
                        :find (searched :uranus :conjunction)}
   :neptune-opposition {:name "Neptune at opposition" :period 367.4867
                        :at #(phenomena/phenomenon :neptune-opposition %)}
   :neptune-conjunction {:name "Neptune in conjunction with the Sun" :period 367.4867
                         :find (searched :neptune :conjunction)}
   ;; The middle of each, geometric -- see `jupiter-moons/shadow-transit`.
   :io-shadow {:name "Io's shadow on Jupiter" :period (jm/synodic-periods :io) :find #(:mid (jm/shadow-transit :io %1 %2))}
   :europa-shadow {:name "Europa's shadow on Jupiter" :period (jm/synodic-periods :europa)
                   :find #(:mid (jm/shadow-transit :europa %1 %2))}
   :ganymede-shadow {:name "Ganymede's shadow on Jupiter" :period (jm/synodic-periods :ganymede)
                     :find #(:mid (jm/shadow-transit :ganymede %1 %2))}
   :callisto-shadow {:name "Callisto's shadow on Jupiter" :period (jm/synodic-periods :callisto)
                     :find #(:mid (jm/shadow-transit :callisto %1 %2))}))

(defn- search [event mjd direction]
  (let [{:keys [at period sparse]} (events event)
        year (fn [t] (time/mjd->julian-epoch t))
        later? (if (pos? direction) #(> % (+ mjd 1e-4)) #(< % (- mjd 1e-4)))
        dy (* direction (/ period 365.25))]
    (loop [y (year mjd) i 0]
      (when (< i (if sparse 400 8))
        (let [t (at y)]
          (if (and t (later? t))
            t
            (recur (+ y dy) (inc i))))))))

(defn- step [event mjd direction]
  (if-let [find (:find (events event))]
    (find mjd direction)
    (search event mjd direction)))

(defn next-event
  "MJD of the first `event` (a key of `events`) after `mjd`."
  [event mjd]
  (step event mjd 1))

(defn previous-event
  "MJD of the last `event` before `mjd`."
  [event mjd]
  (step event mjd -1))
