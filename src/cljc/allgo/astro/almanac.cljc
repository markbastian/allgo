(ns allgo.astro.almanac
  "The next or previous of each recurring event Meeus's chapters predict:
  lunar phases and eclipses, the seasons, the Moon's perigee and apogee,
  the Earth's perihelion, and the planets' oppositions, conjunctions and
  greatest elongations.

  Each of those chapters answers a different question -- the event
  nearest a given year -- so this turns them into the one a calendar
  asks: what comes after now? The event nearest the present is either
  the next one or the last, and stepping by its mean period finds the
  other. Eclipses are the exception, since most lunations have none: they
  are sought lunation by lunation until one turns up.

  Times are MJD (TT)."
  (:require [allgo.astro.eclipse :as eclipse]
            [allgo.astro.lunar-events :as lunar]
            [allgo.astro.phenomena :as phenomena]
            [allgo.astro.planet-orbits :as orbits]
            [allgo.astro.solar :as solar]
            [allgo.astro.time :as time]
            [clojure.math :as math]))

(def ^:private synodic-month 29.530588861)

(def events
  "Each event: `:name`, `:period` in days (roughly -- only used to step
  from one to the next), and `:at`, which gives the MJD of the event
  nearest a decimal year, or nil if there is none that lunation."
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
   :mercury-east  {:name "Mercury at greatest eastern elongation" :period 115.8775
                   :at #(first (phenomena/phenomenon :mercury-greatest-east-elongation %))}
   :mercury-west  {:name "Mercury at greatest western elongation" :period 115.8775
                   :at #(first (phenomena/phenomenon :mercury-greatest-west-elongation %))}
   :venus-inferior {:name "Venus at inferior conjunction" :period 583.9214
                    :at #(phenomena/phenomenon :venus-inferior-conjunction %)}
   :mars-opposition {:name "Mars at opposition" :period 779.9361
                     :at #(phenomena/phenomenon :mars-opposition %)}
   :jupiter-opposition {:name "Jupiter at opposition" :period 398.8840
                        :at #(phenomena/phenomenon :jupiter-opposition %)}
   :saturn-opposition {:name "Saturn at opposition" :period 378.0919
                       :at #(phenomena/phenomenon :saturn-opposition %)}
   :uranus-opposition {:name "Uranus at opposition" :period 369.6560
                       :at #(phenomena/phenomenon :uranus-opposition %)}
   :neptune-opposition {:name "Neptune at opposition" :period 367.4867
                        :at #(phenomena/phenomenon :neptune-opposition %)}))

(defn- step [event mjd direction]
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

(defn next-event
  "MJD of the first `event` (a key of `events`) after `mjd`."
  [event mjd]
  (step event mjd 1))

(defn previous-event
  "MJD of the last `event` before `mjd`."
  [event mjd]
  (step event mjd -1))
