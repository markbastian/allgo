(ns allgo.astro.calendar
  "The date of Easter, and the Julian, Gregorian, Jewish and Moslem
  calendars (Meeus, *Astronomical Algorithms*, chapters 7 to 9).

  All of these are arithmetic: a calendar is a rule for counting days,
  chosen so that its months and years stay near the Moon and the Sun, and
  once the rule is written down it can be computed without any astronomy
  at all. Easter is the clearest case -- its definition sounds like
  astronomy, the first Sunday after the first full moon of spring, but the
  full moon meant is an ecclesiastical one from a table, not the Moon in
  the sky, and the dates below are the tables' dates.

  Years are astronomical: 0 is 1 BC, -1 is 2 BC. Julian Day arithmetic
  lives in `allgo.astro.time`, which takes Julian dates before 1582
  October 15 and Gregorian from then on; the functions here convert a
  date between the two calendars when that switch is not what is wanted."
  (:require [allgo.astro.time :as time]
            [clojure.math :as math]))

(defn- fdiv [a b] (long (math/floor (/ (double a) b))))

;; ------------------------------------------------------------ Easter (8)

(defn easter
  "`[month day]` of Easter Sunday in the Gregorian calendar, for any year
  from 1583 -- Meeus's rendering of the algorithm published anonymously
  in *Nature* in 1876, which needs no exceptions and no table."
  [year]
  (let [a (mod year 19)
        b (quot year 100) c (mod year 100)
        d (quot b 4)      e (mod b 4)
        f (quot (+ b 8) 25)
        g (quot (+ (- b f) 1) 3)
        h (mod (+ (* 19 a) (- b d g) 15) 30)
        i (quot c 4)      k (mod c 4)
        l (mod (- (+ 32 (* 2 e) (* 2 i)) h k) 7)
        m (quot (+ a (* 11 h) (* 22 l)) 451)
        n (+ h l (* -7 m) 114)]
    [(quot n 31) (inc (mod n 31))]))

(defn julian-easter
  "`[month day]` of Easter in the Julian calendar, as the Orthodox churches
  still reckon it. The cycle repeats every 532 years."
  [year]
  (let [a (mod year 4)
        b (mod year 7)
        c (mod year 19)
        d (mod (+ (* 19 c) 15) 30)
        e (mod (+ (* 2 a) (* 4 b) (- d) 34) 7)
        f (+ d e 114)]
    [(quot f 31) (inc (mod f 31))]))

;; ------------------------------------------ Julian and Gregorian dates (7)

(defn- jdn->date
  "Calendar date from the integer part of JD + 0.5 plus 1524, the common
  tail of Meeus's conversions."
  [b]
  (let [c (fdiv (- (* 100 b) 12210) 36525)
        d (fdiv (* 36525 c) 100)
        e (fdiv (* 10000 (- b d)) 306001)
        day (- b d (fdiv (* 306001 e) 10000))
        month (if (< e 14) (dec e) (- e 13))]
    [(if (> month 2) (- c 4716) (- c 4715)) month day]))

(defn julian-date->mjd
  "MJD of a date in the Julian calendar, whatever the year -- unlike
  `time/calendar->mjd`, which switches to Gregorian at the reform."
  ([year month day] (julian-date->mjd year month day 0.0))
  ([year month day hour]
   (let [[y m] (if (<= month 2) [(dec year) (+ month 12)] [year month])]
     (+ (fdiv (* 36525 (+ y 4716)) 100) (fdiv (* 306 (inc m)) 10)
        day -1524.5 (/ hour 24.0) -2400000.5))))

(defn gregorian-date->mjd
  "MJD of a date in the Gregorian calendar, whatever the year -- the
  proleptic calendar before 1582."
  ([year month day] (gregorian-date->mjd year month day 0.0))
  ([year month day hour]
   (let [[y m] (if (<= month 2) [(dec year) (+ month 12)] [year month])
         a (fdiv y 100)]
     (+ (fdiv (* 36525 (+ y 4716)) 100) (fdiv (* 306 (inc m)) 10)
        (+ (- 2 a) (fdiv a 4))
        day -1524.5 (/ hour 24.0) -2400000.5))))

(defn mjd->julian-date
  "`[year month day]` in the Julian calendar, with a fractional day."
  [mjd]
  (let [jd (+ mjd 2400000.5 0.5)
        z  (long (math/floor jd))
        [y m d] (jdn->date (+ z 1524))]
    [y m (+ d (- jd z))]))

(defn mjd->gregorian-date
  "`[year month day]` in the Gregorian calendar, with a fractional day."
  [mjd]
  (let [jd (+ mjd 2400000.5 0.5)
        z  (long (math/floor jd))
        a  (fdiv (- (* 100 z) 186721625) 3652425)
        [y m d] (jdn->date (+ z 1 a (- (fdiv a 4)) 1524))]
    [y m (+ d (- jd z))]))

(defn gregorian->julian
  "The Julian-calendar date of a Gregorian one (Meeus example 9.c)."
  [year month day]
  (let [[y m d] (mjd->julian-date (gregorian-date->mjd year month day))]
    [y m (long d)]))

(defn julian->gregorian
  "The Gregorian-calendar date of a Julian one."
  [year month day]
  (let [[y m d] (mjd->gregorian-date (julian-date->mjd year month day))]
    [y m (long d)]))

;; ----------------------------------------------------- the Jewish year (9)

(defn- pesach-day
  "Day of March of Pesach -- Nisan 15 -- in the Julian calendar through
  1582 and the Gregorian after, from Gauss's formula as Meeus gives it.
  Past 31 it runs on into April."
  [year]
  (let [C (fdiv year 100)
        S (if (>= year 1583) (fdiv (- (* 3 C) 5) 4) 0)
        a (mod (+ (* 12 year) 12) 19)
        b (mod year 4)
        Q (+ -1.904412361576 (* 1.554241796621 a) (* 0.25 b)
             (* -0.003177794022 year) S)
        iq (long (math/floor Q))
        j (mod (+ iq (* 3 year) (* 5 b) 2 (- S)) 7)
        r (- Q iq)]
    (cond
      (#{2 4 6} j)                                   (+ iq 23)
      (and (= j 1) (> a 6) (>= r 0.63287037))        (+ iq 24)
      (and (= j 0) (> a 11) (>= r 0.897723765))      (+ iq 23)
      :else                                          (+ iq 22))))

(defn jewish-year
  "The Jewish year that begins in the autumn of Christian `year`, and when
  it falls:

    :year      the Jewish year (A.M.) beginning at :new-year, `year` + 3761
    :pesach    `[month day]` of Pesach in the spring of `year`, which
               belongs to the Jewish year before
    :new-year  `[month day]` of Rosh Hashanah, 163 days after Pesach
    :months    12, or 13 in the seven embolismic years of the 19-year cycle
    :days      the length of the year then beginning, 353 to 385

  Dates are Julian before 1583 and Gregorian after."
  [year]
  (let [A (+ year 3760)
        D (pesach-day year)
        D1 (pesach-day (inc year))]
    {:year     (inc A)
     :pesach   (if (> D 31) [4 (- D 31)] [3 D])
     :new-year (let [d (- D 21)] (if (> d 30) [10 (- d 30)] [9 d]))
     :months   (if (#{0 3 6 8 11 14 17} (mod (inc A) 19)) 13 12)
     :days     (+ (if (time/leap-year? (inc year)) 366 365) (- D1 D))}))

;; ---------------------------------------------------- the Moslem year (9)

(defn moslem-leap-year?
  "Eleven years in each thirty have 355 days rather than 354, which keeps
  the purely lunar year within a day of twelve lunations."
  [year]
  (> (mod (+ (* 11 (mod year 30)) 3) 30) 18))

(defn moslem->julian
  "The Julian-calendar `[year month day]` of a date in the (tabular)
  Moslem calendar."
  [year month day]
  (let [N  (+ day (fdiv (+ (* 295001 (dec month)) 9900) 10000))
        Q  (fdiv year 30)
        R  (mod year 30)
        A  (fdiv (+ (* 11 R) 3) 30)
        W  (+ (* 404 Q) (* 354 R) 208 A)
        Q1 (fdiv W 1461)
        Q2 (mod W 1461)
        G  (+ 621 (* 28 Q) (* 4 Q1))
        K  (fdiv (* Q2 10000) 3652422)
        E  (fdiv (* 3652422 K) 10000)
        J  (+ (- Q2 E) N -1)
        X  (+ G K)
        [X J] (cond
                (and (> J 366) (zero? (mod X 4))) [(inc X) (- J 366)]
                (and (> J 365) (pos? (mod X 4)))  [(inc X) (- J 365)]
                :else [X J])
        [m d] (time/day-of-year->date J (time/julian-leap-year? X))]
    [X m d]))

(defn moslem->gregorian
  "The Gregorian-calendar `[year month day]` of a Moslem date."
  [year month day]
  (apply julian->gregorian (moslem->julian year month day)))

(defn julian->moslem
  "The Moslem `[year month day]` of a date in the Julian calendar."
  [year month day]
  (let [W  (if (zero? (mod year 4)) 1 2)
        N  (+ (- (fdiv (* 275 month) 9) (* W (fdiv (+ month 9) 12))) day -30)
        A  (- year 623)
        B  (fdiv A 4)
        C  (let [c1 (* 365.25001 (mod A 4))
                 c2 (math/floor c1)]
             (long (if (> (- c1 c2) 0.5) (inc c2) c2)))
        D' (+ (* 1461 B) 170 C)
        Q  (fdiv D' 10631)
        R  (mod D' 10631)
        J  (fdiv R 354)
        K  (mod R 354)
        O  (fdiv (+ (* 11 J) 14) 30)
        H  (+ (* 30 Q) J 1)
        JJ (+ (- K O) N -1)
        days (if (moslem-leap-year? H) 355 354)
        [H JJ] (if (> JJ days) [(inc H) (- JJ days)] [H JJ])]
    (if (= JJ 355)
      [H 12 30]
      (let [S (fdiv (* 10 (dec JJ)) 295)]
        [H (inc S) (fdiv (- (* 10 JJ) (* 295 S)) 10)]))))

(defn gregorian->moslem
  "The Moslem `[year month day]` of a Gregorian date."
  [year month day]
  (apply julian->moslem (gregorian->julian year month day)))
