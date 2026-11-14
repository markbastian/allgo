(ns allgo.probability.seismic
  "Seismic monitoring as an open-universe model: the structure of NET-VISA
  (Arora et al., 2013), as Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, section 15.2.3, describes it, much
  simplified.

  The Earth here is a ring -- a location is an angle -- with stations at
  fixed places around it. In a window of time there may be events, each
  with a location, a time and a magnitude, and how many is unknown. At
  each station an event produces one detection or none: the chance
  rising with magnitude and falling with distance, the detection's time
  the event's plus the travel time and some noise. Every station also
  reports false alarms, a Poisson process uniform over the window -- most
  detections, as at the real network. The evidence is the detections,
  their stations and times, with nothing to say which event made which;
  the query is the bulletin, the events that most probably happened.

  For one event the posterior is exact, on a grid of locations, times and
  magnitudes: against no event at all, each candidate's odds are its
  prior times, at every station, (1 - p) + p times the sum over that
  station's detections of the arrival's density over the false alarms'
  -- the event made none of them, or exactly one. `bulletin` builds a
  bulletin from that, one event at a time: the most probable, if it is
  more probable than none, its detections then taken away. It is greedy
  where NET-VISA searches, but it shares what the section singles out:
  detections are explained by events or put down to noise, and the more
  of them an event explains the better it is placed."
  (:require [allgo.numerics.special :as special]
            [allgo.search :as search]
            [clojure.math :as math]))

(defn ring-distance
  "The angle between two places on the ring."
  [a b]
  (let [d (mod (abs (- a b)) (* 2.0 math/PI))]
    (min d (- (* 2.0 math/PI) d))))

(defn model
  "A monitoring network and its physics:

    :stations    their places, angles
    :window      the length of the window of time
    :speed       how far (radians) a wave travels in a unit of time
    :timing      the standard deviation of an arrival's time
    :false-rate  each station's false alarms, expected per window
    :event-rate  events expected per window
    :detect      [a b c]: P(detection) = 1/(1 + e^-(a + b M - c d)) for
                 magnitude M at distance d
    :magnitudes  the magnitudes on the grid, with prior weights falling
                 tenfold a unit (Gutenberg and Richter)"
  [stations & {:as opts}]
  (merge {:stations (vec stations) :window 100.0 :speed 0.2 :timing 0.25
          :false-rate 3.0 :event-rate 0.5 :detect [-6.0 2.0 2.0]
          :magnitudes [3.0 3.5 4.0 4.5 5.0 5.5 6.0]}
         opts))

(defn p-detect
  "The probability that `station` detects an event of magnitude `mag` at
  `place`."
  [{[a b c] :detect} station place mag]
  (/ 1.0 (+ 1.0 (math/exp (- (- (+ a (* b mag)) (* c (ring-distance station place))))))))

(defn arrival
  "When an event at `place` and `time` arrives at `station`."
  [{:keys [speed]} station place time]
  (+ time (/ (ring-distance station place) speed)))

(defn- magnitude-prior [{:keys [magnitudes]}]
  (let [w (map #(math/pow 10.0 (- (first magnitudes) %)) magnitudes)
        z (reduce + w)]
    (zipmap magnitudes (map #(/ % z) w))))

(defn- log-odds
  "log of an event's posterior odds against none, its prior aside: at
  each station, the event made none of its detections, or one."
  [{:keys [stations window timing false-rate] :as m} detections place time mag]
  (let [noise (/ false-rate window)]
    (reduce + (map-indexed
               (fn [s station]
                 (let [p (p-detect m station place mag)
                       due (arrival m station place time)
                       ratio (reduce + (map (fn [t] (/ (/ (special/normal-pdf (/ (- t due) timing)) timing) noise))
                                            (get detections s)))]
                   (math/log (+ (- 1.0 p) (* p ratio)))))
               stations))))

(defn grid
  "The candidate events: `places` angles round the ring, `times` across
  the window, each magnitude."
  [{:keys [window magnitudes]} places times]
  (for [i (range places) j (range times) mag magnitudes]
    {:place (* 2.0 math/PI (/ i places)) :time (* window (/ j times)) :mag mag}))

(defn posterior
  "The posterior over one event or none given `detections` (a map from
  station index to its detections' times), on the candidates
  `candidates`: `{:none p :events [[candidate p] ...]}`."
  [{:keys [event-rate] :as m} detections candidates]
  (let [prior-mag (magnitude-prior m)
        per-place (/ 1.0 (count (distinct (map (juxt :place :time) candidates))))
        lw (mapv (fn [{:keys [place time mag] :as c}]
                   [c (+ (math/log event-rate) (math/log per-place) (math/log (prior-mag mag))
                         (log-odds m detections place time mag))])
                 candidates)
        top (reduce max 0.0 (map second lw))
        none (math/exp (- top))
        total (+ none (reduce + (map #(math/exp (- (second %) top)) lw)))]
    {:none (/ none total)
     :events (mapv (fn [[c w]] [c (/ (math/exp (- w top)) total)]) lw)}))

(defn explained
  "The detections the event `e` most probably made: at each station, the
  one whose arrival is likeliest to be the event's, if its being the
  event's -- p times the arrival's density -- is likelier than the event
  going unseen there and the detection being a false alarm."
  [{:keys [stations window timing false-rate] :as m} detections {:keys [place time mag]}]
  (into {} (keep (fn [s]
                   (let [station (stations s)
                         p (p-detect m station place mag)
                         due (arrival m station place time)
                         score (fn [t] (* p (/ (special/normal-pdf (/ (- t due) timing)) timing)))
                         best (when (seq (get detections s)) (apply max-key score (get detections s)))]
                     (when (and best (> (score best) (* (- 1.0 p) (/ false-rate window))))
                       [s best])))
                 (range (count stations)))))

(defn- bulletin-step
  "The next event of the bulletin and the detections left, or nil when no
  event is more probable than none, or the likeliest explains fewer than
  `min-detections`."
  [m candidates min-detections {:keys [detections events]}]
  (let [{:keys [none] es :events} (posterior m detections candidates)
        [best p] (apply max-key second es)
        made (explained m detections best)]
    (when (and (> (- 1.0 none) 0.5) (>= (count made) min-detections))
      {:detections (reduce (fn [ds [s t]] (update ds s (fn [ts] (vec (remove #{t} ts))))) detections made)
       :events (conj events (assoc best :p p :detections made))})))

(defn bulletin
  "The events that most probably happened, most probable first, each with
  the detections it explains: one event at a time, its detections taken
  away before the next is sought. As real bulletins do, an event must be
  defined by detections at `min-detections` stations (default 3): among
  a hundred false alarms a few always line up by chance."
  ([m detections candidates] (bulletin m detections candidates 3))
  ([m detections candidates min-detections]
   (:events (last (search/step-seq (partial bulletin-step m candidates min-detections)
                                   {:detections detections :events []})))))
