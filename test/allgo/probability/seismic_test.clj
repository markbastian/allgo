(ns allgo.probability.seismic-test
  "Seismic monitoring as an open-universe model, after the NET-VISA
  system that Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, section 15.2.3, describes: events found among
  detections nine in ten of which are false alarms, as at the real
  network. The section gives no numbers to check, so the detections are
  simulated from the model and the bulletin judged against the events
  that made them."
  (:require [allgo.probability.seismic :as seismic]
            [allgo.random :as random]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(def ^:private stations (mapv #(* 2.0 math/PI (/ % 12.0)) (range 12)))

;; nine false alarms a station a window, about nine in ten detections
(def ^:private m (seismic/model stations :false-rate 9.0 :event-rate 1.0))

;; three degrees by one time unit; the events below fall on it
(def ^:private candidates (vec (seismic/grid m 120 100)))

(defn- simulate
  "Each station's detections: its false alarms, and each event's arrival
  when the station detects it."
  [m rng events]
  (into {} (for [[s station] (map-indexed vector (:stations m))]
             [s (vec (concat (for [_ (range 100) :when (< (rng) (/ (:false-rate m) 100.0))] (* (rng) (:window m)))
                             (for [{:keys [place time mag]} events
                                   :when (< (rng) (seismic/p-detect m station place mag))]
                               (random/gaussian rng (seismic/arrival m station place time) (:timing m)))))])))

(defn- close-to [found event]
  (and (< (seismic/ring-distance (:place found) (:place event)) 0.12)
       (< (abs (- (:time found) (:time event))) 2.0)))

(deftest noise-alone
  (testing "from false alarms alone, nothing like an event: at most a weak one that a chance alignment of a few detections suggests, the price of a model that lets small events go unseen"
    (doseq [seed (range 1 6)]
      (doseq [e (seismic/bulletin m (simulate m (random/rng seed) []) candidates)]
        (is (< (:mag e) 4.5))
        (is (<= (count (:detections e)) 6))))))

(deftest events-among-false-alarms
  (testing "two events found, and where and when they were, though nine detections in ten are false"
    (doseq [seed (range 1 4)]
      (let [events [{:place 1.0 :time 30.0 :mag 5.0} {:place (* 0.7 2.0 math/PI) :time 60.0 :mag 5.0}]
            detections (simulate m (random/rng seed) events)
            found (seismic/bulletin m detections candidates)
            real (* 2 12)]
        (is (> (reduce + (map count (vals detections))) (* 4 real)) "mostly false alarms")
        (is (= 2 (count found)))
        (is (every? (fn [e] (some #(close-to % e) found)) events))
        (is (every? #(>= (count (:detections %)) 7) found) "each defined by more detections than any chance alignment")))))

(deftest one-event-posterior
  (let [event {:place 2.0 :time 40.0 :mag 5.0}
        detections (simulate m (random/rng 7) [event])
        {:keys [none events]} (seismic/posterior m detections candidates)]
    (testing "the posterior over one event or none sums to one, and puts nearly all of it on an event"
      (is (< (abs (- 1.0 (+ none (reduce + (map second events))))) 1e-9))
      (is (< none 1e-4)))
    (testing "near where and when it happened"
      (is (> (reduce + (for [[c p] events :when (close-to c event)] p)) 0.95)))))

(deftest more-stations-place-it-better
  (testing "the more detections an event explains, the more closely it is placed: the location's spread on a denser network is smaller"
    (let [event {:place 2.0 :time 40.0 :mag 4.5}
          spread (fn [n]
                   (let [mm (seismic/model (mapv #(* 2.0 math/PI (/ % n)) (range n)) :false-rate 9.0 :event-rate 1.0)
                         {:keys [events]} (seismic/posterior mm (simulate mm (random/rng 3) [event]) (seismic/grid mm 120 100))
                         z (reduce + (map second events))]
                     (math/sqrt (/ (reduce + (for [[c p] events] (* p (math/pow (seismic/ring-distance (:place c) 2.0) 2)))) z))))]
      (is (< (spread 24) (spread 6))))))
