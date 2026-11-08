(ns allgo.robotics.reactive
  "Reactive control: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, section 26.9 -- a policy written
  directly, with no model of the world and no planning.

    the leg     an augmented finite state machine, the section's controller
                for one leg of a six-legged walker: lift the leg to its
                height h, swing it forward, set it down, push back -- and
                if the swing is blocked, retract, lift higher and swing
                again. h is the machine's own variable, standing for no
                measurement of the ground: the leg finds out by trying
    the gait    the tripod gait: the right front, right rear and left
                middle legs step while the other three hold the body up,
                then the other three, a group waiting for all its legs
                to be down before the other lifts -- so at least three
                legs are always on the ground
    subsumption Brooks's way of building up behaviors (\"A Robust Layered
                Control System for a Mobile Robot\", IEEE J. Robotics and
                Automation 2, 1986): layers of such machines, a higher
                one suppressing the outputs of those below when it has
                something to say -- here, backing away from what the
                bumper touches over walking on

  The walker steps along a strip of ground, each leg's next foothold an
  obstacle of some height, which a swing clears only if the leg is lifted
  higher.")

(def legs [:rf :rm :rr :lf :lm :lr])
(def tripods [#{:rf :rr :lm} #{:lf :lr :rm}])

(defn leg-step
  "One tick of a leg's machine, `{:state :h :retries}`, given the height
  of the `obstacle` in its swing's way and the heights `h0` and `dh`:
  :lift raises the leg to h; :swing clears the obstacle if h is above it,
  else the leg is :stuck; :stuck retracts it and raises h by `dh`, back
  to :lift; :down sets it down; :stance pushes back, and is where the
  machine rests until the gait tells it to lift again (`go?`). Setting
  down resets h."
  [{:keys [state h retries] :as leg} obstacle go? h0 dh]
  (case state
    :stance (if go? (assoc leg :state :lift) leg)
    :lift (assoc leg :state :swing)
    :swing (if (> h obstacle) (assoc leg :state :down) (assoc leg :state :stuck))
    :stuck (assoc leg :state :lift :h (+ h dh) :retries (inc retries))
    :down (assoc leg :state :stance :h h0)))

(defn grounded?
  "Whether the leg is bearing weight."
  [{:keys [state]}]
  (= state :stance))

(defn walker
  "A walker with every leg standing, lifting to `h0` first."
  [h0]
  {:legs (zipmap legs (repeat {:state :stance :h h0 :retries 0})) :group 0 :done #{} :strides 0 :ticks 0})

(defn walk-seq
  "The walker under the tripod gait, a tick a step: the current tripod's
  legs run their machines -- each meeting the obstacle `(terrain leg
  stride)` -- while the other three stand; a leg that has stepped waits,
  standing, for the rest of its tripod -- the machines' synchronizing
  message -- and when all three are down again the body moves a stride
  and the other tripod goes. Each state `{:legs :group :done :strides
  :ticks}`."
  [w terrain h0 dh]
  (iterate (fn [{:keys [legs group done strides ticks] :as w}]
             (let [moving (tripods group)
                   stepped (into {} (map (fn [[k leg]]
                                           [k (if (and (moving k) (not (done k)))
                                                (leg-step leg (terrain k strides) true h0 dh)
                                                leg)])
                                         legs))
                   done (into done (filter #(= :down (:state (legs %))) moving))]
               (if (= done moving)
                 {:legs stepped :group (- 1 group) :done #{} :strides (inc strides) :ticks (inc ticks)}
                 (assoc w :legs stepped :done done :ticks (inc ticks)))))
           w))

(defn subsumption
  "The command of the highest of `layers` -- each (fn [sensors]) ->
  command or nil, lowest first -- that has one: a layer above suppresses
  every layer below whenever it speaks."
  [layers sensors]
  (some #(% sensors) (reverse layers)))

(defn walk-layer
  "The bottom layer: always walk forward."
  [_]
  :walk)

(defn avoid-layer
  "A layer over walking: when the bumper touches something, back up and
  turn away from it."
  [{:keys [bumper]}]
  (case bumper :left :back-up-turn-right :right :back-up-turn-left nil))
