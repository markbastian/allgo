(ns allgo.demo.clock
  "Simulation time that keeps to the wall clock, whatever the display's
  refresh rate.

  A demo that steps its simulation a fixed sixtieth of a second each time
  requestAnimationFrame calls it runs at real speed only on a 60 Hz
  display: twice as fast on a 120 Hz one, and slower whenever a frame is
  dropped. This counts out fixed steps against the time that has actually
  passed instead, carrying the remainder into the next frame, so the
  simulation keeps its timestep -- and so its behavior -- and only how
  many steps a frame takes changes.")

(defn fixed-steps
  "Returns `(fn [])`, to call once a rendered frame, giving how many steps
  of `dt` seconds are due.

  At most `max-steps` a frame: a simulation slower than real time would
  otherwise fall further behind each frame trying to catch up. A gap of
  more than a quarter second -- the loop stopped, or the tab hidden -- is
  dropped rather than made up, and the first frame after it takes one
  step."
  ([dt] (fixed-steps dt 4))
  ([dt max-steps]
   (let [state (atom {:last nil :owed 0.0})]
     (fn []
       (let [now (/ (js/performance.now) 1000.0)
             {:keys [last owed]} @state
             gap (if last (- now last) ##Inf)
             owed (if (> gap 0.25) dt (+ owed gap))
             n (min max-steps (js/Math.floor (/ owed dt)))]
         (reset! state {:last now
                        ;; What could not be stepped within the cap is
                        ;; forgotten, not owed: falling behind once should
                        ;; not mean hurrying for seconds after.
                        :owed (if (= n max-steps) 0.0 (- owed (* n dt)))})
         n)))))
