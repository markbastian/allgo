(ns procedurals.fps
  "A frame-rate readout pinned to a demo viewport.

  Reports two independent numbers: the rate requestAnimationFrame actually
  delivers, and the mean time the demo spends building a frame. They are not
  reciprocals -- at a comfortable 60fps the period is 16.7ms however little
  work you do, so the second number is the one carrying headroom.")

(def ^:private window-ms 500)

(defn- render! [^js el frames elapsed work]
  (set! (.-textContent el)
        (str (js/Math.round (/ (* 1000 frames) elapsed))
             " fps · " (.toFixed (/ work frames) 1) " ms")))

(defn meter!
  "Append a readout to `container`, returning `(fn [work-ms])` to call once
  per rendered frame with how long that frame's work took."
  [^js container]
  (let [el    (js/document.createElement "div")
        state (atom nil)]
    (set! (.-className el) "fps-meter")
    (.appendChild container el)
    (fn [work-ms]
      (let [now (js/performance.now)
            {:keys [frames since work last]} @state]
        (if (or (nil? last) (> (- now last) window-ms))
          ;; First frame, or resuming after the loop was stopped or the tab
          ;; was backgrounded. Either way the old window is meaningless.
          (reset! state {:frames 0 :since now :work 0.0 :last now})
          (let [frames  (inc frames)
                work    (+ work work-ms)
                elapsed (- now since)]
            (if (>= elapsed window-ms)
              (do (render! el frames elapsed work)
                  (reset! state {:frames 0 :since now :work 0.0 :last now}))
              (reset! state {:frames frames :since since :work work :last now}))))))))
