(ns allgo.demo.progressive
  "A picture computed a pixel at a time, drawn coarse first and sharpened
  over as many frames as it takes, so the page never stalls.

  For the demos whose every pixel is a question asked of a function --
  textures, volumes -- and too slow to answer a whole canvas in a frame.
  Each pass samples one pixel in a block and paints the whole block with
  it, the blocks halving each pass, and a frame works only until its
  budget is spent.")

(defn job
  "A picture to make: `sample` is `(fn [px py] -> [r g b])` in [0, 1], over
  a `w` by `h` canvas, in passes of the block sizes in `blocks`."
  [^js ctx w h sample blocks]
  {:img (.createImageData ctx w h)
   :sample sample :w w :h h :blocks blocks
   :at (atom {:level 0 :row 0 :x 0})})

(defn done? [{:keys [at blocks]}]
  (>= (:level @at) (count blocks)))

(defn work!
  "Samples blocks until `budget-ms` is spent; true when the picture is
  finished. The budget is checked every block, not every row: a row of a
  slow picture can be most of a second."
  [{:keys [^js img sample w h blocks at]} budget-ms]
  (let [data (.-data img)
        t0 (js/performance.now)]
    (loop []
      (let [{:keys [level row x]} @at]
        (cond
          (>= level (count blocks)) true
          (> (- (js/performance.now) t0) budget-ms) false
          :else
          (let [b (nth blocks level)
                [r g bl] (sample (+ x (/ b 2)) (+ row (/ b 2)))
                r (* 255 r) g (* 255 g) bl (* 255 bl)]
            (dotimes [dy (min b (- h row))]
              (dotimes [dx (min b (- w x))]
                (let [p (* 4 (+ x dx (* w (+ row dy))))]
                  (aset data p r)
                  (aset data (+ p 1) g)
                  (aset data (+ p 2) bl)
                  (aset data (+ p 3) 255))))
            (reset! at (cond
                         (< (+ x b) w) {:level level :row row :x (+ x b)}
                         (< (+ row b) h) {:level level :row (+ row b) :x 0}
                         :else {:level (inc level) :row 0 :x 0}))
            (recur)))))))
