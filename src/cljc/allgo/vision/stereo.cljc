(ns allgo.vision.stereo
  "Binocular stereopsis: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 25.

  Two cameras side by side see a point at slightly different places in
  their images: the disparity, larger the nearer the point. With the
  cameras' image rows aligned, a point's match lies on the same row, so
  finding disparity is a one-dimensional search: for each pixel of the
  left image, the shift along the row whose window in the right image
  matches best by the sum of squared differences -- the chapter's
  approach, winner take all. Depth follows from `allgo.vision.camera`'s
  f b / d.

  A right camera to the right of the left one sees everything shifted
  left: the pixel at x in the left image is at x - d in the right."
  (:require [allgo.vision.camera :as camera]
            [allgo.vision.flow :as flow]
            [allgo.vision.image :as img]))

(defn disparity-at
  "The disparity at `(x, y)` of the left image: the shift from 0 to
  `max-d` whose window of half-width `half` in the right image matches
  best, smallest on a tie."
  [left right x y half max-d]
  (apply min-key #(flow/ssd left right x y half (- %) 0) (reverse (range (inc max-d)))))

(defn disparity-seq
  "The disparity map, a row at a time: each state the image of
  disparities so far, the rows not yet reached zero."
  [left right half max-d]
  (let [{:keys [w h]} left]
    (reductions (fn [d y]
                  (reduce (fn [d x] (assoc-in d [:px (+ (* y w) x)] (disparity-at left right x y half max-d)))
                          d (range w)))
                {:w w :h h :px (vec (repeat (* w h) 0))}
                (range h))))

(defn disparity-map
  "The disparity of every pixel of the left image."
  [left right half max-d]
  (last (disparity-seq left right half max-d)))

(defn depth-map
  "Depths from a disparity map, for cameras `baseline` apart with focal
  length `f`: nil where the disparity is zero -- as far as can be told,
  at infinity."
  [disparities f baseline]
  (img/pixels disparities #(when (pos? %) (camera/stereo-depth f baseline %))))

(defn stereo-pair
  "A synthetic stereo pair, `[left right]`, `w` by `h`, of `layers` back
  to front: each `{:texture (fn [x y]) :disparity d}`, and `:box [x0 y0
  x1 y1]` where it covers less than the whole picture -- the box and
  texture as the left camera sees them, the right seeing each layer
  shifted left by its disparity. The truth, `{[x y] d}` as the left camera
  sees it, is `(stereo-truth layers w h)`."
  [layers w h]
  (let [view (fn [shift?]
               (img/image w h (fn [x y]
                                (let [cx (+ x 0.5) cy (+ y 0.5)]
                                  (reduce (fn [v {:keys [texture disparity box]}]
                                            (let [sx (if shift? (+ cx disparity) cx)]
                                              (if (or (nil? box)
                                                      (let [[x0 y0 x1 y1] box] (and (<= x0 sx) (< sx x1) (<= y0 cy) (< cy y1))))
                                                (texture sx cy)
                                                v)))
                                          0.0 layers)))))]
    [(view false) (view true)]))

(defn stereo-truth
  "The true disparity of each pixel of the left image of `stereo-pair`."
  [layers w h]
  (img/image w h (fn [x y]
                   (let [cx (+ x 0.5) cy (+ y 0.5)]
                     (:disparity (last (filter (fn [{:keys [box]}]
                                                 (or (nil? box) (let [[x0 y0 x1 y1] box] (and (<= x0 cx) (< cx x1) (<= y0 cy) (< cy y1)))))
                                               layers)))))))
