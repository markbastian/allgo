(ns allgo.procedural.terrain
  (:require #?(:clj [clojure.java.io :as io]))
  (:import #?(:clj (java.awt Color))
           #?(:clj (java.awt.image BufferedImage))
           #?(:clj (javax.imageio ImageIO))))

;https://danielbeard.wordpress.com/2010/08/07/terrain-generation-and-smoothing/

(defn double-indices-1d [cells]
  (zipmap (map #(* 2 %) (keys cells)) (vals cells)))

(defn double-indices-2d [{:keys [cells] :as m}]
  (-> m
      (assoc :cells (zipmap (map (fn [[x y]] [(* 2 x) (* 2 y)]) (keys cells)) (vals cells)))
      (update :width * 0.5)))

(defn diamond [{:keys [cells width dim] :as m}]
  (let [dim (+ dim (dec dim))]
    (-> m
        (assoc :dim dim)
        (update :cells into (for [i (range 1 dim 2) j (range 1 dim 2)
                                  :let
                                  [is ((juxt inc dec dec inc) i)
                                   js ((juxt inc inc dec dec) j)
                                   c (map vector is js)
                                   vs (map cells c)
                                   x (/ (reduce + vs) (count vs))]]
                              [[i j] (+ x (* width (dec (* 2.0 (rand)))))])))))

(defn square [{:keys [cells width dim] :as m}]
  (update m :cells into
          (for [i (range 0 dim) j (range 0 dim)
                :when (not (cells [i j]))
                :let
                [is ((juxt inc identity dec identity) i)
                 js ((juxt identity inc identity dec) j)
                 c (map vector is js)
                 vs (filter identity (map cells c))
                 x (/ (reduce + vs) (count vs))]]
            [[i j] (+ x (* width (dec (* 2.0 (rand)))))])))

(def step (comp square diamond double-indices-2d))

(defn init-grid
  "Seeds the 2x2 corner grid the diamond-square algorithm grows from.
  corners are given in [top-left top-right bottom-right bottom-left] order."
  [{:keys [width corners] :or {corners [0.0 0.0 0.0 0.0]}}]
  (let [[c00 c01 c11 c10] corners]
    {:cells {[0 0] c00, [0 1] c01, [1 1] c11, [1 0] c10}
     :width width
     :dim   2}))

(defn generate
  "Runs the diamond-square algorithm for `iterations` steps, producing a
  heightmap grid of size (2^iterations + 1) on a side."
  [{:keys [iterations] :as opts}]
  (nth (iterate step (init-grid opts)) iterations))

(defn cells->grid
  "Converts the {:cells {[i j] height ...} :dim n} sparse map produced by
  `generate` into a dim x dim vector-of-vectors, indexed [row][col]."
  [{:keys [cells dim]}]
  (vec (for [i (range dim)]
         (vec (for [j (range dim)]
                (get cells [i j]))))))

(defn heights
  "Flat seq of every height value in a grid produced by `generate`."
  [{:keys [cells]}]
  (vals cells))

#?(:clj (defn create-image-map [{:keys [cells dim]}]
          (let [img (BufferedImage. dim dim BufferedImage/TYPE_INT_RGB)
                lo (apply min (map second cells))
                hi (apply max (map second cells))]
            (doseq [[[i j] v] cells
                    :let [c (float (/ (- v lo) (- hi lo)))]]
              (.setRGB img i j (.getRGB (Color. c c c))))
            (ImageIO/write img "png" (io/file "img.png")))))

(comment
  (-> (generate {:width 1.0 :iterations 8 :corners [0.0 0.5 0.0 0.5]})
      create-image-map))
