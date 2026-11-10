(ns allgo.spatial.pairs
  "What the broad phases share: whether two axis-aligned boxes overlap,
  and a list of overlapping pairs laid out as a compressed adjacency
  list -- the form `allgo.spatial.hash/adjacency` gives too.

  Boxes are flat arrays, three numbers a box: `mins` and `maxs`."
  (:require [allgo.array :as a]))

(defn overlap?
  "Do boxes `i` and `j` overlap on all three axes?"
  [^doubles mins ^doubles maxs ^long i ^long j]
  (let [a (* 3 i) b (* 3 j)]
    (and (<= (aget mins a) (aget maxs b))
         (>= (aget maxs a) (aget mins b))
         (<= (aget mins (+ a 1)) (aget maxs (+ b 1)))
         (>= (aget maxs (+ a 1)) (aget mins (+ b 1)))
         (<= (aget mins (+ a 2)) (aget maxs (+ b 2)))
         (>= (aget maxs (+ a 2)) (aget mins (+ b 2))))))

(defn adjacency
  "The overlapping `pairs` of `n` objects, each `[i j]`, as a compressed
  adjacency list: `{:starts :ids :pairs}`, object `i`'s partners
  `ids[starts[i] .. starts[i+1])`, each pair recorded once under the
  higher-numbered object."
  [pairs ^long n]
  (let [by-hi (reduce (fn [m [i j]] (update m j (fnil conj []) i)) {} pairs)
        ^ints starts (a/i32 (inc n))]
    (loop [i 0 acc (transient []) total 0]
      (if (= i n)
        (do (aset starts n (int total))
            {:starts starts
             :ids (let [v (persistent! acc)
                        ^ints out (a/i32 (count v))]
                    (dotimes [k (count v)] (aset out k (int (v k))))
                    out)
             :pairs total})
        (let [_ (aset starts i (int total))
              partners (sort (get by-hi i))]
          (recur (inc i) (reduce conj! acc partners) (+ total (count partners))))))))
