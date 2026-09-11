(ns allgo.array
  "Flat numeric arrays, and the two lines of reader conditional that make
  one on either platform.

  Nine namespaces here had grown a private copy of this apiece. They are
  small enough that the duplication was easy to live with and easy to get
  subtly wrong: the copies had already diverged over whether a constructor
  took a length or a collection, and over whether there was a fill at all.

  The arrays themselves are the point. A simulation that tracks thousands
  of things -- particles, boids, grid cells -- keeps them in one flat array
  per quantity and indexes it, because a vector of maps allocates an object
  per thing per frame and a flat array allocates nothing. A simulation
  that tracks tens of things uses ordinary Clojure values instead, where
  clarity is worth more than the allocation; `allgo.physics.rigid` and
  `allgo.physics.height-field` both make that trade the other way, on
  purpose.")

(defn f32
  "A float array: `n` zeros, or the contents of a collection.

  Both readings, as `clojure.core/float-array` gives its argument."
  ([n-or-coll]
   #?(:clj (float-array n-or-coll)
      :cljs (if (number? n-or-coll)
              (js/Float32Array. n-or-coll)
              (js/Float32Array. (into-array n-or-coll)))))
  ([n fill]
   #?(:clj (float-array n fill)
      :cljs (doto (js/Float32Array. n) (.fill fill)))))

(defn i32
  "An int array: `n` zeros, or the contents of a collection."
  ([n-or-coll]
   #?(:clj (int-array n-or-coll)
      :cljs (if (number? n-or-coll)
              (js/Int32Array. n-or-coll)
              (js/Int32Array. (into-array n-or-coll)))))
  ([n fill]
   #?(:clj (int-array n fill)
      :cljs (doto (js/Int32Array. n) (.fill fill)))))

(defn f64
  "A double array. In ClojureScript this is an ordinary JavaScript array
  rather than a `Float64Array`, which is what `cljs.core/double-array`
  gives and what `aget` and `aset` expect there."
  ([n-or-coll] (double-array n-or-coll))
  ([n fill] (double-array n fill)))

(defn fill!
  "Sets every element of a float array, and returns it."
  [^floats a v]
  (let [v (float v)]
    (dotimes [i (alength a)] (aset a i v))
    a))

(defn ifill!
  "Sets every element of an int array, and returns it."
  [^ints a v]
  (let [v (int v)]
    (dotimes [i (alength a)] (aset a i v))
    a))

(defn copy!
  "Copies `src` into `dst`, element for element, and returns `dst`.

  Not `System/arraycopy`, because this has to hold on both platforms --
  and a typed array's own `set` is not reachable from a Clojure array."
  [^floats dst ^floats src]
  (dotimes [i (alength dst)] (aset dst i (aget src i)))
  dst)
