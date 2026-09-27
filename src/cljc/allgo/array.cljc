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
  purpose."
  #?(:clj (:import (java.util Arrays))))

#?(:cljs
   (defn- typed
     "A typed array of `ctor`'s kind: `n` zeros, or the contents of a
     collection. The one place the three constructors differ is which
     constructor, so in ClojureScript they share this. The JVM's cannot:
     there `float-array`, `int-array` and `double-array` build different
     primitive arrays, and it is their type hints that keep `aget` and
     `aset` unboxed."
     ([ctor n-or-coll]
      (if (number? n-or-coll)
        (new ctor n-or-coll)
        (new ctor (into-array n-or-coll))))
     ([ctor n fill]
      (.fill (new ctor n) fill))))

(defn f32
  "A float array: `n` zeros, or the contents of a collection.

  Both readings, as `clojure.core/float-array` gives its argument."
  ([n-or-coll]
   #?(:clj (float-array n-or-coll) :cljs (typed js/Float32Array n-or-coll)))
  ([n fill]
   #?(:clj (float-array n fill) :cljs (typed js/Float32Array n fill))))

(defn i32
  "An int array: `n` zeros, or the contents of a collection."
  ([n-or-coll]
   #?(:clj (int-array n-or-coll) :cljs (typed js/Int32Array n-or-coll)))
  ([n fill]
   #?(:clj (int-array n fill) :cljs (typed js/Int32Array n fill))))

(defn f64
  "A double array: `n` zeros, or the contents of a collection.

  A real `Float64Array` in ClojureScript rather than the ordinary array
  `cljs.core/double-array` returns. `aget`, `aset` and `alength` work on
  either, and the typed one stores eight bytes a number where the plain
  one stores a pointer to a boxed one."
  (^doubles [n-or-coll]
   #?(:clj (double-array n-or-coll) :cljs (typed js/Float64Array n-or-coll)))
  (^doubles [n fill]
   #?(:clj (double-array n fill) :cljs (typed js/Float64Array n fill))))

;; The rest are each platform's own bulk operations -- `java.util.Arrays`
;; and `System/arraycopy`, a typed array's `fill` and `set` -- which run as
;; one native loop rather than one `aset` at a time. They take the typed
;; arrays the constructors above make; a plain JavaScript array has no
;; `set`.

(defn fill!
  "Sets every element of a float array, and returns it."
  [^floats a v]
  #?(:clj (Arrays/fill a (float v))
     :cljs (.fill a v))
  a)

(defn ifill!
  "Sets every element of an int array, and returns it."
  [^ints a v]
  #?(:clj (Arrays/fill a (int v))
     :cljs (.fill a v))
  a)

(defn copy!
  "Copies float array `src` into `dst`, which is the same length, and
  returns `dst`."
  [^floats dst ^floats src]
  #?(:clj (System/arraycopy src 0 dst 0 (alength src))
     :cljs (.set dst src))
  dst)
