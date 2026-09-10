(ns procedurals.perlin
  (:require #?(:clj [clojure.math :as math] :cljs [cljs.math :as math])))

;I can't remember where I originally coded this from but here are some resources:
;  * https://mrl.nyu.edu/~perlin/noise/
;  * https://flafla2.github.io/2014/08/09/perlinnoise.html

(def permutation
  [151 160 137 91 90 15 131 13 201 95 96 53 194 233 7 225 140 36 103 30 69 142 8 99 37 240 21 10 23 190 6 148 247 120
   234 75 0 26 197 62 94 252 219 203 117 35 11 32 57 177 33 88 237 149 56 87 174 20 125 136 171 168 68 175 74 165 71 134
   139 48 27 166 77 146 158 231 83 111 229 122 60 211 133 230 220 105 92 41 55 46 245 40 244 102 143 54 65 25 63 161 1
   216 80 73 209 76 132 187 208 89 18 169 200 196 135 130 116 188 159 86 164 100 109 198 173 186 3 64 52 217 226 250 124
   123 5 202 38 147 118 126 255 82 85 212 207 206 59 227 47 16 58 17 182 189 28 42 223 183 170 213 119 248 152 2 44 154
   163 70 221 153 101 155 167 43 172 9 129 22 39 253 19 98 108 110 79 113 224 232 178 185 112 104 218 246 97 228 251 34
   242 193 238 210 144 12 191 179 162 241 81 51 145 235 249 14 239 107 49 192 214 31 181 199 106 157 184 84 204 176 115
   121 50 45 127 4 150 254 138 236 205 93 222 114 67 29 24 72 243 141 128 195 78 66 215 61 156 180])

(defonce p (into permutation permutation))

(defn fade [t] (* t t t (+ (* t (- (* t 6.0) 15.0)) 10.0)))

(defn lerp [a b t] (+ a (* t (- b a))))

(defn grad [hash x y z]
  (case (bit-and hash 0xF)
    0x0 (+ x y)
    0x1 (- y x)
    0x2 (- x y)
    0x3 (- x (- y))
    0x4 (+ x z)
    0x5 (- z x)
    0x6 (- x z)
    0x7 (- x (- z))
    0x8 (+ y z)
    0x9 (- z y)
    0xA (- y z)
    0xB (- y (- z))
    0xC (+ y x)
    0xD (- x y)
    0xE (- y x)
    0xF (- y (- x))))

(defn- floor [x] (math/floor x))

(defn- mask-idx [x0] #?(:clj (bit-and (long x0) 0xFF) :cljs (bit-and x0 0xFF)))

(defn perlin [x y z]
  (letfn [(f [a] (reduce #(+ (p %1) %2) (conj a 0)))]
    (let [x0 (floor x) y0 (floor y) z0 (floor z)
          xi (mask-idx x0) yi (mask-idx y0) zi (mask-idx z0)
          aaa (f [xi yi zi])
          aba (f [xi (inc yi) zi])
          aab (f [xi yi (inc zi)])
          abb (f [xi (inc yi) (inc zi)])
          baa (f [(inc xi) yi zi])
          bba (f [(inc xi) (inc yi) zi])
          bab (f [(inc xi) yi (inc zi)])
          bbb (f [(inc xi) (inc yi) (inc zi)])
          xf (- x x0) yf (- y y0) zf (- z z0)
          u (fade xf) v (fade yf) w (fade zf)]
      (* 0.5 (inc (lerp (lerp (lerp (grad aaa xf yf zf) (grad baa (dec xf) yf zf) u)
                              (lerp (grad aba xf (dec yf) zf) (grad bba (dec xf) (dec yf) zf) u) v)
                        (lerp (lerp (grad aab xf yf (dec zf)) (grad bab (dec xf) yf (dec zf)) u)
                              (lerp (grad abb xf (dec yf) (dec zf)) (grad bbb (dec xf) (dec yf) (dec zf)) u) v) w))))))

(defn operlin [x y z persistence octaves]
  (loop [total 0.0 freq 1.0 amp 1.0 max-val 0.0 octave 0]
    (if (< octave octaves)
      (recur (+ total (* amp (perlin (* x freq) (* y freq) (* z freq))))
             (* freq 2.0)
             (* amp persistence)
             (+ max-val amp)
             (inc octave))
      (/ total max-val))))
