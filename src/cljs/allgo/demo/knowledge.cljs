(ns allgo.demo.knowledge
  "Canvas demo for the knowledge representation of Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 10
  (chapter 12 of the 3rd edition), each algorithm a step at a time.

  Semantic network: inheritance climbing figure 12.5's links -- Mary's
  two legs from Persons, John's one his own -- and the Nixon diamond,
  where Quakers and Republicans disagree. Possible worlds: figure 12.4's
  eight worlds and who can tell them apart, a knowledge sentence checked
  by labeling every world where each part of it holds. Event calculus: a
  wumpus hunt as events, its fluents holding and ceasing as time runs.
  Defaults: the Nixon diamond in default logic, its rules tried in every
  order, and in circumscription, its preferred models. Description
  logic: concepts placed in a taxonomy from their definitions alone.
  Truth maintenance: the 2048 Olympic sites -- the chapter's example, the
  consequences ours -- with Athletics moved from Bucharest to Sibiu, and
  the conclusions that survive kept."
  (:require [allgo.logic.categories :as cat]
            [allgo.logic.defaults :as df]
            [allgo.logic.description :as dl]
            [allgo.logic.events :as ev]
            [allgo.logic.modal :as md]
            [allgo.logic.tms :as tms]
            ["lil-gui" :default GUI]
            [clojure.math :as math]
            [clojure.string :as str]))

(def ^:private ^js controls
  #js {:method "semantic network" :query "Mary: Legs" :sentence "K_L (K_S I | K_S ~I)" :speed 3})

(def ^:private method-names
  #js ["semantic network" "possible worlds" "event calculus" "defaults" "description logic" "truth maintenance"])

;; ---------------------------------------------------------------------------
;; Semantic networks

(def ^:private figure-12-5
  {:member {:Mary #{:FemalePersons} :John #{:MalePersons}}
   :subset {:FemalePersons #{:Persons} :MalePersons #{:Persons} :Persons #{:Mammals}}
   :links {:Mary {:SisterOf :John} :John {:Legs 1}}
   :boxed {:Persons {:Legs 2}}
   :double-boxed {:Persons {:HasMother :FemalePersons}}})

(def ^:private diamond
  {:member {:Nixon #{:Quakers :Republicans}}
   :subset {:Quakers #{:Persons} :Republicans #{:Persons}}
   :boxed {:Quakers {:Pacifist true} :Republicans {:Pacifist false}}})

(def ^:private queries
  {"Mary: Legs" [figure-12-5 :Mary :Legs]
   "John: Legs" [figure-12-5 :John :Legs]
   "Nixon: Pacifist" [diamond :Nixon :Pacifist]})

(def ^:private layouts
  {figure-12-5 {:Mammals [0.5 0.12] :Persons [0.5 0.32] :FemalePersons [0.28 0.54] :MalePersons [0.72 0.54]
                :Mary [0.28 0.8] :John [0.72 0.8]}
   diamond {:Persons [0.5 0.15] :Quakers [0.28 0.45] :Republicans [0.72 0.45] :Nixon [0.5 0.78]}})

;; ---------------------------------------------------------------------------
;; Possible worlds: figure 12.4, bottom

(defn- both-ways [pairs]
  (reduce (fn [m [a b]] (-> m (update a (fnil conj #{}) b) (update b (fnil conj #{}) a))) {} pairs))

(defn- all-of [ws] (for [x ws y ws :when (not= x y)] [x y]))

(def ^:private superman-pairs (concat (all-of [:w0 :w2 :w4 :w6]) (all-of [:w1 :w3 :w5 :w7])))
(def ^:private lois-pairs (concat (all-of [:w0 :w1 :w2 :w3]) [[:w4 :w5] [:w6 :w7]]))

(def ^:private worlds-model
  (md/model {:w0 [:I :R] :w1 [:R] :w2 [:I] :w3 [] :w4 [:I :R] :w5 [:R] :w6 [:I] :w7 []}
            {:S (both-ways superman-pairs) :L (both-ways lois-pairs)}))

(def ^:private sentences
  {"K_S I" [:knows :S :I]
   "K_S R" [:knows :S :R]
   "K_L R" [:knows :L :R]
   "K_L (K_S I | K_S ~I)" [:knows :L [:or [:knows :S :I] [:knows :S [:not :I]]]]
   "K_S (K_L R | K_L ~R)" [:knows :S [:or [:knows :L :R] [:knows :L [:not :R]]]]
   "K_L R | K_L ~R" [:or [:knows :L :R] [:knows :L [:not :R]]]})

(defn- show-modal [s]
  (cond
    (keyword? s) (name s)
    :else (let [[op & args] s
                wrap #(if (#{:and :or :implies} (when (vector? %) (first %))) (str "(" (show-modal %) ")") (show-modal %))]
            (case op
              :knows (str "K_" (name (first args)) " " (wrap (second args)))
              :not (str "~" (wrap (first args)))
              :and (str/join " & " (map wrap args))
              :or (str/join " | " (map wrap args))
              :implies (str (wrap (first args)) " => " (wrap (second args)))))))

(def ^:private world-at
  {:w0 [0.12 0.3] :w1 [0.36 0.3] :w2 [0.12 0.62] :w3 [0.36 0.62]
   :w4 [0.62 0.3] :w5 [0.86 0.3] :w6 [0.62 0.62] :w7 [0.86 0.62]})

;; ---------------------------------------------------------------------------
;; Event calculus: a wumpus hunt, our own

(def ^:private hunt
  {:events [{:id "Start" :interval [0 0]}
            {:id "Shoot" :category :Shootings :interval [5 6]}
            {:id "Scream" :interval [6 6]}
            {:id "Grab" :interval [11 12]}
            {:id "Climb" :interval [16 17]}]
   :initiates (fn [e _] (case (:id e) "Start" [:HaveArrow :WumpusAlive :InCave] "Grab" [:HaveGold] []))
   :terminates (fn [e _] (case (:id e) "Shoot" [:HaveArrow] "Scream" [:WumpusAlive] "Climb" [:InCave] []))})

(def ^:private hunt-fluents [:HaveArrow :WumpusAlive :HaveGold :InCave])
(def ^:private hunt-end 20)

;; ---------------------------------------------------------------------------
;; Defaults

(def ^:private nixon-facts ["Republican(Nixon) & Quaker(Nixon)"])
(def ^:private nixon-defaults ["Republican(x) : ~Pacifist(x) / ~Pacifist(x)" "Quaker(x) : Pacifist(x) / Pacifist(x)"])
(def ^:private nixon-kb ["Republican(Nixon) & Quaker(Nixon)"
                         "Republican(x) & ~Abnormal2(x) => ~Pacifist(x)"
                         "Quaker(x) & ~Abnormal3(x) => Pacifist(x)"])

;; ---------------------------------------------------------------------------
;; Description logic

(def ^:private tbox
  {:Man [:and :Adult :Male]
   :Woman [:and :Adult :Female]
   :Bachelor [:and :Unmarried :Adult :Male]
   :Father [:and :Man [:at-least 1 :Son]]
   :Mother [:and :Woman [:at-least 1 :Son]]
   :Parent [:at-least 1 :Son]
   :Busy [:and :Man [:at-least 3 :Son] [:at-most 2 :Daughter]
          [:all :Son [:and :Unemployed :Married [:all :Spouse :Doctor]]]
          [:all :Daughter [:and :Professor [:fills :Department :Physics :Math]]]]})

(defn- show-concept [c]
  (cond
    (keyword? c) (name c)
    :else (let [[op & args] c]
            (case op
              :and (str "And(" (str/join ", " (map show-concept args)) ")")
              :all (str "All(" (name (first args)) ", " (show-concept (second args)) ")")
              :at-least (str "AtLeast(" (first args) ", " (name (second args)) ")")
              :at-most (str "AtMost(" (first args) ", " (name (second args)) ")")
              :fills (str "Fills(" (str/join ", " (map name args)) ")")
              :one-of (str "OneOf(" (str/join ", " (map name args)) ")")))))

;; ---------------------------------------------------------------------------
;; Truth maintenance: the Olympic sites

(def ^:private swim "Site(Swimming, Pitesti)")
(def ^:private bucharest "Site(Athletics, Bucharest)")
(def ^:private sibiu "Site(Athletics, Sibiu)")
(def ^:private arad "Site(Equestrian, Arad)")

(def ^:private olympic-justifications
  [{:antecedents #{bucharest} :consequent "UsesNationalStadium"}
   {:antecedents #{sibiu} :consequent "BuildStadium(Sibiu)"}
   {:antecedents #{bucharest swim} :consequent "ShortTransfers"}
   {:antecedents #{sibiu arad} :consequent "ShortTransfers"}
   {:antecedents #{swim arad} :consequent "WesternVenues"}
   {:antecedents #{"ShortTransfers"} :consequent "OneAthletesVillage"}
   {:antecedents #{"OneAthletesVillage" "WesternVenues"} :consequent "VillageNearArad"}])

(def ^:private olympic-at
  {swim [0.17 0.25] bucharest [0.17 0.45] sibiu [0.17 0.65] arad [0.17 0.85]
   "UsesNationalStadium" [0.5 0.35] "BuildStadium(Sibiu)" [0.5 0.55] "ShortTransfers" [0.5 0.75]
   "WesternVenues" [0.5 0.92] "OneAthletesVillage" [0.8 0.6] "VillageNearArad" [0.8 0.85]})

(defn- olympic-script []
  (let [t0 (tms/jtms olympic-justifications)
        run (fn [[states t] [label f node]]
              (let [ss (map #(assoc % :label label :target node) (f t node))]
                [(into states ss) (last ss)]))]
    (first (reduce run [[] t0]
                   [["TELL " tms/enable-seq swim] ["TELL " tms/enable-seq arad] ["TELL " tms/enable-seq bucharest]
                    ["RETRACT " tms/retract-seq bucharest] ["TELL " tms/enable-seq sibiu]]))))

(def ^:private car
  (tms/atms [:BatteryDead :StarterBroken :NoFuel] [:LightsOn]
            (concat (tms/horn "BatteryDead => NoCrank" "StarterBroken => NoCrank"
                              "NoCrank => WontStart" "NoFuel => WontStart")
                    [{:antecedents #{:LightsOn :BatteryDead} :consequent :contradiction}])))

;; ---------------------------------------------------------------------------

(defn- init [method]
  (case method
    "semantic network" (let [[net obj rel] (queries (.-query controls))]
                         {:net net :states (vec (cat/inherit-seq net obj rel)) :object obj :relation rel})
    "possible worlds" (let [s (sentences (.-sentence controls))]
                        {:sentence s :states (vec (md/check-seq worlds-model s))})
    "event calculus" {:states (vec (take (inc hunt-end) (ev/timeline-seq hunt 0)))}
    "defaults" {:states (vec (df/extensions-seq nixon-facts nixon-defaults))
                :preferred (df/preferred-models nixon-kb #{:Abnormal2 :Abnormal3})
                :prioritized (df/preferred-models nixon-kb [#{:Abnormal3} #{:Abnormal2}])}
    "description logic" {:states (vec (dl/classify-seq tbox :Adult :Male :Female))}
    "truth maintenance" {:states (olympic-script) :atms (vec (tms/atms-seq car))}))

(defn init! [^js container]
  (let [canvas (js/document.createElement "canvas")
        ctx (.getContext canvas "2d")
        caption (js/document.createElement "div")
        running? (atom false)
        model (atom nil)]
    (set! (.-style canvas) "display:block")
    (set! (.-style caption)
          (str "position:absolute;left:12px;bottom:8px;right:12px;padding:3px 8px;border-radius:4px;"
               "background:rgba(10,14,22,0.75);font:11px ui-monospace,Menlo,monospace;color:#c9d2e0;"
               "pointer-events:none;white-space:pre-wrap"))
    (doto container (.appendChild canvas) (.appendChild caption))
    (letfn [(restart! [] (reset! model (assoc (init (.-method controls)) :i 0 :clock 0 :last nil)))
            (advance! []
              (let [now (js/performance.now)
                    {:keys [clock last states i]} @model
                    clock (+ clock (* (.-speed controls) 0.001 (- now (or last now))))
                    n (count states)]
                ;; A step per tick of the clock; at the end, a pause, then
                ;; from the start again.
                (if (>= clock 1)
                  (swap! model assoc :clock 0 :last now :i (if (>= i (+ n (* 2 (.-speed controls)))) 0 (inc i)))
                  (swap! model assoc :clock clock :last now))))
            (current [] (let [{:keys [states i]} @model] (nth states (min i (dec (count states))))))
            (text! [s x y color & [size align]]
              (set! (.-fillStyle ctx) color) (set! (.-font ctx) (str (or size 12) "px ui-monospace,Menlo,monospace"))
              (set! (.-textAlign ctx) (or align "left"))
              (.fillText ctx s x y)
              (set! (.-textAlign ctx) "left"))
            (box! [label x y fill stroke & [size]]
              (set! (.-font ctx) (str (or size 13) "px ui-monospace,Menlo,monospace"))
              (let [w (+ 16 (.-width (.measureText ctx label))) h 26]
                (set! (.-fillStyle ctx) fill) (set! (.-strokeStyle ctx) stroke) (set! (.-lineWidth ctx) 1.5)
                (.beginPath ctx) (.roundRect ctx (- x (/ w 2)) (- y (/ h 2)) w h 6) (.fill ctx) (.stroke ctx)
                (text! label x (+ y 4) "#e6edf7" (or size 13) "center")))
            (line! [[x1 y1] [x2 y2] color & [dash width]]
              (set! (.-strokeStyle ctx) color) (set! (.-lineWidth ctx) (or width 1.5))
              (.setLineDash ctx (clj->js (or dash [])))
              (.beginPath ctx) (.moveTo ctx x1 y1) (.lineTo ctx x2 y2) (.stroke ctx)
              (.setLineDash ctx #js []))
            (arrow! [[x1 y1] [x2 y2] color label]
              (let [dx (- x2 x1) dy (- y2 y1) len (math/sqrt (+ (* dx dx) (* dy dy)))
                    ux (/ dx len) uy (/ dy len)
                    [sx sy] [(+ x1 (* ux 16)) (+ y1 (* uy 16))] [ex ey] [(- x2 (* ux 16)) (- y2 (* uy 16))]]
                (line! [sx sy] [ex ey] color)
                (set! (.-fillStyle ctx) color)
                (.beginPath ctx) (.moveTo ctx ex ey)
                (.lineTo ctx (- ex (* 8 ux) (* 4 uy)) (+ (- ey (* 8 uy)) (* 4 ux)))
                (.lineTo ctx (+ (- ex (* 8 ux)) (* 4 uy)) (- ey (* 8 uy) (* 4 ux))) (.fill ctx)
                (when label (text! label (+ (/ (+ sx ex) 2) 6) (/ (+ sy ey) 2) "#7f8da3" 11))))
            (draw! []
              (let [cw (.-width canvas) ch (.-height canvas)
                    method (.-method controls)]
                (set! (.-fillStyle ctx) "#0d1118") (.fillRect ctx 0 0 cw ch)
                (set! (.-textContent caption)
                      (str method "  |  "
                           (case method
                             "semantic network" (draw-network! cw ch)
                             "possible worlds" (draw-worlds! cw ch)
                             "event calculus" (draw-events! cw ch)
                             "defaults" (draw-defaults! cw ch)
                             "description logic" (draw-taxonomy! cw ch)
                             "truth maintenance" (draw-tms! cw ch))))))
            (draw-network! [cw ch]
              (let [{:keys [net object relation]} @model
                    st (current)
                    at (into {} (map (fn [[k [fx fy]]] [k [(* fx cw) (+ 40 (* fy (- ch 120)))]])) (layouts net))
                    done? (contains? st :values)]
                (doseq [[from tos] (:member net) to tos] (arrow! (at from) (at to) "#5b6b85" "MemberOf"))
                (doseq [[from tos] (:subset net) to tos] (arrow! (at from) (at to) "#5b6b85" "SubsetOf"))
                (doseq [[from links] (:links net) [rel v] links :when (keyword? v) :when (at v)]
                  (line! (at from) (at v) "#8f6bd8" [4 4]) (text! (name rel) (/ (+ (first (at from)) (first (at v))) 2) (- (second (at from)) 10) "#8f6bd8" 11 "center"))
                (doseq [[node [x y]] at]
                  (let [visited? (contains? (:seen st) node)
                        here? (= node (:node st))
                        found (get (:found st) node)
                        winner? (and done? (contains? (:values st) node))]
                    (box! (name node) x y (cond winner? "#1f5a37" here? "#5a4a1f" visited? "#202a3a" :else "#161c27")
                          (cond winner? "#7dffa0" here? "#ffd479" :else "#3a4760"))
                    (doseq [[i [rel v]] (map-indexed vector (concat (get-in net [:boxed node]) (remove (comp keyword? val) (get-in net [:links node]))))]
                      (text! (str (when (get-in net [:boxed node rel]) "each member: ") (name rel) " " v)
                             x (+ y 30 (* 14 i)) (if (= rel relation) (if (some? found) "#ffd479" "#c9d2e0") "#7f8da3") 11 "center"))
                    (doseq [[rel c] (get-in net [:double-boxed node])]
                      (text! (str "each member's " (name rel) " in " (name c)) x (- y 22) "#7f8da3" 11 "center"))))
                (str (name object) "'s " (name relation) ": "
                     (if done?
                       (let [vs (distinct (vals (:values st)))]
                         (cond (empty? vs) "none found" (= 1 (count vs)) (str (first vs) ", from " (str/join ", " (map name (keys (:values st)))))
                               :else (str "ambiguous -- " (str/join ", " (map (fn [[k v]] (str (name k) " says " v)) (:values st))))))
                       (str "looking at " (name (:node st))))
                     "\nup the MemberOf and SubsetOf links from the object; a value below another overrides it")))
            (draw-worlds! [cw ch]
              (let [{:keys [sentence states i]} @model
                    st (current)
                    at (into {} (map (fn [[k [fx fy]]] [k [(* fx cw) (+ 30 (* fy (- ch 100)))]])) world-at)
                    final? (>= i (dec (count states)))]
                (doseq [[a b] superman-pairs :when (neg? (compare a b))] (line! (at a) (at b) "rgba(125,196,255,0.55)" nil 2))
                (doseq [[a b] lois-pairs :when (neg? (compare a b))]
                  ;; Nudged aside, so that a pair both can't tell apart shows both lines.
                  (let [[x1 y1] (at a) [x2 y2] (at b) l (max 1 (math/sqrt (+ (* (- x2 x1) (- x2 x1)) (* (- y2 y1) (- y2 y1)))))
                        nx (* 5 (/ (- y1 y2) l)) ny (* 5 (/ (- x2 x1) l))]
                    (line! [(+ x1 nx) (+ y1 ny)] [(+ x2 nx) (+ y2 ny)] "rgba(255,143,200,0.85)" [3 5] 2)))
                (doseq [[w [x y]] at]
                  (let [on? (contains? (:worlds st) w)
                        atoms ((:worlds worlds-model) w)]
                    (set! (.-fillStyle ctx) (if on? "#1f5a37" "#161c27"))
                    (set! (.-strokeStyle ctx) (if on? "#7dffa0" "#3a4760")) (set! (.-lineWidth ctx) 2)
                    (.beginPath ctx) (.arc ctx x y 34 0 (* 2 math/PI)) (.fill ctx) (.stroke ctx)
                    (text! (name w) x (- y 6) "#e6edf7" 13 "center")
                    (text! (str (if (atoms :I) "I" "~I") "," (if (atoms :R) "R" "~R")) x (+ y 12) "#9fb4d8" 12 "center")))
                (text! "solid blue: worlds Superman can't tell apart     dotted pink: Lois" 40 (- ch 70) "#9fb4d8" 11)
                (text! (str "step " (min (inc i) (count states)) " of " (count states) ":  " (show-modal (:sentence st))) 40 24 "#ffd479" 14)
                (str (show-modal sentence) " holds in " (if final? (let [ws (:worlds (last states))] (if (seq ws) (str/join " " (map name (sort ws))) "no world")) "...")
                     "\nI: Superman is Clark Kent, R: rain is forecast; each step labels the worlds where one part of the sentence holds, innermost first")))
            (draw-events! [cw _ch]
              (let [st (current)
                    t (:t st)
                    x0 160 w (- cw 200) dx (/ w (inc hunt-end))
                    tx (fn [t] (+ x0 (* t dx)))]
                (doseq [{:keys [id interval]} (:events hunt)]
                  (let [[a b] interval]
                    (set! (.-fillStyle ctx) (if (<= a t) "#5a4a1f" "#202a3a"))
                    (.fillRect ctx (tx a) 170 (max 6 (* dx (- b a))) 22)
                    (text! id (tx a) 160 "#ffd479" 12)))
                (doseq [[row f] (map-indexed vector hunt-fluents)]
                  (let [y (+ 240 (* 60 row))]
                    (text! (name f) 20 (+ y 25) "#c9d2e0" 13)
                    (doseq [s (range (inc hunt-end))]
                      (let [v (ev/holds hunt f s)]
                        (set! (.-fillStyle ctx) (cond (> s t) "#161c27" (true? v) "#2f8f57" (false? v) "#8f3a3a" :else "#2a3140"))
                        (.fillRect ctx (+ (tx s) 1) y (- dx 2) 40)))))
                (line! [(+ (tx t) (/ dx 2)) 150] [(+ (tx t) (/ dx 2)) (+ 250 (* 60 (count hunt-fluents)))] "#ffd479")
                (doseq [s (range 0 (inc hunt-end) 5)] (text! (str s) (tx s) (+ 265 (* 60 (count hunt-fluents))) "#7f8da3" 11))
                (text! (str "Before(Shoot, Grab): " (ev/before? [5 6] [11 12]) "    Meet(Shoot, Scream): " (ev/meet? [5 6] [6 6]))
                       20 (+ 300 (* 60 (count hunt-fluents))) "#9fb4d8" 12)
                (str "t = " t "  |  holds: " (str/join ", " (map name (sort (:holds st))))
                     "\ngreen T(f, t), red ~T(f, t), grey neither -- nothing has initiated or terminated it yet; a fluent holds after an event initiates it until one terminates it")))
            (draw-defaults! [_cw _ch]
              (let [{:keys [states i preferred prioritized]} @model
                    shown (take (inc i) states)]
                (text! "Default logic: Republican(Nixon) & Quaker(Nixon)" 40 40 "#ffd479" 17)
                (doseq [[k d] (map-indexed vector nixon-defaults)] (text! d 60 (+ 70 (* 24 k)) "#c9d2e0" 15))
                (text! "processes tried, depth first:" 40 140 "#9fb4d8" 15)
                (doseq [[k st] (map-indexed vector shown)]
                  (text! (str (if (empty? (:process st)) "(nothing applied yet)" (str/join "  then  " (map :text (:process st))))
                              (cond (:failed? st) "  -- failed" (:closed? st) "  -- closed: an extension" :else ""))
                         60 (+ 170 (* 26 k)) (cond (:failed? st) "#ff8f70" (:closed? st) "#7dffa0" :else "#c9d2e0") 15))
                (let [exts (:extensions (last shown)) y 300]
                  (text! (str "extensions: " (count exts)) 40 y "#9fb4d8" 15)
                  (doseq [[k e] (map-indexed vector exts)]
                    (text! (str "{" (str/join ", " (map #(if (keyword? %) (name %) (str "~" (name (second %)))) (:conclusions e))) "}")
                           60 (+ y 28 (* 26 k)) "#7dffa0" 16)))
                (let [y 420
                      show-model (fn [m] (str/join ", " (sort (keep (fn [[k v]] (when v (name k))) m))))]
                  (text! "Circumscription of Abnormal2 and Abnormal3: the preferred models" 40 y "#ffd479" 17)
                  (doseq [[k m] (map-indexed vector preferred)] (text! (show-model m) 60 (+ y 30 (* 26 k)) "#c9d2e0" 15))
                  (text! "prioritized, religion first (Abnormal3 minimized before Abnormal2):" 40 (+ y 100) "#9fb4d8" 15)
                  (doseq [[k m] (map-indexed vector prioritized)] (text! (show-model m) 60 (+ y 130 (* 26 k)) "#7dffa0" 15)))
                "is Nixon a pacifist?\ntwo extensions, two preferred models: the reasoner stays agnostic until told which rule to prefer"))
            (draw-taxonomy! [cw ch]
              (let [st (current)
                    parents (:parents st)
                    depth (fn depth [n] (if (empty? (parents n)) 0 (inc (apply max (map depth (parents n))))))
                    levels (group-by depth (sort (keys parents)))
                    at (into {} (for [[d ns] levels [k n] (map-indexed vector ns)]
                                  [n [(* cw (/ (inc k) (inc (count ns)))) (+ 160 (* d (/ (- ch 300) 3)))]]))]
                (doseq [[n ps] parents p ps] (arrow! (at n) (at p) "#5b6b85" nil))
                (doseq [[n [x y]] at]
                  (box! (name n) x y (if (= n (:name st)) "#5a4a1f" "#202a3a") (if (= n (:name st)) "#ffd479" "#3a4760")))
                (when-let [n (:name st)]
                  (text! (str (name n) " = " (if-let [d (tbox n)] (show-concept d) "a primitive"))
                         20 (- ch 60) "#ffd479" 12))
                (str "placed " (count parents) " of " (+ 3 (count tbox))
                     (when-let [n (:name st)] (str "  |  " (name n) " under " (if (seq (parents n)) (str/join ", " (map name (sort (parents n)))) "nothing")))
                     "\neach concept under the most specific ones whose definitions subsume its own -- no subset links are ever stated")))
            (draw-tms! [cw ch]
              (let [{:keys [atms i states]} @model
                    st (current)
                    at (into {} (map (fn [[k [fx fy]]] [k [(* fx (* 0.68 cw)) (+ 20 (* fy (- ch 140)))]])) olympic-at)]
                (doseq [{:keys [antecedents consequent]} olympic-justifications a antecedents]
                  (arrow! (at a) (at consequent) (if (and ((:in st) a) ((:in st) consequent)) "#2f8f57" "#3a4760") nil))
                (doseq [[n [x y]] at]
                  (let [in? ((:in st) n) changed? (some #{n} (:changed st))]
                    (box! n x y (cond changed? "#5a4a1f" in? "#1f5a37" :else "#161c27") (cond changed? "#ffd479" in? "#7dffa0" :else "#3a4760") 11)))
                (text! (str (:label st) (:target st)) 20 24 "#ffd479" 13)
                (let [x (* 0.72 cw) {:keys [labels nogoods]} (nth atms (min (quot i 2) (dec (count atms))))]
                  (text! "ATMS: why won't the car start?" x 170 "#ffd479" 13)
                  (text! "the lights are on" x 192 "#9fb4d8" 11)
                  (doseq [[k n] (map-indexed vector [:WontStart :NoCrank :BatteryDead :StarterBroken :NoFuel])]
                    (text! (str (name n) ": " (str/join " " (map #(str "{" (str/join "," (map name %)) "}") (get labels n))))
                           x (+ 220 (* 22 k)) "#c9d2e0" 12))
                  (text! (str "nogoods: " (str/join " " (map #(str "{" (str/join "," (map name %)) "}") nogoods))) x 340 "#ff8f70" 12))
                (str "step " (min (inc i) (count states)) " of " (count states) "  |  " (:label st) (:target st)
                     (when (:phase st) (str "  |  " (case (:phase st) :out "marking out what depended on " :support "checking for other support: " :in "bringing in what follows from ") (:node st)))
                     "\nJTMS: green in, yellow just changed; retracting Bucharest takes out what rested on it, then Sibiu's justifications bring the rest back without re-deriving")))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (advance!)
                (draw!)))
            (resize! []
              (let [cw (.-clientWidth container) ch (.-clientHeight container)]
                (when (and (pos? cw) (pos? ch))
                  (set! (.-width canvas) cw)
                  (set! (.-height canvas) ch))))]
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "method" method-names) restart!)
        (.onChange (.add gui controls "query" (clj->js (keys queries))) restart!)
        (.onChange (.add gui controls "sentence" (clj->js (keys sentences))) restart!)
        (.add gui controls "speed" 0.5 12 0.5)
        (.add gui #js {:restart restart!} "restart"))
      (.observe (js/ResizeObserver. resize!) container)
      (resize!)
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (resize!) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "knowledge")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
