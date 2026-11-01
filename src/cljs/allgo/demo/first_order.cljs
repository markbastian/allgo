(ns allgo.demo.first-order
  "Demo for `allgo.logic.first-order`: inference in first-order logic,
  after Russell and Norvig, *Artificial Intelligence: A Modern Approach*,
  4th edition, chapter 9, a step at a time.

  Pick a knowledge base and a query, and watch it answered three ways.
  Forward chaining applies each rule to the facts so far, a round at a
  time, until a round adds nothing; the facts it inferred are marked with
  their round. Backward chaining works from the query down: the stack of
  goal lists waiting is shown top first, and each step replaces the top
  one with what each rule matching its first goal would need instead --
  depth first, so Append(x, y, [1, 2]) finds its three answers one after
  another. Resolution negates the query, turns everything into clauses,
  and resolves until the empty clause appears; then the proof is shown,
  each line with the two it came from.

  The crime is the chapter's, as aima-python (MIT) states it; the rabbits
  are aima-python's own."
  (:require [allgo.logic.first-order :as fo]
            ["lil-gui" :default GUI]
            [clojure.string :as str]))

(def ^:private kbs
  {"crime" {:rules ["American(x) & Weapon(y) & Sells(x, y, z) & Hostile(z) => Criminal(x)"
                    "Owns(Nono, M1)" "Missile(M1)"
                    "Missile(x) & Owns(Nono, x) => Sells(West, x, Nono)"
                    "Missile(x) => Weapon(x)" "Enemy(x, America) => Hostile(x)"
                    "American(West)" "Enemy(Nono, America)"]
            :query "Criminal(x)"}
   "rabbits" {:rules ["Farmer(Mac)" "Rabbit(Pete)" "Mother(MrsMac, Mac)" "Mother(MrsRabbit, Pete)"
                      "Rabbit(r) & Farmer(f) => Hates(f, r)" "Mother(m, c) => Loves(m, c)"
                      "Mother(m, r) & Rabbit(r) => Rabbit(m)" "Farmer(f) => Human(f)"
                      "Mother(m, h) & Human(h) => Human(m)"]
              :query "Human(x)"}
   "append" {:rules ["Append(Nil, y, y)" "Append(x, y, z) => Append(Cons(a, x), y, Cons(a, z))"]
             :query "Append(x, y, Cons(1, Cons(2, Nil)))"}})

(def ^:private ^js controls
  #js {:kb     "crime"
       :query  "Criminal(x)"
       :method "backward chaining"
       :speed  2})

(defn- sentence-of
  "A definite clause as a sentence, its variables universal, for
  resolution."
  [{:keys [premises conclusion]}]
  (let [body (if (seq premises)
               [:implies (if (= 1 (count premises)) (first premises) (into [:and] premises)) conclusion]
               conclusion)
        vs (fo/variables body)]
    (if (seq vs) [:forall (vec vs) body] body)))

(defn- esc [s] (-> (str s) (str/replace "&" "&amp;") (str/replace "<" "&lt;") (str/replace ">" "&gt;")))

(defn init! [^js container]
  (let [panel (js/document.createElement "div")
        running? (atom false)
        model (atom nil)]
    (set! (.-style panel)
          (str "position:absolute;inset:44px 250px 12px 12px;overflow:auto;display:grid;"
               "grid-template-columns:minmax(0,1fr) minmax(0,1.4fr);gap:16px;"
               "font:13px/1.5 ui-monospace,Menlo,monospace;color:#c9d2e0"))
    (.appendChild container panel)
    (letfn [(restart! []
              (let [{:keys [rules]} (kbs (.-kb controls))
                    kb (map fo/definite rules)
                    q (try (fo/parse (.-query controls)) (catch :default _ nil))
                    states (when q
                             (case (.-method controls)
                               "forward chaining" (fo/fol-fc-ask-seq kb q)
                               "backward chaining" (fo/fol-bc-ask-seq kb q {:max-depth 30})
                               "resolution" (take 2000 (fo/resolution-seq (map sentence-of kb) q))))]
                (reset! model {:kb kb :rules rules :query q :rest states :st (first states)
                               :steps 0 :answers [] :done? (nil? q) :clock 0 :last-time nil})))
            (advance! []
              (let [now (js/performance.now)
                    {:keys [clock last-time]} @model
                    clock (+ clock (* (.-speed controls) 0.001 (- now (or last-time now))))]
                (swap! model assoc :last-time now)
                (loop [clock clock]
                  (let [{:keys [rest done? query]} @model]
                    (if (and (not done?) (>= clock 1.0))
                      (let [nxt (next rest)
                            st (first (or nxt rest))
                            vs (fo/variables query)]
                        (swap! model (fn [m]
                                       (cond-> (assoc m :rest (or nxt rest) :st st :steps (inc (:steps m))
                                                      :done? (nil? nxt))
                                         (:answer st) (update :answers conj (fo/show-theta (:answer st) vs)))))
                        (recur (- clock 1.0)))
                      (swap! model assoc :clock (min clock 1.0)))))))
            (render! []
              (let [{:keys [rules query st steps answers done?]} @model
                    method (.-method controls)
                    kb-html (str "<h3 style='margin:0 0 6px;color:#7dd3ff'>Knowledge base</h3>"
                                 (str/join "<br>" (map esc rules))
                                 "<h3 style='margin:14px 0 6px;color:#7dd3ff'>Query</h3>"
                                 (if query (esc (fo/show query)) "<span style='color:#ff8080'>could not read it</span>"))
                    right
                    (case method
                      "forward chaining"
                      (str "<h3 style='margin:0 0 6px;color:#7dd3ff'>Facts, round " (:round st) "</h3>"
                           (str/join "<br>" (map #(esc (fo/show %)) (:facts st)))
                           (when (seq (:new st))
                             (str "<br><span style='color:#ffd479'>"
                                  (str/join "<br>" (map #(str "+ " (esc (fo/show %))) (:new st))) "</span>"))
                           (when-let [r (:rule st)]
                             (str "<p style='color:#8b95a8'>this step: "
                                  (esc (str/join " & " (map fo/show (:premises r)))) " => " (esc (fo/show (:conclusion r))) "</p>"))
                           (when (seq (:answers st))
                             (str "<p style='color:#7dffa0'>answers: "
                                  (esc (str/join "  " (distinct (map #(fo/show-theta % (fo/variables query)) (:answers st))))) "</p>")))
                      "backward chaining"
                      (str "<h3 style='margin:0 0 6px;color:#7dd3ff'>Goals waiting, top first</h3>"
                           (str/join "<br>"
                                     (for [[goals theta depth] (take 14 (rseq (:stack st)))]
                                       (str (apply str (repeat depth "&nbsp;"))
                                            (if (seq goals)
                                              (esc (str/join ", " (map #(fo/show (fo/subst theta %)) goals)))
                                              "<span style='color:#7dffa0'>(nothing left: an answer)</span>"))))
                           (when (> (count (:stack st)) 14) (str "<br>... " (- (count (:stack st)) 14) " more"))
                           (when-let [g (:goal st)] (str "<p style='color:#8b95a8'>this step: " (esc (fo/show g)) "</p>"))
                           (when (seq answers)
                             (str "<p style='color:#7dffa0'>answers: " (esc (str/join "  " (distinct answers))) "</p>")))
                      "resolution"
                      (if (:proof st)
                        (str "<h3 style='margin:0 0 6px;color:#7dffa0'>Refuted: the query follows</h3>"
                             (str/join "<br>"
                                       (map-indexed (fn [i {:keys [clause from goal?]}]
                                                      (str i ". " (if (empty? clause) "<b>empty clause</b>"
                                                                      (esc (str/join " | " (map fo/show clause))))
                                                           (cond (seq from) (str " <span style='color:#8b95a8'>from " (str/join ", " from) "</span>")
                                                                 goal? " <span style='color:#ffd479'>negated query</span>"
                                                                 :else " <span style='color:#8b95a8'>given</span>")))
                                                    (fo/proof-lines (:proof st)))))
                        (str "<h3 style='margin:0 0 6px;color:#7dd3ff'>Resolving</h3>"
                             (count (:usable st)) " clauses used, " (count (:queue st)) " waiting"
                             (when-let [g (:given st)] (str "<p>given: " (esc (str/join " | " (map fo/show (:clause g)))) "</p>"))
                             (when (seq (:made st))
                               (str "<p style='color:#ffd479'>"
                                    (str/join "<br>" (map #(esc (str/join " | " (map fo/show (:clause %)))) (take 8 (:made st)))) "</p>")))))]
                (set! (.-innerHTML panel)
                      (str "<div>" kb-html
                           "<p style='margin-top:14px;color:#8b95a8'>" method ", step " steps
                           (when done?
                             (case method
                               "resolution" (cond (:proved? st) " -- proved" (false? (:proved? st)) " -- no proof" :else " -- gave up")
                               " -- finished"))
                           "</p></div><div>" right "</div>"))))
            (frame []
              (when @running?
                (js/requestAnimationFrame frame)
                (advance!)
                (render!)))]
      (let [gui (GUI. #js {:container container})]
        (.onChange (.add gui controls "kb" (clj->js (keys kbs)))
                   (fn [v] (set! (.-query controls) (:query (kbs v))) (restart!)))
        (.onFinishChange (.listen (.add gui controls "query")) restart!)
        (.onChange (.add gui controls "method" #js ["forward chaining" "backward chaining" "resolution"]) restart!)
        (.add gui controls "speed" 0.5 30 0.5)
        (.add gui #js {:restart restart!} "restart"))
      (restart!)
      {:start (fn [] (when-not @running? (reset! running? true) (frame)))
       :stop  (fn [] (reset! running? false))})))

(defonce ^:private controller (atom nil))

(defn start! []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "first-order")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop! []
  (when-let [c @controller] ((:stop c))))
