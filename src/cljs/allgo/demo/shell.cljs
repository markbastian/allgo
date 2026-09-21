(ns allgo.demo.shell
  "Chrome for the demo page: a picker that shows one demo at a time.

  Demos that own a render loop are started on first selection and stopped when
  you navigate away, so exactly one loop is ever running. That laziness is also
  load-bearing for correctness -- a card hidden with `display:none` measures
  0x0, so a renderer built while hidden would get a degenerate viewport."
  (:require [allgo.demo.arm :as arm]
            [allgo.demo.boids-3d :as boids-3d]
            [allgo.demo.boids-viewer :as boids-2d]
            [allgo.demo.boids-voronoi :as boids-voronoi]
            [allgo.demo.boids-voronoi-3d :as boids-voronoi-3d]
            [allgo.demo.bricks :as bricks]
            [allgo.demo.broad-phase :as broad-phase]
            [allgo.demo.cave :as cave]
            [allgo.demo.cloth :as cloth]
            [allgo.demo.delaunay-viewer :as delaunay]
            [allgo.demo.dla :as dla]
            [allgo.demo.dungeon :as dungeon]
            [allgo.demo.dungeon-boids :as dungeon-boids]
            [allgo.demo.erosion :as erosion]
            [allgo.demo.fire :as fire]
            [allgo.demo.flip :as flip]
            [allgo.demo.fluid :as fluid]
            [allgo.demo.hex :as hex]
            [allgo.demo.human-arm :as human-arm]
            [allgo.demo.island :as island]
            [allgo.demo.joints :as joints]
            [allgo.demo.kepler :as kepler]
            [allgo.demo.lorenz :as lorenz]
            [allgo.demo.orbit-determination :as od]
            [allgo.demo.planet :as planet]
            [allgo.demo.ragdoll :as ragdoll]
            [allgo.demo.rigid :as rigid]
            [allgo.demo.satellite :as satellite]
            [allgo.demo.skinning :as skinning]
            [allgo.demo.soft-body :as soft-body]
            [allgo.demo.solar-system :as solar]
            [allgo.demo.spatial-hash :as spatial-hash]
            [allgo.demo.sphere-fluid :as sphere-fluid]
            [allgo.demo.terrain-webgl :as terrain]
            [allgo.demo.water :as water]
            [reagent.core :as r]
            [reagent.dom.client :as rdom]))

(def ^:private storage-key "allgo.demo")

(def ^:private lifecycles
  {"caves"    {:start cave/start!    :stop cave/stop!}
   "delaunay" {:start delaunay/start! :stop delaunay/stop!}
   "hex"      {:start hex/start!      :stop hex/stop!}
   "soft-body" {:start soft-body/start! :stop soft-body/stop!}
   "skinning" {:start skinning/start! :stop skinning/stop!}
   "cloth"    {:start cloth/start!    :stop cloth/stop!}
   "broad-phase" {:start broad-phase/start! :stop broad-phase/stop!}
   "fluid"    {:start fluid/start!    :stop fluid/stop!}
   "flip"     {:start flip/start!     :stop flip/stop!}
   "fire"     {:start fire/start!     :stop fire/stop!}
   "water"    {:start water/start!    :stop water/stop!}
   "rigid"    {:start rigid/start!    :stop rigid/stop!}
   "bricks"   {:start bricks/start!   :stop bricks/stop!}
   "ragdoll"  {:start ragdoll/start!  :stop ragdoll/stop!}
   "arm"      {:start arm/start!      :stop arm/stop!}
   "human-arm" {:start human-arm/start! :stop human-arm/stop!}
   "joints"   {:start joints/start!   :stop joints/stop!}
   "spatial-hash" {:start spatial-hash/start! :stop spatial-hash/stop!}
   "terrain"  {:start terrain/start!  :stop terrain/stop!}
   "erosion"  {:start erosion/start!  :stop erosion/stop!}
   "island"   {:start island/start!   :stop island/stop!}
   "dla"      {:start dla/start!      :stop dla/stop!}
   "planet"   {:start planet/start!   :stop planet/stop!}
   "sphere-fluid" {:start sphere-fluid/start! :stop sphere-fluid/stop!}
   "dungeon"  {:start dungeon/start! :stop dungeon/stop!}
   "boids"    {:start boids-2d/start! :stop boids-2d/stop!}
   "dungeon-boids" {:start dungeon-boids/start! :stop dungeon-boids/stop!}
   "boids-voronoi" {:start boids-voronoi/start! :stop boids-voronoi/stop!}
   "boids-voronoi-3d" {:start boids-voronoi-3d/start! :stop boids-voronoi-3d/stop!}
   "boids-3d" {:start boids-3d/start! :stop boids-3d/stop!}
   "lorenz"   {:start lorenz/start! :stop lorenz/stop!}
   "kepler"   {:start kepler/start! :stop kepler/stop!}
   "satellite" {:start satellite/start! :stop satellite/stop!}
   "orbit-determination" {:start od/start! :stop od/stop!}
   "solar-system" {:start solar/start! :stop solar/stop!}})

(defn- demos
  "The page itself is the source of truth for which demos exist.

  Order and grouping come from the document too, so adding a demo is
  adding a card -- there is no list here to keep in step with it."
  []
  (vec (for [^js el (array-seq (js/document.querySelectorAll "section.demo-card[data-demo]"))]
         {:id    (.. el -dataset -demo)
          :group (or (.. el -dataset -group) "Other")
          :label (.. el -dataset -label)
          :blurb (.. el -dataset -blurb)
          :el    el})))

(defn- grouped
  "`[[group items] ...]`, each group once, in the order it first appears.

  Not `partition-by`, which only gathers cards that are already next to
  each other: two \"Spatial Queries\" cards with a couple of robotics ones
  between them came out as two runs under one name, so the picker drew
  the heading twice and handed React the same key twice with it. Where a
  card sits in the document is not something the group list should
  depend on."
  [items]
  (let [order (distinct (map :group items))
        by-group (group-by :group items)]
    (mapv (fn [g] [g (vec (by-group g))]) order)))

(defn- matches?
  "Whether a demo answers to what has been typed.

  Matched against the group as well as the name, so that typing \"fluid\"
  finds the fluids and typing \"flock\" finds all five flocking demos
  whatever they are called."
  [q {:keys [label blurb group]}]
  (let [q (.toLowerCase (.trim (or q "")))]
    (or (empty? q)
        (some #(.includes (.toLowerCase (str %)) q) [label blurb group]))))

(defn- remember! [id]
  (try (.setItem js/localStorage storage-key id) (catch :default _ nil)))

(defn- recall []
  (try (.getItem js/localStorage storage-key) (catch :default _ nil)))

(defn- lifecycle! [stage id]
  (when-let [f (get-in lifecycles [id stage])] (f)))

(defn- activate! [demos id]
  (doseq [{other :id ^js el :el} demos]
    (let [active? (= other id)]
      (set! (.-hidden el) (not active?))
      (when active?
        ;; Re-trigger the entrance animation. The re-add has to land on a later
        ;; frame or the browser coalesces it into a no-op.
        (.remove (.-classList el) "is-entering")
        (js/requestAnimationFrame #(.add (.-classList el) "is-entering")))
      (when-not active? (lifecycle! :stop other))))
  (lifecycle! :start id)
  (remember! id))

(defn- icon [cls d]
  [:svg {:class cls :width 15 :height 15 :viewBox "0 0 24 24" :fill "none"
         :stroke "currentColor" :stroke-width 2.5
         :stroke-linecap "round" :stroke-linejoin "round" :aria-hidden true}
   [:polyline {:points d}]])

(defn picker [items]
  (let [state    (r/atom {:open? false
                          :query ""
                          :selected (or (some #{(recall)} (map :id items))
                                        (:id (first items)))})
        nodes    (atom {})
        root     (atom nil)
        close!   #(swap! state assoc :open? false :query "")
        outside! (fn [e] (when (and (:open? @state) @root
                                    (not (.contains @root (.-target e))))
                           (close!)))]
    (r/create-class
     {:component-did-mount
      (fn [_]
        (js/document.addEventListener "mousedown" outside!)
        (activate! items (:selected @state)))

      :component-will-unmount
      (fn [_] (js/document.removeEventListener "mousedown" outside!))

      :reagent-render
      (fn [items]
        (let [{:keys [open? selected query]} @state
              by-id   (into {} (map (juxt :id identity)) items)
              current (get by-id selected (first items))
              ;; Arrow keys walk what is on screen, not what exists, or
              ;; they step into rows the filter has taken away -- and in
              ;; the order it is on screen, which is the grouped order
              ;; rather than the order the cards were declared in.
              groups  (grouped (filterv #(matches? query %) items))
              visible (into [] (mapcat second) groups)
              idx     (fn [id] (or (first (keep-indexed #(when (= (:id %2) id) %1) visible)) 0))
              focus!  (fn [i] (when (seq visible)
                                (some-> (@nodes (mod i (count visible))) .focus)))
              trigger! #(some-> @root (.querySelector ".demo-picker__trigger") .focus)
              choose! (fn [id]
                        (swap! state assoc :selected id :open? false :query "")
                        (activate! items id)
                        (trigger!))
              open!   (fn []
                        (swap! state assoc :open? true)
                        ;; The panel does not exist until this render lands.
                        (js/setTimeout
                         #(some-> @root (.querySelector ".demo-picker__search") .focus)
                         0))]
          [:div.demo-picker {:data-open (str open?) :ref #(reset! root %)}
           [:button.demo-picker__trigger
            {:type          :button
             :aria-haspopup "listbox"
             :aria-expanded (str open?)
             :on-click      #(if open? (close!) (open!))
             :on-key-down   (fn [e]
                              (when (contains? #{"ArrowDown" "ArrowUp"} (.-key e))
                                (.preventDefault e)
                                (open!)))}
            [:span.demo-picker__text
             [:span.demo-picker__group (:group current)]
             [:span.demo-picker__name (:label current)]
             [:span.demo-picker__blurb (:blurb current)]]
            [icon "demo-picker__chevron" "6 9 12 15 18 9"]]
           (when open?
             [:div.demo-picker__panel
              [:div.demo-picker__searchbar
               [:input.demo-picker__search
                {:type        "text"
                 :value       query
                 :placeholder (str "Search " (count items) " demos")
                 :aria-label  "Search demos"
                 :on-change   #(swap! state assoc :query (.. % -target -value))
                 :on-key-down (fn [e]
                                (case (.-key e)
                                  "ArrowDown" (do (.preventDefault e) (focus! 0))
                                  "ArrowUp"   (do (.preventDefault e) (focus! (dec (count visible))))
                                  ;; Enter takes the only thing left, which
                                  ;; is what a search box that narrows to
                                  ;; one result ought to do.
                                  "Enter"     (when (= 1 (count visible))
                                                (.preventDefault e)
                                                (choose! (:id (first visible))))
                                  "Escape"    (do (close!) (trigger!))
                                  nil))}]]
              [:div.demo-picker__list {:role "listbox"}
               (if (empty? visible)
                 [:p.demo-picker__empty "Nothing matches that."]
                 (doall
                  (for [[group run] groups]
                    ^{:key group}
                    ;; A listbox may only contain options and groups, so
                    ;; the block carries the role and the heading is left
                    ;; out of the tree -- otherwise a screen reader either
                    ;; loses the options or reads every heading twice.
                    [:div.demo-picker__group-block {:role "group" :aria-label group}
                     [:div.demo-picker__heading {:aria-hidden true} group]
                     (doall
                      (for [{:keys [id] :as item} run
                            :let [i (idx id)]]
                        ^{:key id}
                        [:button.demo-picker__option
                         {:type          :button
                          :role          "option"
                          :aria-selected (= id selected)
                          :ref           #(swap! nodes assoc i %)
                          :on-click      #(choose! id)
                          :on-key-down   (fn [e]
                                           (case (.-key e)
                                             "ArrowDown" (do (.preventDefault e) (focus! (inc i)))
                                             "ArrowUp"   (do (.preventDefault e) (focus! (dec i)))
                                             "Escape"    (do (close!) (trigger!))
                                             nil))}
                         [icon "demo-picker__check" "20 6 9 17 4 12"]
                         [:span.demo-picker__text
                          [:span.demo-picker__name (:label item)]
                          [:span.demo-picker__blurb (:blurb item)]]]))])))]])]))})))

(defonce ^:private picker-root
  (some-> (js/document.getElementById "demo-picker") rdom/create-root))

(defn main []
  (when picker-root
    (rdom/render picker-root [picker (demos)])))

(main)

(defn on-js-reload []
  (main))
