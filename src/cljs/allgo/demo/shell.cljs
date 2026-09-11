(ns allgo.demo.shell
  "Chrome for the demo page: a picker that shows one demo at a time.

  Demos that own a render loop are started on first selection and stopped when
  you navigate away, so exactly one loop is ever running. That laziness is also
  load-bearing for correctness -- a card hidden with `display:none` measures
  0x0, so a renderer built while hidden would get a degenerate viewport."
  (:require [allgo.demo.boids-3d :as boids-3d]
            [allgo.demo.boids-viewer :as boids-2d]
            [allgo.demo.boids-voronoi :as boids-voronoi]
            [allgo.demo.boids-voronoi-3d :as boids-voronoi-3d]
            [allgo.demo.broad-phase :as broad-phase]
            [allgo.demo.cave :as cave]
            [allgo.demo.cloth :as cloth]
            [allgo.demo.delaunay-viewer :as delaunay]
            [allgo.demo.dungeon :as dungeon]
            [allgo.demo.dungeon-boids :as dungeon-boids]
            [allgo.demo.fire :as fire]
            [allgo.demo.flip :as flip]
            [allgo.demo.fluid :as fluid]
            [allgo.demo.hex :as hex]
            [allgo.demo.kepler :as kepler]
            [allgo.demo.lorenz :as lorenz]
            [allgo.demo.orbit-determination :as od]
            [allgo.demo.satellite :as satellite]
            [allgo.demo.skinning :as skinning]
            [allgo.demo.soft-body :as soft-body]
            [allgo.demo.solar-system :as solar]
            [allgo.demo.spatial-hash :as spatial-hash]
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
   "spatial-hash" {:start spatial-hash/start! :stop spatial-hash/stop!}
   "terrain"  {:start terrain/start!  :stop terrain/stop!}
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
  "The page itself is the source of truth for which demos exist."
  []
  (vec (for [^js el (array-seq (js/document.querySelectorAll "section.demo-card[data-demo]"))]
         {:id    (.. el -dataset -demo)
          :label (.. el -dataset -label)
          :blurb (.. el -dataset -blurb)
          :el    el})))

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

(defn- entry [{:keys [label blurb]}]
  [:span.demo-picker__text
   [:span.demo-picker__name label]
   [:span.demo-picker__blurb blurb]])

(defn picker [items]
  (let [state    (r/atom {:open? false :selected (or (some #{(recall)} (map :id items))
                                                     (:id (first items)))})
        nodes    (atom {})
        root     (atom nil)
        close!   #(swap! state assoc :open? false)
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
        (let [{:keys [open? selected]} @state
              idx     (fn [id] (first (keep-indexed #(when (= (:id %2) id) %1) items)))
              current (nth items (idx selected))
              choose! (fn [id]
                        (swap! state assoc :selected id :open? false)
                        (activate! items id)
                        (some-> @root (.querySelector ".demo-picker__trigger") .focus))
              focus!  (fn [i] (some-> (@nodes (mod i (count items))) .focus))]
          [:div.demo-picker {:data-open (str open?) :ref #(reset! root %)}
           [:button.demo-picker__trigger
            {:type          :button
             :aria-haspopup "listbox"
             :aria-expanded (str open?)
             :on-click      #(swap! state update :open? not)
             :on-key-down   (fn [e]
                              (when (contains? #{"ArrowDown" "ArrowUp"} (.-key e))
                                (.preventDefault e)
                                (swap! state assoc :open? true)
                                (js/setTimeout #(focus! (idx selected)) 0)))}
            [entry current]
            [icon "demo-picker__chevron" "6 9 12 15 18 9"]]
           (when open?
             [:div.demo-picker__panel {:role "listbox"}
              (doall
               (for [[i {:keys [id] :as item}] (map-indexed vector items)]
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
                                      "Escape"    (do (close!)
                                                      (some-> @root (.querySelector ".demo-picker__trigger") .focus))
                                      nil))}
                  [icon "demo-picker__check" "20 6 9 17 4 12"]
                  [entry item]]))])]))})))

(defonce ^:private picker-root
  (some-> (js/document.getElementById "demo-picker") rdom/create-root))

(defn main []
  (when picker-root
    (rdom/render picker-root [picker (demos)])))

(main)

(defn on-js-reload []
  (main))
