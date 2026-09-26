(ns allgo.app
  "The pages that run one thing full window rather than a gallery.

  Each has its own HTML page and names itself on the body --
  `<body data-app=\"motorcycle\">` -- and loads the same bundle; this
  starts whichever it names, and does nothing on a page that names none.
  `play` is the general one: `play/?demo=<id>` runs any demo in
  `allgo.demo.registry`, which is where every gallery tile goes that has
  no page of its own."
  (:require [allgo.demo.motorcycle :as motorcycle]
            [allgo.demo.registry :as registry]
            [allgo.demo.solar-system :as solar-system]
            [allgo.game.crossbows :as crossbows]))

(defn- play!
  "Runs the demo `?demo=` names, full window.

  The demo draws into the element whose id is its own, so that element is
  made here before it is started. Its name comes from `demos.json`, the
  same list the gallery is drawn from, and goes in the page's title bar;
  it arrives after the demo has started, which does not wait for it."
  []
  (let [id (.get (js/URLSearchParams. (.. js/window -location -search)) "demo")
        ^js host (js/document.getElementById "play")
        ^js title (js/document.getElementById "play-title")]
    (if-let [{:keys [start]} (get registry/lifecycles id)]
      (let [^js el (doto (js/document.createElement "div")
                     (-> .-id (set! id))
                     (-> .-className (set! "play-viewport")))]
        (.appendChild host el)
        (start)
        ;; A phone's width is not enough for the scene and the controls
        ;; both; the panel starts closed, a tap on its title away.
        (when (.-matches (js/window.matchMedia "(max-width: 640px)"))
          (some-> (.querySelector el ".lil-gui.lil-root > .lil-title") (.click)))
        (-> (js/fetch "../demos.json")
            (.then (fn [^js r] (.json r)))
            (.then (fn [^js catalog]
                     ;; `aget`, not `.-label`: these are the JSON's own
                     ;; names, which a release build would otherwise rename.
                     (when-let [^js d (some #(when (= id (aget % "id")) %)
                                            (array-seq (aget catalog "demos")))]
                       (set! (.-title js/document) (str (aget d "label") " · allgo"))
                       (set! (.-textContent title) (aget d "label"))
                       (set! (.-title title) (aget d "blurb")))))
            (.catch (fn [_] nil))))
      (set! (.-textContent host)
            (str "There is no demo called \"" id "\".")))))

(def ^:private apps
  {"crossbows" crossbows/start!
   "motorcycle" motorcycle/start!
   "solar-system" solar-system/start!
   "play" play!})

(defonce started
  ;; `?thumb` in the address takes the page's own chrome away -- title
  ;; bar, controls, readouts -- leaving only the scene, which is what the
  ;; gallery's thumbnails are screenshots of. The stylesheets do the
  ;; hiding; this only says so.
  ;;
  ;; `getAttribute`, not `(.. body -dataset -app)`: a release build
  ;; renames `app`, which is this page's name and no browser's, and the
  ;; page then starts nothing. A method the browser defines keeps its name.
  (when-let [start! (some-> js/document.body (.getAttribute "data-app") apps)]
    (when (.has (js/URLSearchParams. (.. js/window -location -search)) "thumb")
      (.add (.-classList js/document.documentElement) "thumb"))
    (start!)))
