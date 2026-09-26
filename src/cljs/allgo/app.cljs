(ns allgo.app
  "The pages that are one thing rather than a gallery of demos.

  Each has its own HTML page and names itself on the body --
  `<body data-app=\"motorcycle\">` -- and loads the same bundle as the
  demo page; this starts whichever it names, full window, and does
  nothing on a page that names none."
  (:require [allgo.demo.motorcycle :as motorcycle]
            [allgo.demo.solar-system :as solar-system]
            [allgo.game.crossbows :as crossbows]))

(def ^:private apps
  {"crossbows" crossbows/start!
   "motorcycle" motorcycle/start!
   "solar-system" solar-system/start!})

(defonce started
  ;; `getAttribute`, not `(.. body -dataset -app)`: a release build
  ;; renames `app`, which is this page's name and no browser's, and the
  ;; page then starts nothing. A method the browser defines keeps its name.
  (when-let [start! (some-> js/document.body (.getAttribute "data-app") apps)]
    (start!)))
