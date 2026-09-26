(ns allgo.app
  "The pages that are one thing rather than a gallery of demos.

  Each has its own HTML page and names itself on the body --
  `<body data-app=\"motorcycle\">` -- and loads the same bundle as the
  demo page; this starts whichever it names, full window, and does
  nothing on a page that names none."
  (:require [allgo.demo.motorcycle :as motorcycle]
            [allgo.game.crossbows :as crossbows]))

(def ^:private apps
  {"crossbows" crossbows/start!
   "motorcycle" motorcycle/start!})

(defonce started
  (when-let [start! (some-> js/document.body .-dataset .-app apps)]
    (start!)))
