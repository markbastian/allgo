(ns allgo.desktop.swingui
  (:require [allgo.procedural.terrain :as proc])
  (:import [java.awt BorderLayout]
           [javax.swing Box BoxLayout JFrame JSlider]
           [javax.swing.event ChangeListener]
           [org.jfree.chart ChartFactory ChartPanel]
           [org.jfree.data.xy XYSeries XYSeriesCollection])
  (:gen-class))

(defn add-data [series data]
  (.clear series)
  (doseq [h (->> data :cells (sort-by first))]
    (.add series (first h) (second h))))

(defn run-app [exit-behavior]
  (let [data-series (XYSeries. "Height")
        iterations-slider (JSlider. 1 10 5)
        left-slider (JSlider. -100 100 0)
        right-slider (JSlider. -100 100 0)
        roughness-slider (JSlider. 0 1000 200)
        f #(doto data-series (add-data (proc/generate
                                        {:width      (-> roughness-slider .getValue (/ 100.0))
                                         :cells      {0 (-> left-slider .getValue (/ 100.0))
                                                      1 (-> right-slider .getValue (/ 100.0))}
                                         :iterations (-> iterations-slider .getValue)})))
        cl (reify ChangeListener (stateChanged [_ _] (f)))]
    (doto (JFrame. "Fractal Line")
      (.setSize 800 600)
      (.setDefaultCloseOperation exit-behavior)
      (.add (ChartPanel.
             (ChartFactory/createXYLineChart
              "Fractal Terrain"
              "Width"
              "Height"
              (XYSeriesCollection. (f)))) BorderLayout/CENTER)
      (.add (doto (Box. BoxLayout/Y_AXIS)
              (.add (doto iterations-slider (.addChangeListener cl)))
              (.add (doto left-slider (.addChangeListener cl)))
              (.add (doto right-slider (.addChangeListener cl)))
              (.add (doto roughness-slider (.addChangeListener cl)))) BorderLayout/SOUTH)
      (.setVisible true))))

(defn -main [] (run-app JFrame/EXIT_ON_CLOSE))
