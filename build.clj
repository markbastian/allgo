(ns build
  (:require [clojure.tools.build.api :as b]))

(def class-dir "target/classes")
(def uber-file "target/procedurals.jar")
(def main-ns 'procedurals.swingui)

(defn- basis [] (b/create-basis {:project "deps.edn"}))

(defn clean [_]
  (b/delete {:path "target"}))

(defn uber [_]
  (clean nil)
  (let [basis (basis)]
    (b/copy-dir {:src-dirs   ["src/clj" "src/cljc" "resources"]
                 :target-dir class-dir})
    (b/compile-clj {:basis      basis
                    :src-dirs   ["src/clj" "src/cljc"]
                    :class-dir  class-dir})
    (b/uber {:class-dir class-dir
             :uber-file uber-file
             :basis     basis
             :main      main-ns})))
