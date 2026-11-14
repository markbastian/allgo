(ns build
  "Builds the library jar for io.github.markbastian/allgo.

      clojure -T:build jar      ; target/allgo-<version>.jar plus its pom
      clojure -T:build clean

  The jar is src/cljc alone, whose only dependency is Clojure. The demos,
  desktop viewers and demo site sit under the :demo alias and stay out."
  (:require [clojure.tools.build.api :as b]))

(def lib 'io.github.markbastian/allgo)
(def version (format "0.1.%s" (b/git-count-revs nil)))
(def class-dir "target/jar")
(def jar-file (format "target/%s-%s.jar" (name lib) version))
(def src-dirs ["src/cljc"])

;; Top-level deps only, so the pom lists Clojure and nothing from :demo.
(def basis (delay (b/create-basis {:project "deps.edn"})))

(defn clean [_]
  (b/delete {:path class-dir})
  (b/delete {:path jar-file}))

(defn jar [_]
  (clean nil)
  (b/write-pom {:class-dir class-dir
                :lib       lib
                :version   version
                :basis     @basis
                :src-dirs  src-dirs
                :scm       {:url                 "https://github.com/markbastian/allgo"
                            :connection          "scm:git:https://github.com/markbastian/allgo.git"
                            :developerConnection "scm:git:ssh://git@github.com/markbastian/allgo.git"
                            :tag                 (b/git-process {:git-args "rev-parse HEAD"})}
                :pom-data  [[:description "Algorithms in Clojure and ClojureScript."]
                            [:url "https://github.com/markbastian/allgo"]
                            [:licenses
                             [:license
                              [:name "Eclipse Public License 1.0"]
                              [:url "https://www.eclipse.org/legal/epl-v10.html"]]]]})
  (b/copy-dir {:src-dirs   src-dirs
               :target-dir class-dir})
  (b/jar {:class-dir class-dir
          :jar-file  jar-file})
  (println "Built" jar-file))
