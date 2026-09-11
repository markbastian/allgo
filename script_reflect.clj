;; Load every namespace with reflection warnings on and fail if any appear.
;; Reflection is silent, legal, and can cost two orders of magnitude, so it
;; is worth a build step rather than a note in a docstring.
(require '[clojure.string :as str])
(let [nses (->> (file-seq (java.io.File. "src"))
                (map #(.getPath ^java.io.File %))
                (filter #(or (str/ends-with? % ".clj") (str/ends-with? % ".cljc")))
                ;; src/cljs is ClojureScript-only; reflection does not apply.
                (remove #(str/includes? % "src/cljs/"))
                (map #(-> %
                          (str/replace #"^src/(clj|cljc)/" "")
                          (str/replace #"\.cljc?$" "")
                          (str/replace "/" ".")
                          (str/replace "_" "-")
                          symbol))
                sort)
      out (java.io.StringWriter.)]
  (binding [*warn-on-reflection* true *err* out]
    (doseq [n nses]
      (try (require n :reload)
           (catch Exception e (println "FAILED to load" n ":" (.getMessage e))))))
  (let [warnings (->> (str/split-lines (str out))
                      (filter #(str/includes? % "Reflection warning"))
                      ;; Only our own namespaces: third-party code gets
                      ;; loaded too and we cannot fix theirs.
                      (filter #(str/includes? % "allgo/"))
                      distinct)]
    (println (count nses) "namespaces checked")
    (if (seq warnings)
      (do (println (count warnings) "reflection warnings:")
          (doseq [w warnings] (println " " w))
          (System/exit 1))
      (println "no reflection warnings"))))
