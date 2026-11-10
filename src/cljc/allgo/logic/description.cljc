(ns allgo.logic.description
  "Description logic: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 10 (chapter 12 of the 3rd
  edition, whose text this follows) -- the subset of the CLASSIC language
  its figure gives, less SameAs.

    concepts        `:Thing`, a concept name (a keyword), `[:and c ...]`,
                    `[:all role c]`, `[:at-least n role]`, `[:at-most n
                    role]`, `[:fills role individual ...]`, `[:one-of
                    individual ...]`; a TBox maps defined names to their
                    definitions -- `{:Bachelor [:and :Unmarried :Adult
                    :Male]}` -- and names it does not define are primitive
    normal form     a concept expanded and its constraints gathered: the
                    primitives it needs, for each role the least and most
                    fillers allowed, the fillers it must have and the
                    concept every filler must be, and the individuals it
                    may be one of. A concept whose constraints conflict --
                    at least three sons and at most two -- is incoherent
    subsumption     C subsumes D when every D is a C: structurally, when
                    each of C's constraints follows from D's, which takes
                    time polynomial in the two descriptions -- the point of
                    the language. (An incoherent D is subsumed by all.)
    classification  a TBox's names placed in a taxonomy, one a step: each
                    under the most specific names that subsume it; and an
                    individual's membership of a concept, from what is
                    told of it, its role fillers taken as all it has

  The chapter notes that CLASSIC states no subset links: one category
  falls under another only because their descriptions say so."
  (:require [allgo.search :as search]
            [clojure.set :as set]))

(def ^:private top {:prims #{} :roles {} :one-of nil})

(defn- role-constraint [] {:at-least 0 :at-most nil :fills #{} :all nil})

(defn- conjoin
  "The normal form of two normal forms together: their primitives, their
  one-of sets intersected, and their role constraints merged -- the
  tighter of each bound, the fillers of both, and the two `:all`s
  conjoined in turn, so the two definitions recurse into each other."
  [a b]
  (letfn [(merge-role [ra rb]
            {:at-least (max (long (:at-least ra)) (long (:at-least rb)))
             :at-most (cond (nil? (:at-most ra)) (:at-most rb) (nil? (:at-most rb)) (:at-most ra)
                            :else (min (long (:at-most ra)) (long (:at-most rb))))
             :fills (into (:fills ra) (:fills rb))
             :all (cond (nil? (:all ra)) (:all rb) (nil? (:all rb)) (:all ra)
                        :else (both (:all ra) (:all rb)))})
          (both [a b]
            {:prims (into (:prims a) (:prims b))
             :roles (merge-with merge-role (:roles a) (:roles b))
             :one-of (cond (nil? (:one-of a)) (:one-of b) (nil? (:one-of b)) (:one-of a)
                           :else (set/intersection (:one-of a) (:one-of b)))})]
    (both a b)))

(defn- with-role [role f] {:prims #{} :roles {role (f (role-constraint))} :one-of nil})

(defn normalize
  "The normal form of concept `c` under `tbox`: `{:prims #{names} :roles
  {role {:at-least :at-most :fills :all}} :one-of #{individuals}|nil}`,
  each role's `:all` itself a normal form, nil for no constraint."
  [tbox c]
  (cond
    (= c :Thing) top
    (map? c) c
    (keyword? c) (if-let [d (tbox c)] (normalize tbox d) (assoc top :prims #{c}))
    :else
    (let [[op & args] c]
      (case op
        :and (reduce conjoin top (map #(normalize tbox %) args))
        :all (let [[role d] args] (with-role role #(assoc % :all (normalize tbox d))))
        :at-least (let [[n role] args] (with-role role #(assoc % :at-least n)))
        :at-most (let [[n role] args] (with-role role #(assoc % :at-most n)))
        :fills (let [[role & inds] args] (with-role role #(assoc % :fills (set inds) :at-least (count (set inds)))))
        :one-of (assoc top :one-of (set args))))))

(defn incoherent?
  "Whether nothing can be a normal-form concept `n`: a role allowed fewer
  fillers than it must have, no individual left to be, or a filler that
  must be of an incoherent concept."
  [n]
  (or (and (some? (:one-of n)) (empty? (:one-of n)))
      (some (fn [[_ {:keys [at-least at-most fills all]}]]
              (or (and at-most (> (max (long at-least) (count fills)) (long at-most)))
                  (and all (pos? (max (long at-least) (count fills))) (incoherent? all))
                  (and all (seq fills) (some? (:one-of all)) (not (set/subset? fills (:one-of all))))))
            (:roles n))))

(defn- subsumes-normal? [c d]
  (or (incoherent? d)
      (and (set/subset? (:prims c) (:prims d))
           (or (nil? (:one-of c)) (and (some? (:one-of d)) (set/subset? (:one-of d) (:one-of c))))
           (every? (fn [[role rc]]
                     (let [rd (merge (role-constraint) (get-in d [:roles role]))
                           least-d (max (long (:at-least rd)) (count (:fills rd)))]
                       (and (<= (long (:at-least rc)) least-d)
                            (or (nil? (:at-most rc)) (and (:at-most rd) (<= (long (:at-most rd)) (long (:at-most rc)))))
                            (set/subset? (:fills rc) (:fills rd))
                            (or (nil? (:all rc))
                                (= 0 (:at-most rd))
                                (and (:all rd) (subsumes-normal? (:all rc) (:all rd)))
                                (subsumes-normal? (:all rc) top)))))
                   (:roles c)))))

(defn subsumes?
  "Whether concept `c` subsumes concept `d` under `tbox`: every D a C."
  [tbox c d]
  (subsumes-normal? (normalize tbox c) (normalize tbox d)))

(defn equivalent? [tbox c d] (and (subsumes? tbox c d) (subsumes? tbox d c)))

;; ---------------------------------------------------------------------------
;; Classification

(defn- classify-step
  [{:keys [tbox todo parents] :as st}]
  (when-let [concept (first todo)]
    (let [placed (keys parents)
          above (filter #(subsumes? tbox % concept) placed)
          below (filter #(subsumes? tbox concept %) placed)
          ;; Equivalent names sit together; the most specific of the rest
          ;; above are its parents.
          strict-above (remove #(subsumes? tbox concept %) above)
          most-specific (remove (fn [a] (some #(and (not= % a) (subsumes? tbox a %)) strict-above)) strict-above)
          parents (reduce (fn [ps b]
                            ;; A concept now below this one loses the parents
                            ;; this one sits between.
                            (let [bp (ps b)]
                              (if (or (subsumes? tbox b concept)
                                      (some #(and (subsumes? tbox concept %) (not (subsumes? tbox % concept))) bp))
                                ps
                                (assoc ps b (conj (set (remove #(subsumes? tbox % concept) bp)) concept)))))
                          (assoc parents concept (set most-specific))
                          below)]
      (assoc st :todo (rest todo) :name concept :parents parents
             :equivalent (set (filter #(subsumes? tbox concept %) above))))))

(defn classify-seq
  "The names of `tbox` -- and any `extra` primitive names -- placed in a
  taxonomy, a name a step: each state's `:parents` maps every name placed
  so far to the most specific names placed above it, `:name` the one just
  placed and `:equivalent` the names it turned out to mean the same as."
  [tbox & extra]
  (search/step-seq classify-step (classify-step {:tbox tbox :todo (distinct (concat extra (keys tbox))) :parents {}})))

(defn taxonomy
  "Each name of `tbox` and `extra` with its most specific subsumers."
  [tbox & extra]
  (:parents (last (apply classify-seq tbox extra))))

(defn- instance-normal? [abox n ind]
  (let [{:keys [prims fills]} (get abox ind)]
    (and (set/subset? (:prims n) (set prims))
         (or (nil? (:one-of n)) (contains? (:one-of n) ind))
         (every? (fn [[role {:keys [at-least at-most all] rf :fills}]]
                   (let [have (set (get fills role))]
                     (and (<= (long at-least) (count have))
                          (or (nil? at-most) (<= (count have) (long at-most)))
                          (set/subset? rf have)
                          (or (nil? all) (every? #(instance-normal? abox all %) have)))))
                 (:roles n)))))

(defn instance-of?
  "Whether individual `ind` is a `c`, from `abox`, `{individual {:prims
  [names] :fills {role [individuals]}}}`, its fillers taken as all it
  has."
  [tbox abox c ind]
  (instance-normal? abox (normalize tbox c) ind))
