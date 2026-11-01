(ns allgo.logic.categories
  "Categories, semantic networks and inheritance: Russell and Norvig,
  *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 10
  (chapter 12 of the 3rd edition, whose text this follows).

    taxonomy      objects are members of categories (MemberOf links),
                  categories subsets of others (SubsetOf)
    properties    a category's boxed links say something of each of its
                  members -- persons have two legs -- and its double-boxed
                  links the category every value of a relation of a member
                  falls in -- a person's mother is a female person; an
                  object's own links say something of it alone
    inheritance   a property of an object is found by following MemberOf
                  and SubsetOf links up from it. A value on the object
                  itself, or on a category below another that has one,
                  overrides the more general one: persons have two legs,
                  but John has one. When values found on categories
                  neither of which is below the other disagree -- the
                  Nixon diamond, a Quaker and a Republican -- the network
                  cannot say, and `inherit` returns them all
    reification   an n-ary assertion as an event object with a binary
                  link for each argument, Fly(Shankar, NewYork, NewDelhi,
                  Yesterday) as Fly17 with Agent, Origin, Destination and
                  During
    partitions    Disjoint, ExhaustiveDecomposition and Partition on
                  categories given by their members

  A network is `{:member {object #{categories}} :subset {category
  #{categories}} :links {node {relation value}} :boxed {category
  {relation value}} :double-boxed {category {relation category}}}`.")

(defn- up-links [{:keys [member subset]} node]
  (concat (get member node) (get subset node)))

(defn ancestors-of
  "Every category `node` is in or below, nearest first."
  [net node]
  (loop [frontier [node] seen #{node} out []]
    (if-let [n (first frontier)]
      (let [ps (remove seen (up-links net n))]
        (recur (into (subvec frontier 1) ps) (into seen ps) (into out ps)))
      out)))

(defn- inherit-step
  [{:keys [net relation frontier seen found] :as st}]
  (when-let [n (first frontier)]
    (let [v (if (= n (:object st))
              (get-in net [:links n relation])
              (get-in net [:boxed n relation]))
          ps (remove seen (up-links net n))]
      (assoc st :node n :frontier (into (subvec frontier 1) ps) :seen (into seen ps)
             :found (if (some? v) (assoc found n v) found)))))

(defn- decide
  "The values found that nothing more specific overrides: a value on a
  node is overridden by one on any node below it."
  [net found]
  (let [specific (remove (fn [[n _]] (some (fn [[m _]] (and (not= m n) (some #{n} (ancestors-of net m)))) found))
                         found)]
    (into {} specific)))

(defn inherit-seq
  "Inheritance of `relation` for `object`, a node a step: up the MemberOf
  and SubsetOf links breadth first from the object, each state's `:node`
  the one just looked at and `:found` the values seen so far by where
  they were found. The last state's `:values` are those found that
  nothing below them overrides, by node: one value, or several
  disagreeing ones when the network is ambiguous."
  [net object relation]
  (let [states (take-while some? (iterate #(some-> % inherit-step)
                                          (inherit-step {:net net :object object :relation relation
                                                         :frontier [object] :seen #{object} :found {}})))
        end (last states)]
    (concat (butlast states) [(assoc end :values (decide net (:found end)))])))

(defn inherit
  "The value of `relation` for `object`: the value, nil if none is found,
  or `{:ambiguous {node value}}` when the most specific values disagree."
  [net object relation]
  (let [vs (:values (last (inherit-seq net object relation)))
        distinct-vs (distinct (vals vs))]
    (cond (empty? vs) nil
          (= 1 (count distinct-vs)) (first distinct-vs)
          :else {:ambiguous vs})))

(defn value-category
  "The category every value of `relation` for `object` falls in, from a
  double-boxed link on a category it is in, nearest first; nil if none."
  [net object relation]
  (some #(get-in net [:double-boxed % relation]) (ancestors-of net object)))

(defn member?
  "Whether `object` is in `category`, directly or through subsets."
  [net object category]
  (boolean (some #{category} (ancestors-of net object))))

(defn reify-event
  "The network links of an n-ary assertion made an object: `event` in
  `category`, linked to each argument by its role."
  [event category roles args]
  {:member {event #{category}} :links {event (zipmap roles args)}})

(defn merge-networks
  "Networks joined."
  [& nets]
  (apply merge-with (fn [a b] (merge-with (fn [x y] (if (set? x) (into x y) (merge x y))) a b)) nets))

;; ---------------------------------------------------------------------------
;; Disjointness and decomposition, of categories given by their members

(defn disjoint?
  "Disjoint(s): no two of the categories `s`, each a set of members,
  share a member."
  [s]
  (let [s (vec s)]
    (every? empty? (for [i (range (count s)) j (range (inc i) (count s))]
                     (filter (s i) (s j))))))

(defn exhaustive-decomposition?
  "ExhaustiveDecomposition(s, c): every member of `c` is in some category
  of `s`, and every member of those in `c`."
  [s c]
  (= (set c) (reduce into #{} s)))

(defn partition?
  "Partition(s, c): disjoint and exhaustive."
  [s c]
  (and (disjoint? s) (exhaustive-decomposition? s c)))
