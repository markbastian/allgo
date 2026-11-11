(ns allgo.demo.registry
  "Every demo, by the id `resources/public/demos.json` knows it by.

  Each demo draws into the element whose id is its own, and is started and
  stopped by these. `allgo.app` uses them to run any one of them full
  window on the player page, `play/?demo=<id>`; the gallery that links
  there is `index.html`, and needs none of this -- it is plain HTML that
  loads without the bundle."
  (:require [allgo.demo.antialias :as antialias]
            [allgo.demo.arm :as arm]
            [allgo.demo.boids-3d :as boids-3d]
            [allgo.demo.boids-viewer :as boids-2d]
            [allgo.demo.boids-voronoi :as boids-voronoi]
            [allgo.demo.boids-voronoi-3d :as boids-voronoi-3d]
            [allgo.demo.bricks :as bricks]
            [allgo.demo.broad-phase :as broad-phase]
            [allgo.demo.car-paths :as car-paths]
            [allgo.demo.cave :as cave]
            [allgo.demo.city :as city]
            [allgo.demo.city-agents :as city-agents]
            [allgo.demo.city-plan :as city-plan]
            [allgo.demo.classifiers :as classifiers]
            [allgo.demo.cloth :as cloth]
            [allgo.demo.constraints :as constraints]
            [allgo.demo.contract-net :as contract-net]
            [allgo.demo.delaunay-viewer :as delaunay]
            [allgo.demo.dla :as dla]
            [allgo.demo.dungeon :as dungeon]
            [allgo.demo.dungeon-boids :as dungeon-boids]
            [allgo.demo.erosion :as erosion]
            [allgo.demo.fire :as fire]
            [allgo.demo.first-order :as first-order]
            [allgo.demo.flip :as flip]
            [allgo.demo.fluid :as fluid]
            [allgo.demo.game-search :as game-search]
            [allgo.demo.generative :as generative]
            [allgo.demo.genetic :as genetic]
            [allgo.demo.geometry :as geometry]
            [allgo.demo.gesture :as gesture]
            [allgo.demo.grid-navigation :as grid-navigation]
            [allgo.demo.grid-world :as grid-world]
            [allgo.demo.hex :as hex]
            [allgo.demo.human-arm :as human-arm]
            [allgo.demo.island :as island]
            [allgo.demo.joints :as joints]
            [allgo.demo.kepler :as kepler]
            [allgo.demo.knowledge :as knowledge]
            [allgo.demo.language :as language]
            [allgo.demo.local-search :as local-search]
            [allgo.demo.localization :as localization]
            [allgo.demo.lorenz :as lorenz]
            [allgo.demo.many-robots :as many-robots]
            [allgo.demo.motion-planning :as motion-planning]
            [allgo.demo.motorcycle :as motorcycle]
            [allgo.demo.orbit-determination :as od]
            [allgo.demo.planet :as planet]
            [allgo.demo.probabilistic :as probabilistic]
            [allgo.demo.ragdoll :as ragdoll]
            [allgo.demo.rigid :as rigid]
            [allgo.demo.robot-control :as robot-control]
            [allgo.demo.satellite :as satellite]
            [allgo.demo.search :as search]
            [allgo.demo.skinning :as skinning]
            [allgo.demo.slam :as slam]
            [allgo.demo.soft-body :as soft-body]
            [allgo.demo.solar-system :as solar]
            [allgo.demo.spatial-hash :as spatial-hash]
            [allgo.demo.sphere-fluid :as sphere-fluid]
            [allgo.demo.surflets :as surflets]
            [allgo.demo.terrain-webgl :as terrain]
            [allgo.demo.textures :as textures]
            [allgo.demo.three-body :as three-body]
            [allgo.demo.uncertain-decisions :as uncertain-decisions]
            [allgo.demo.vision :as vision]
            [allgo.demo.volumes :as volumes]
            [allgo.demo.water :as water]
            [allgo.demo.wumpus :as wumpus]))

(def lifecycles
  "id -> `{:start :stop}`. Adding a demo is a line here and an entry in
  `demos.json`."
  {"caves"    {:start cave/start!    :stop cave/stop!}
   "delaunay" {:start delaunay/start! :stop delaunay/stop!}
   "hex"      {:start hex/start!      :stop hex/stop!}
   "search"   {:start search/start!   :stop search/stop!}
   "local-search" {:start local-search/start! :stop local-search/stop!}
   "game-search" {:start game-search/start! :stop game-search/stop!}
   "constraints" {:start constraints/start! :stop constraints/stop!}
   "wumpus"   {:start wumpus/start!   :stop wumpus/stop!}
   "first-order" {:start first-order/start! :stop first-order/stop!}
   "grid-world" {:start grid-world/start! :stop grid-world/stop!}
   "classifiers" {:start classifiers/start! :stop classifiers/stop!}
   "generative" {:start generative/start! :stop generative/stop!}
   "vision" {:start vision/start! :stop vision/stop!}
   "knowledge" {:start knowledge/start! :stop knowledge/stop!}
   "contract-net" {:start contract-net/start! :stop contract-net/stop!}
   "uncertain-decisions" {:start uncertain-decisions/start! :stop uncertain-decisions/stop!}
   "probabilistic" {:start probabilistic/start! :stop probabilistic/stop!}
   "language" {:start language/start! :stop language/stop!}
   "localization" {:start localization/start! :stop localization/stop!}
   "motion-planning" {:start motion-planning/start! :stop motion-planning/stop!}
   "car-paths" {:start car-paths/start! :stop car-paths/stop!}
   "grid-navigation" {:start grid-navigation/start! :stop grid-navigation/stop!}
   "many-robots" {:start many-robots/start! :stop many-robots/stop!}
   "slam" {:start slam/start! :stop slam/stop!}
   "robot-control" {:start robot-control/start! :stop robot-control/stop!}
   "soft-body" {:start soft-body/start! :stop soft-body/stop!}
   "skinning" {:start skinning/start! :stop skinning/stop!}
   "cloth"    {:start cloth/start!    :stop cloth/stop!}
   "broad-phase" {:start broad-phase/start! :stop broad-phase/stop!}
   "fluid"    {:start fluid/start!    :stop fluid/stop!}
   "flip"     {:start flip/start!     :stop flip/stop!}
   "fire"     {:start fire/start!     :stop fire/stop!}
   "water"    {:start water/start!    :stop water/stop!}
   "rigid"    {:start rigid/start!    :stop rigid/stop!}
   "bricks"   {:start bricks/start!   :stop bricks/stop!}
   "ragdoll"  {:start ragdoll/start!  :stop ragdoll/stop!}
   "motorcycle" {:start motorcycle/start! :stop motorcycle/stop!}
   "arm"      {:start arm/start!      :stop arm/stop!}
   "human-arm" {:start human-arm/start! :stop human-arm/stop!}
   "joints"   {:start joints/start!   :stop joints/stop!}
   "spatial-hash" {:start spatial-hash/start! :stop spatial-hash/stop!}
   "terrain"  {:start terrain/start!  :stop terrain/stop!}
   "erosion"  {:start erosion/start!  :stop erosion/stop!}
   "island"   {:start island/start!   :stop island/stop!}
   "dla"      {:start dla/start!      :stop dla/stop!}
   "planet"   {:start planet/start!   :stop planet/stop!}
   "textures" {:start textures/start! :stop textures/stop!}
   "volumes"  {:start volumes/start!  :stop volumes/stop!}
   "geometry" {:start geometry/start! :stop geometry/stop!}
   "genetic"  {:start genetic/start!  :stop genetic/stop!}
   "antialias" {:start antialias/start! :stop antialias/stop!}
   "surflets" {:start surflets/start! :stop surflets/stop!}
   "gesture"  {:start gesture/start!  :stop gesture/stop!}
   "city-plan" {:start city-plan/start! :stop city-plan/stop!}
   "city"     {:start city/start!     :stop city/stop!}
   "city-agents" {:start city-agents/start! :stop city-agents/stop!}
   "sphere-fluid" {:start sphere-fluid/start! :stop sphere-fluid/stop!}
   "dungeon"  {:start dungeon/start! :stop dungeon/stop!}
   "boids"    {:start boids-2d/start! :stop boids-2d/stop!}
   "dungeon-boids" {:start dungeon-boids/start! :stop dungeon-boids/stop!}
   "boids-voronoi" {:start boids-voronoi/start! :stop boids-voronoi/stop!}
   "boids-voronoi-3d" {:start boids-voronoi-3d/start! :stop boids-voronoi-3d/stop!}
   "boids-3d" {:start boids-3d/start! :stop boids-3d/stop!}
   "lorenz"   {:start lorenz/start! :stop lorenz/stop!}
   "kepler"   {:start kepler/start! :stop kepler/stop!}
   "three-body" {:start three-body/start! :stop three-body/stop!}
   "satellite" {:start satellite/start! :stop satellite/stop!}
   "orbit-determination" {:start od/start! :stop od/stop!}
   "solar-system" {:start solar/start! :stop solar/stop!}})
