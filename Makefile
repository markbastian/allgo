# allgo -- common tasks.
#
# The default JDK on this machine is 19, and the Closure compiler that
# shadow-cljs shells out to needs 21 or newer: on 19 a release build dies
# with a NoClassDefFoundError deep in the compiler. So every target that
# touches the JVM runs under JAVA_HOME below rather than whatever `java`
# happens to be first on PATH. Override it if your default is already 21+:
#
#     make release JDK=/path/to/jdk
#
# Deliberately not named JAVA_HOME: make imports the environment, and a
# JAVA_HOME already exported by your shell would win over `?=` here -- so
# this would silently pin whichever JDK it names, which on this machine is
# the 19 that cannot build.
JDK ?= $(shell /usr/libexec/java_home -v 21 2>/dev/null)
ifeq ($(strip $(JDK)),)
$(error No JDK 21+ found. Install one, or pass JDK=/path/to/jdk)
endif

JAVA    := PATH="$(JDK)/bin:$$PATH" JAVA_HOME="$(JDK)"
SHADOW  := $(JAVA) npx shadow-cljs
BUNDLE  := resources/public/js/compiled/allgo.js

.DEFAULT_GOAL := help
.PHONY: help test lint check dev serve compile release repl clean

help: ## List the targets
	@grep -hE '^[a-z-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk -F':.*?## ' '{printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'

test: ## Run the Clojure test suite
	@$(JAVA) clojure -M:test

lint: ## Lint every source and test namespace
	@clj-kondo --lint src test

check: lint test ## Lint and test

dev: ## Build the demo page once, unminified, with source maps
	@$(SHADOW) compile app

serve: ## Watch and serve the demo page on http://localhost:3000
	@$(SHADOW) watch app

release: ## Build the minified demo bundle (advanced optimisations)
	@$(SHADOW) release app
	@printf '\nMinified bundle: %s (%s)\n' "$(BUNDLE)" "$$(ls -lh $(BUNDLE) | awk '{print $$5}')"
	@echo 'Serve resources/public as-is. Run "make dev" to get the'
	@echo 'readable build with source maps back -- release overwrites it.'

repl: ## Start an nREPL with CIDER middleware
	@$(JAVA) clojure -M:nrepl

clean: ## Remove build output and caches
	@rm -rf target .cpcache .shadow-cljs resources/public/js/compiled
