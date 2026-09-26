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
.PHONY: help test lint reflect check dev serve bundle compress release repl thumbs clean

help: ## List the targets
	@grep -hE '^[a-z-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk -F':.*?## ' '{printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'

test: ## Run the Clojure test suite
	@$(JAVA) clojure -M:test

lint: ## Lint every source and test namespace
	@clj-kondo --lint src test

reflect: ## Fail if any namespace uses reflection
	@$(JAVA) clojure -M script_reflect.clj

check: lint reflect test ## Lint, reflection check, and tests

dev: ## Build the demo page once, unminified, with source maps
	@$(SHADOW) compile app

serve: ## Watch and serve the demo page on http://localhost:3000
	@$(SHADOW) watch app

release: compress ## Build the minified bundle, pre-compressed
	@echo 'Serve resources/public as-is. Run "make dev" to get the'
	@echo 'readable build with source maps back -- release overwrites it.'

compress: bundle ## Pre-compress the bundle as .gz and .br
	@gzip -9 -c $(BUNDLE) > $(BUNDLE).gz
	@if command -v brotli >/dev/null 2>&1; then brotli -q 11 -c $(BUNDLE) > $(BUNDLE).br; \
	 else echo 'note: brotli not installed, skipping .br (brew install brotli)'; fi
	@printf '\n%-12s %10s\n' ENCODING SIZE
	@printf '%-12s %10s\n' identity "$$(ls -lh $(BUNDLE) | awk '{print $$5}')"
	@printf '%-12s %10s\n' gzip "$$(ls -lh $(BUNDLE).gz | awk '{print $$5}')"
	@[ -f $(BUNDLE).br ] && printf '%-12s %10s\n' br "$$(ls -lh $(BUNDLE).br | awk '{print $$5}')" || true
	@echo
	@echo 'These are only used if the server is told to. See the README:'
	@echo 'nginx gzip_static/brotli_static, Caddy precompressed. Most'
	@echo 'managed hosts compress on the fly and ignore them.'

bundle: ## Build the minified bundle without compressing it
	@# The release build emits one self-contained file, but it writes into
	@# the same directory the development build fills with ~15MB of
	@# per-namespace cljs-runtime files. Those are stale the moment a
	@# release is built and get deployed with it unless cleared first.
	@rm -rf $(dir $(BUNDLE))
	@$(SHADOW) release app

repl: ## Start an nREPL with CIDER middleware
	@$(JAVA) clojure -M:nrepl

thumbs: ## Screenshot every demo for the gallery (needs `make serve` running)
	@node script/thumbnails.mjs

clean: ## Remove build output and caches
	@rm -rf target .cpcache .shadow-cljs resources/public/js/compiled
