.PHONY: all clean debug release test test_sql setup configure lint format

PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
EXTENSION_NAME := raquet
TARGET_DUCKDB_VERSION := v1.2.0
EXT_VERSION := $(shell grep '^version' $(PROJ_DIR)Cargo.toml | head -1 | sed 's/.*"\(.*\)"/\1/')

# Path to the DuckDB CLI used for SQL tests.
# The ./duckdb/ directory is gitignored so you can place any build
# (including pre-release) there.  Override on the command line if needed:
#   make test_sql DUCKDB_BINARY=/usr/local/bin/duckdb
DUCKDB_BINARY ?= $(PROJ_DIR)duckdb/duckdb

all: debug

configure:
	@which cargo > /dev/null 2>&1 || \
	  (echo "ERROR: Rust toolchain not found. Install from https://rustup.rs/" && exit 1)
	@echo "Rust toolchain: $$(cargo --version)"

# Download the release DuckDB CLI into ./duckdb/ (skipped when already present).
# Override TARGET_DUCKDB_VERSION to fetch a different release, e.g.:
#   make setup TARGET_DUCKDB_VERSION=v1.3.0-dev
setup:
	@$(PROJ_DIR)scripts/setup_duckdb.sh "$(TARGET_DUCKDB_VERSION)" "$(PROJ_DIR)"

debug: configure
	cargo build
	@python3 $(PROJ_DIR)scripts/append_metadata.py target/debug/lib$(EXTENSION_NAME).so \
		--duckdb-version $(TARGET_DUCKDB_VERSION) --ext-version $(EXT_VERSION)

release: configure
	cargo build --release
	@python3 $(PROJ_DIR)scripts/append_metadata.py target/release/lib$(EXTENSION_NAME).so \
		--duckdb-version $(TARGET_DUCKDB_VERSION) --ext-version $(EXT_VERSION)

test: configure
	cargo test

test_sql: debug
	@$(PROJ_DIR)scripts/run_sql_tests.sh "$(DUCKDB_BINARY)" "$(PROJ_DIR)"

clean:
	cargo clean

format:
	cargo fmt

lint:
	cargo clippy
