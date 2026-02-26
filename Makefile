.PHONY: all clean debug release test test_sql

# Default target
all: release

# Directory structure
PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
BUILD_DIR := build
DUCKDB_DIR := duckdb

# Build configurations
release:
	mkdir -p $(BUILD_DIR)/release && \
	cd $(BUILD_DIR)/release && \
	cmake -DCMAKE_BUILD_TYPE=Release \
		-DDUCKDB_EXTENSION_CONFIGS="$(PROJ_DIR)extension_config.cmake" \
		-DEXTENSION_STATIC_BUILD=1 \
		../../$(DUCKDB_DIR) && \
	cmake --build . --config Release -- -j$(shell nproc 2>/dev/null || sysctl -n hw.ncpu)

debug:
	mkdir -p $(BUILD_DIR)/debug && \
	cd $(BUILD_DIR)/debug && \
	cmake -DCMAKE_BUILD_TYPE=Debug \
		-DDUCKDB_EXTENSION_CONFIGS="$(PROJ_DIR)extension_config.cmake" \
		-DEXTENSION_STATIC_BUILD=1 \
		../../$(DUCKDB_DIR) && \
	cmake --build . --config Debug -- -j$(shell nproc 2>/dev/null || sysctl -n hw.ncpu)

# Testing
test: release
	cd $(BUILD_DIR)/release && \
	./duckdb -c "LOAD raquet; SELECT quadbin_from_tile(0,0,0);"

test_sql: release
	./$(BUILD_DIR)/release/test/unittest "test/sql/**"

# Cleanup
clean:
	rm -rf $(BUILD_DIR)

# Format code
format:
	find src -name "*.cpp" -o -name "*.hpp" | xargs clang-format -i

# Setup: clone DuckDB if not present
setup:
	@if [ ! -d "$(DUCKDB_DIR)" ]; then \
		echo "Cloning DuckDB..."; \
		git clone --depth 1 https://github.com/duckdb/duckdb.git $(DUCKDB_DIR); \
	else \
		echo "DuckDB already exists"; \
	fi
