#!/usr/bin/env bash
# run_sql_tests.sh - Run DuckDB sqllogictest files for the raquet extension.
#
# Usage: ./scripts/run_sql_tests.sh [duckdb_binary] [proj_dir]
#
# The DuckDB binary defaults to ./duckdb/duckdb (populated by 'make setup').
# Run 'make setup' first if the binary is not present, or place a custom
# (pre-release) binary at ./duckdb/duckdb manually.
#
# The script parses each .test file in test/sql/ and executes `query` and
# `statement` blocks, comparing query output against expected values.

set -euo pipefail

DUCKDB_BINARY="${1:-./duckdb/duckdb}"
PROJ_DIR="${2:-$(cd "$(dirname "$0")/.." && pwd)}"

EXT_PATH="${PROJ_DIR}/target/debug/libraquet.so"
TEST_DIR="${PROJ_DIR}/test/sql"

# ---------------------------------------------------------------------------
# 1. Locate the DuckDB CLI
# ---------------------------------------------------------------------------
locate_duckdb() {
    if [[ -x "${DUCKDB_BINARY}" ]]; then
        echo "Using DuckDB: $(${DUCKDB_BINARY} --version 2>&1 | head -1)"
        echo "  Path: ${DUCKDB_BINARY}"
        return
    fi
    # Also search PATH
    if command -v "${DUCKDB_BINARY}" &>/dev/null; then
        echo "Using DuckDB: $(${DUCKDB_BINARY} --version 2>&1 | head -1)"
        return
    fi

    echo "ERROR: DuckDB binary not found at '${DUCKDB_BINARY}'."
    echo ""
    echo "Run 'make setup' to download the default release into ./duckdb/, or"
    echo "place a custom (pre-release) binary there manually, or override:"
    echo "  make test_sql DUCKDB_BINARY=/path/to/duckdb"
    exit 1
}

# ---------------------------------------------------------------------------
# 2. Verify the extension .so exists
# ---------------------------------------------------------------------------
check_extension() {
    if [[ ! -f "${EXT_PATH}" ]]; then
        echo "ERROR: Extension not found at ${EXT_PATH}"
        echo "Run 'make debug' (or 'cargo build') first."
        exit 1
    fi
}

# ---------------------------------------------------------------------------
# 3. Run a single .test file
# ---------------------------------------------------------------------------
PASS=0
FAIL=0

run_test_file() {
    local test_file="$1"
    local test_name
    test_name="$(basename "${test_file}")"

    local -a queries=()
    local -a expecteds=()
    local -a types=()

    local cur_sql=""
    local cur_expected=""
    local cur_type=""
    local reading_expected=0

    while IFS= read -r line || [[ -n "$line" ]]; do
        # Skip comment lines and 'require' directives
        if [[ "$line" =~ ^#.*$ || "$line" =~ ^require\ .*$ ]]; then
            continue
        fi

        # Start of a query block: "query I", "query III", "query R", etc.
        if [[ "$line" =~ ^query\ [A-Z] ]]; then
            cur_type="query"
            cur_sql=""
            cur_expected=""
            reading_expected=0
            continue
        fi

        # Start of a statement block: "statement ok" or "statement error"
        if [[ "$line" =~ ^statement\ (ok|error) ]]; then
            cur_type="statement"
            cur_sql=""
            cur_expected=""
            reading_expected=0
            continue
        fi

        # Separator between SQL and expected output
        if [[ "$line" == "----" ]] && [[ -n "$cur_type" ]]; then
            reading_expected=1
            continue
        fi

        # Empty line ends a test block
        if [[ -z "$line" ]] && [[ -n "$cur_type" ]]; then
            if [[ -n "$cur_sql" ]]; then
                queries+=("$cur_sql")
                expecteds+=("$cur_expected")
                types+=("$cur_type")
            fi
            cur_type=""
            cur_sql=""
            cur_expected=""
            reading_expected=0
            continue
        fi

        # Accumulate SQL or expected-output lines
        if [[ -n "$cur_type" ]]; then
            if [[ $reading_expected -eq 1 ]]; then
                cur_expected="${cur_expected:+${cur_expected}$'\n'}${line}"
            else
                cur_sql="${cur_sql:+${cur_sql}$'\n'}${line}"
            fi
        fi
    done < "${test_file}"

    # Flush the last block (file may not end with a blank line)
    if [[ -n "$cur_type" && -n "$cur_sql" ]]; then
        queries+=("$cur_sql")
        expecteds+=("$cur_expected")
        types+=("$cur_type")
    fi

    local num_tests=${#queries[@]}
    local file_pass=0
    local file_fail=0

    for ((i = 0; i < num_tests; i++)); do
        local sql="${queries[$i]}"
        local expected="${expecteds[$i]}"
        local type="${types[$i]}"

        local init_sql
        init_sql="LOAD '${EXT_PATH}'; ${sql}"

        local actual
        actual="$("${DUCKDB_BINARY}" -unsigned -csv -noheader -separator $'\t' \
            -c "${init_sql}" 2>&1)" || true

        if [[ "$type" == "statement" ]]; then
            if echo "$actual" | grep -qi "error"; then
                echo "  FAIL [${test_name}#$((i+1))]: ${sql%%$'\n'*}"
                echo "       Error: ${actual}"
                file_fail=$((file_fail + 1))
            else
                file_pass=$((file_pass + 1))
            fi
        else
            local actual_trimmed expected_trimmed
            actual_trimmed="$(printf '%s' "$actual" | sed 's/[[:space:]]*$//' | grep -v '^$')"
            expected_trimmed="$(printf '%s' "$expected" | sed 's/[[:space:]]*$//' | grep -v '^$')"

            if [[ "$actual_trimmed" == "$expected_trimmed" ]]; then
                file_pass=$((file_pass + 1))
            else
                echo "  FAIL [${test_name}#$((i+1))]: ${sql%%$'\n'*}"
                echo "       Expected: $(printf '%s' "$expected_trimmed" | head -3)"
                echo "       Got:      $(printf '%s' "$actual_trimmed" | head -3)"
                file_fail=$((file_fail + 1))
            fi
        fi
    done

    PASS=$((PASS + file_pass))
    FAIL=$((FAIL + file_fail))

    if [[ $file_fail -eq 0 ]]; then
        echo "  PASS ${test_name} (${file_pass} checks)"
        return 0
    else
        echo "  FAIL ${test_name} (${file_pass} passed, ${file_fail} failed)"
        return 1
    fi
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
locate_duckdb
check_extension

echo ""
echo "Running SQL tests from ${TEST_DIR}/"
echo "Extension:  ${EXT_PATH}"
echo ""

overall_fail=0
for test_file in "${TEST_DIR}"/*.test; do
    if ! run_test_file "${test_file}"; then
        overall_fail=1
    fi
done

echo ""
echo "Results: ${PASS} passed, ${FAIL} failed"

if [[ $overall_fail -ne 0 ]]; then
    echo "FAILED"
    exit 1
else
    echo "OK"
fi
