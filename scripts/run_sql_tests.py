#!/usr/bin/env python3
"""
Runner for DuckDB sqllogictest-format .test files.

Handles the subset of sqllogictest directives used in this project:
  require <ext>   - load an extension
  statement ok    - run SQL, expect success
  query [types]   - run SQL, compare tabular results
  # comment lines - ignored
"""

import argparse
import glob
import os
import sys
import tempfile

import duckdb


def load_extension(con, ext_name, ext_dir):
    """Load an extension by name.

    For 'raquet' (the primary extension under test), looks for a locally-built
    .duckdb_extension file under ext_dir.  For all other names, attempts DuckDB
    auto-load and logs a warning if that fails.
    """
    if ext_name == "raquet":
        pattern = os.path.join(ext_dir, "extension", "raquet", "raquet.duckdb_extension")
        matches = glob.glob(pattern)
        if matches:
            con.execute(f"LOAD '{matches[0]}'")
            return
        raise RuntimeError(f"Could not find raquet extension under {ext_dir}")
    try:
        con.execute(f"LOAD {ext_name}")
    except Exception as e:
        # parquet is bundled; other extensions may simply not be installed
        print(f"  Warning: could not load extension '{ext_name}': {e}", file=sys.stderr)


def substitute_paths(sql, tmp_dir):
    """Substitute runtime path placeholders in SQL strings."""
    return sql.replace("duckdb_unittest_tempdir", tmp_dir)


def format_value(v):
    """Format a result value to match sqllogictest expected output."""
    if v is None:
        return "NULL"
    return str(v)


def run_test_file(path, ext_dir, tmp_dir):
    """Parse and run a single .test file. Returns (passed, failed, errors)."""
    passed = 0
    failed = 0
    errors = []

    con = duckdb.connect()
    # Make the temp directory path available to tests that use duckdb_unittest_tempdir
    con.execute(f"SET temp_directory='{tmp_dir}'")

    with open(path, encoding="utf-8") as f:
        lines = f.readlines()

    i = 0
    total = len(lines)

    while i < total:
        line = lines[i].rstrip("\n")
        i += 1

        # Skip blank lines and comments
        if not line.strip() or line.startswith("#"):
            continue

        directive_lineno = i  # after increment, i is the 1-based line number of the directive

        # --- require ---
        if line.startswith("require "):
            ext = line.split()[1]
            try:
                load_extension(con, ext, ext_dir)
            except Exception as e:
                errors.append(f"{path}:{directive_lineno}: require {ext} failed: {e}")
                failed += 1
            continue

        # --- statement ok ---
        if line == "statement ok":
            sql_lines = []
            while i < total:
                sql_line = lines[i].rstrip("\n")
                i += 1
                if not sql_line.strip():
                    break
                sql_lines.append(sql_line)
            sql = " ".join(sql_lines)
            # Substitute the unittest tempdir placeholder
            sql = substitute_paths(sql, tmp_dir)
            try:
                con.execute(sql)
                passed += 1
            except Exception as e:
                errors.append(
                    f"{path}:{directive_lineno}: statement ok failed:\n"
                    f"  SQL: {sql}\n"
                    f"  Error: {e}"
                )
                failed += 1
            continue

        # --- query [types] ---
        if line.startswith("query "):
            sql_lines = []
            while i < total:
                sql_line = lines[i].rstrip("\n")
                i += 1
                if sql_line.startswith("----"):
                    break
                sql_lines.append(sql_line)
            sql = " ".join(sql_lines)
            sql = substitute_paths(sql, tmp_dir)

            # Collect expected result lines (until blank line or EOF)
            expected_rows = []
            while i < total:
                result_line = lines[i].rstrip("\n")
                i += 1
                if not result_line.strip():
                    break
                expected_rows.append(result_line)

            try:
                result = con.execute(sql).fetchall()
                # Build actual result rows as tab-separated strings
                actual_rows = ["\t".join(format_value(v) for v in row) for row in result]

                if actual_rows == expected_rows:
                    passed += 1
                else:
                    errors.append(
                        f"{path}:{directive_lineno}: query result mismatch:\n"
                        f"  SQL: {sql}\n"
                        f"  Expected: {expected_rows}\n"
                        f"  Got:      {actual_rows}"
                    )
                    failed += 1
            except Exception as e:
                errors.append(
                    f"{path}:{directive_lineno}: query failed:\n"
                    f"  SQL: {sql}\n"
                    f"  Error: {e}"
                )
                failed += 1
            continue

    con.close()
    return passed, failed, errors


def main():
    parser = argparse.ArgumentParser(description="Run DuckDB sqllogictest .test files")
    parser.add_argument("--test-dir", default="test/sql", help="Directory containing .test files")
    parser.add_argument(
        "--ext-dir",
        default="build/release",
        help="Directory containing the built extension (default: build/release)",
    )
    args = parser.parse_args()

    test_files = sorted(glob.glob(os.path.join(args.test_dir, "*.test")))
    if not test_files:
        print(f"No .test files found in {args.test_dir}", file=sys.stderr)
        sys.exit(1)

    total_passed = 0
    total_failed = 0
    all_errors = []

    with tempfile.TemporaryDirectory() as tmp_dir:
        for test_file in test_files:
            print(f"Running {test_file}...", end=" ", flush=True)
            p, f, errs = run_test_file(test_file, args.ext_dir, tmp_dir)
            total_passed += p
            total_failed += f
            all_errors.extend(errs)
            status = "PASS" if f == 0 else f"FAIL ({f} failures)"
            print(status)

    print(f"\n{total_passed} passed, {total_failed} failed")

    if all_errors:
        print("\nFailures:")
        for err in all_errors:
            print(f"  {err}")
        sys.exit(1)


if __name__ == "__main__":
    main()
