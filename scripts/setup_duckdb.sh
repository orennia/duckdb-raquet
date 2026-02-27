#!/usr/bin/env bash
# setup_duckdb.sh - Download the DuckDB CLI into ./duckdb/
#
# Usage: ./scripts/setup_duckdb.sh [version] [proj_dir]
#
# Populates $(proj_dir)/duckdb/duckdb with the requested CLI binary.
# The duckdb/ directory is gitignored so you can also place a custom
# (pre-release or self-built) binary there manually.

set -euo pipefail

TARGET_DUCKDB_VERSION="${1:-v1.2.0}"
PROJ_DIR="${2:-$(cd "$(dirname "$0")/.." && pwd)}"
DUCKDB_DIR="${PROJ_DIR}/duckdb"
DUCKDB_CLI="${DUCKDB_DIR}/duckdb"

if [[ -x "${DUCKDB_CLI}" ]]; then
    echo "DuckDB already present: $(${DUCKDB_CLI} --version 2>&1 | head -1)"
    echo "  Path: ${DUCKDB_CLI}"
    echo "To replace it, remove ./duckdb/duckdb and re-run 'make setup'."
    exit 0
fi

OS="$(uname -s | tr '[:upper:]' '[:lower:]')"
ARCH="$(uname -m)"
case "${OS}-${ARCH}" in
    linux-x86_64)   ZIP_NAME="duckdb_cli-linux-amd64.zip" ;;
    linux-aarch64)  ZIP_NAME="duckdb_cli-linux-aarch64.zip" ;;
    darwin-*)       ZIP_NAME="duckdb_cli-osx-universal.zip" ;;
    *)
        echo "ERROR: Unsupported platform '${OS}-${ARCH}'."
        echo "Download DuckDB manually from:"
        echo "  https://github.com/duckdb/duckdb/releases/tag/${TARGET_DUCKDB_VERSION}"
        echo "and place the binary at: ${DUCKDB_CLI}"
        exit 1 ;;
esac

DOWNLOAD_URL="https://github.com/duckdb/duckdb/releases/download/${TARGET_DUCKDB_VERSION}/${ZIP_NAME}"
echo "Downloading DuckDB ${TARGET_DUCKDB_VERSION} (${ZIP_NAME})..."

mkdir -p "${DUCKDB_DIR}"
curl -sSL "${DOWNLOAD_URL}" -o "${DUCKDB_DIR}/duckdb.zip"
unzip -qo "${DUCKDB_DIR}/duckdb.zip" -d "${DUCKDB_DIR}"
rm -f "${DUCKDB_DIR}/duckdb.zip"
chmod +x "${DUCKDB_CLI}"

echo "Installed: $(${DUCKDB_CLI} --version 2>&1 | head -1)"
echo "  Path: ${DUCKDB_CLI}"
echo ""
echo "To use a pre-release or custom build, replace ${DUCKDB_CLI} with your binary."
