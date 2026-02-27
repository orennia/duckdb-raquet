#!/usr/bin/env python3
"""append_metadata.py - Append (or replace) the DuckDB extension metadata footer.

Usage: python3 scripts/append_metadata.py <extension.so> [--platform PLATFORM]
                                           [--duckdb-version VERSION]
                                           [--ext-version VERSION]

DuckDB validates a 512-byte footer at the end of every extension file:

  bytes  0..256  8 metadata fields × 32 bytes (written in reverse: field[7]..field[0])
  bytes 256..512 256-byte RSA signature (all zeros for unsigned extensions)

Fields (after DuckDB reverses them on read):
  field[0] = magic value      "4"
  field[1] = platform         e.g. "linux_amd64"
  field[2] = duckdb_capi_ver  e.g. "v1.2.0"
  field[3] = ext_version      e.g. "0.5.0"
  field[4] = abi_type         "C_STRUCT"
  field[5..7] = empty
"""

import argparse
import platform as _platform
import struct
import sys
from pathlib import Path

FIELD_WIDTH   = 32
NUM_FIELDS    = 8
METADATA_SIZE = NUM_FIELDS * FIELD_WIDTH   # 256 bytes
SIGNATURE_SIZE = 256                        # 256 bytes (zeros for unsigned)
FOOTER_SIZE   = METADATA_SIZE + SIGNATURE_SIZE  # 512 bytes


def make_field(value: str) -> bytes:
    """Return value encoded as UTF-8 and zero-padded / truncated to FIELD_WIDTH bytes."""
    b = value.encode("utf-8")[:FIELD_WIDTH]
    return b.ljust(FIELD_WIDTH, b"\x00")


def detect_platform() -> str:
    machine = _platform.machine().lower()
    arch = "amd64" if machine in ("x86_64", "amd64") else "arm64"
    system = _platform.system().lower()
    if system == "linux":
        # DuckDB uses the gcc4 ABI suffix on Linux for pre-built binaries.
        return f"linux_{arch}_gcc4"
    if system == "darwin":
        return f"osx_{arch}"
    if system == "windows":
        return f"windows_{arch}"
    return f"unknown_{arch}"


def looks_like_footer(data: bytes) -> bool:
    """Return True if the last FIELD_WIDTH bytes of the metadata block == b'4\x00...'"""
    if len(data) < FOOTER_SIZE:
        return False
    # Last field (METADATA1 = magic "4") is at bytes[METADATA_SIZE-FIELD_WIDTH:METADATA_SIZE]
    magic = data[METADATA_SIZE - FIELD_WIDTH : METADATA_SIZE]
    return magic[:1] == b"4" and magic[1:] == b"\x00" * (FIELD_WIDTH - 1)


def build_footer(platform: str, duckdb_version: str, ext_version: str) -> bytes:
    """Return the 512-byte footer (metadata + signature)."""
    # Fields written in reverse order (METADATA8 first → METADATA1 last):
    fields = [
        make_field(""),              # field[7] empty
        make_field(""),              # field[6] empty
        make_field(""),              # field[5] empty
        make_field("C_STRUCT"),      # field[4] ABI type
        make_field(ext_version),     # field[3] extension version
        make_field(duckdb_version),  # field[2] DuckDB C API version
        make_field(platform),        # field[1] platform
        make_field("4"),             # field[0] magic value
    ]
    metadata = b"".join(fields)
    assert len(metadata) == METADATA_SIZE, f"metadata length {len(metadata)} != {METADATA_SIZE}"
    signature = b"\x00" * SIGNATURE_SIZE
    return metadata + signature


def append_metadata(path: Path, platform: str, duckdb_version: str, ext_version: str) -> None:
    data = path.read_bytes()

    # Strip any existing footer so the operation is idempotent.
    if len(data) >= FOOTER_SIZE and looks_like_footer(data[-FOOTER_SIZE:]):
        data = data[:-FOOTER_SIZE]

    footer = build_footer(platform, duckdb_version, ext_version)
    path.write_bytes(data + footer)
    print(f"  Metadata appended to {path.name} "
          f"(platform={platform}, duckdb={duckdb_version}, ext={ext_version})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("extension", type=Path, help="Path to the extension .so/.dylib/.dll")
    parser.add_argument("--platform",       default=detect_platform(),
                        help="DuckDB platform string (default: auto-detect)")
    parser.add_argument("--duckdb-version", default="v1.2.0",
                        help="Minimum DuckDB C API version (default: v1.2.0)")
    parser.add_argument("--ext-version",    default="0.0.0",
                        help="Extension version string (default: 0.0.0)")
    args = parser.parse_args()

    if not args.extension.exists():
        sys.exit(f"ERROR: Extension file not found: {args.extension}")

    append_metadata(args.extension, args.platform, args.duckdb_version, args.ext_version)


if __name__ == "__main__":
    main()
