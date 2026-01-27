# DuckDB Community Extensions Submission

This directory contains the files needed to submit the raquet extension to the [DuckDB Community Extensions](https://github.com/duckdb/community-extensions) repository.

## Submission Steps

1. **Fork the community-extensions repo:**
   ```bash
   git clone https://github.com/duckdb/community-extensions.git
   cd community-extensions
   ```

2. **Create the extension directory:**
   ```bash
   mkdir -p extensions/raquet
   cp /path/to/duckdb-raquet/community-extensions/description.yml extensions/raquet/
   ```

3. **Submit Pull Request:**
   ```bash
   git checkout -b add-raquet-extension
   git add extensions/raquet/description.yml
   git commit -m "Add raquet extension - Raster analytics with QUADBIN indexing"
   git push origin add-raquet-extension
   ```

4. **Open PR** at https://github.com/duckdb/community-extensions/pulls

## Requirements

- Extension must build with DuckDB's CI toolchain
- Currently requires DuckDB development branch (main) for `LogicalType::GEOMETRY()` support
- Will be compatible with DuckDB 1.5+ once GEOMETRY type is in stable release

## Testing Before Submission

```bash
# Build and test locally
cd /path/to/duckdb-raquet
make release
./build/release/test/unittest "*raquet*"
```

## Status

- [ ] Waiting for DuckDB 1.5 release (GEOMETRY type in core)
- [ ] Submit PR to community-extensions
- [ ] Extension available via `INSTALL raquet FROM community`
