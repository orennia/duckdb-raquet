#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// raquet_is_metadata_row(block UBIGINT) -> BOOLEAN
// Check if this row is the metadata row (block == 0)
static void RaquetIsMetadataRowFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &block_vec = args.data[0];

    UnaryExecutor::Execute<uint64_t, bool>(block_vec, result, args.size(),
        [](uint64_t block) {
            return block == 0;
        });
}

// raquet_is_data_row(block UBIGINT) -> BOOLEAN
// Check if this row is a data row (block != 0)
static void RaquetIsDataRowFunction(DataChunk &args, ExpressionState &state, Vector &result) {
    auto &block_vec = args.data[0];

    UnaryExecutor::Execute<uint64_t, bool>(block_vec, result, args.size(),
        [](uint64_t block) {
            return block != 0;
        });
}

void RegisterMetadataFunctions(ExtensionLoader &loader) {
    // raquet_is_metadata_row(block) -> BOOLEAN
    ScalarFunction is_metadata("raquet_is_metadata_row",
        {LogicalType::UBIGINT},
        LogicalType::BOOLEAN,
        RaquetIsMetadataRowFunction);
    loader.RegisterFunction(is_metadata);

    // raquet_is_data_row(block) -> BOOLEAN
    ScalarFunction is_data("raquet_is_data_row",
        {LogicalType::UBIGINT},
        LogicalType::BOOLEAN,
        RaquetIsDataRowFunction);
    loader.RegisterFunction(is_data);
}

} // namespace duckdb
