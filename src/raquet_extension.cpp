#define DUCKDB_EXTENSION_MAIN

#include "raquet_extension.hpp"
#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// Forward declarations for function registration
void RegisterQuadbinFunctions(ExtensionLoader &loader);
void RegisterRasterValueFunctions(ExtensionLoader &loader);
void RegisterRasterStatsFunctions(ExtensionLoader &loader);
void RegisterMetadataFunctions(ExtensionLoader &loader);
void RegisterRaquetTableFunctions(ExtensionLoader &loader);

static void LoadInternal(ExtensionLoader &loader) {
    // Register all functions
    RegisterQuadbinFunctions(loader);
    RegisterRasterValueFunctions(loader);
    RegisterRasterStatsFunctions(loader);
    RegisterMetadataFunctions(loader);
    RegisterRaquetTableFunctions(loader);
}

void RaquetExtension::Load(ExtensionLoader &loader) {
    LoadInternal(loader);
}

std::string RaquetExtension::Name() {
    return "raquet";
}

std::string RaquetExtension::Version() const {
#ifdef EXT_VERSION_RAQUET
    return EXT_VERSION_RAQUET;
#else
    return "0.1.0";
#endif
}

} // namespace duckdb

extern "C" {

DUCKDB_CPP_EXTENSION_ENTRY(raquet, loader) {
    duckdb::LoadInternal(loader);
}

}

#ifndef DUCKDB_EXTENSION_MAIN
#error DUCKDB_EXTENSION_MAIN not defined
#endif
