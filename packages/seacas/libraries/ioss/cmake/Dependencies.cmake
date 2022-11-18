if (CMAKE_PROJECT_NAME STREQUAL "Seacas")
TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_OPTIONAL_PACKAGES SEACASExodus Zoltan
  LIB_OPTIONAL_TPLS HDF5 Pamgen CGNS ParMETIS Faodel Cereal DLlib Pthread ADIOS2 Catalyst2 ${SEACAS_GTest_TPL_name} Kokkos DataWarp fmt
)
else()
TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_OPTIONAL_PACKAGES SEACASExodus Pamgen Zoltan Kokkos
  LIB_OPTIONAL_TPLS HDF5 CGNS ParMETIS Faodel Cereal DLlib Pthread DataWarp ADIOS2 Catalyst2 ${SEACAS_GTest_TPL_name}
)
endif()

TRIBITS_TPL_TENTATIVELY_ENABLE(DLlib)
