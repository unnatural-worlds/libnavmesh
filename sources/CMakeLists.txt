
file(GLOB_RECURSE sources "navmesh/*" "include/**/*")
add_library(unnatural-navmesh STATIC ${sources})
target_include_directories(unnatural-navmesh PUBLIC include)
target_link_libraries(unnatural-navmesh PUBLIC cage-core)
cage_ide_category(unnatural-navmesh unnatural/libs)
cage_ide_sort_files(unnatural-navmesh)
