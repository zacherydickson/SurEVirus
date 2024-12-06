cmake_minimum_required(VERSION 3.5)
project(SurEVirus)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

include_directories(htslib-1.21/include)
include_directories(sparsehash-sparsehash-2.0.3/include)

set(ISOLATE_RELEVANT_PAIRS_SRC isolate_relevant_pairs.cpp)
add_executable(isolate_relevant_pairs ${ISOLATE_RELEVANT_PAIRS_SRC})

set(ISOLATE_RELEVANT_PAIRS_FQ_SRC isolate_relevant_pairs_fq.cpp)
add_executable(isolate_relevant_pairs_fq ${ISOLATE_RELEVANT_PAIRS_FQ_SRC})

set(EXTRACT_CLIPS_SRC extract_clips.cpp)
add_executable(extract_clips ${EXTRACT_CLIPS_SRC})

set(EXTRACT_REGIONS_SRC extract_regions.cpp)
add_executable(extract_regions ${EXTRACT_REGIONS_SRC})

set(CLUSTER_JUNCTIONS_SRC cluster_junctions.cpp)
add_executable(cluster_junctions ${CLUSTER_JUNCTIONS_SRC})

set(EDGE_MAPPER_SRC libs/ssw.c libs/ssw_cpp.cpp edge_mapper.cpp)
add_executable(edge_mapper ${EDGE_MAPPER_SRC})

set(FILTER_SRC libs/ssw.c libs/ssw_cpp.cpp filter.cpp)
add_executable(filter ${FILTER_SRC})

find_library(HTS_LIB hts PATHS htslib-1.21/lib NO_DEFAULT_PATH)
find_package(ZLIB)
target_link_libraries(cluster_junctions "${HTS_LIB}")
target_link_libraries(isolate_relevant_pairs "${HTS_LIB}")
target_link_libraries(isolate_relevant_pairs_fq "${HTS_LIB}")
target_link_libraries(isolate_relevant_pairs_fq ZLIB::ZLIB)
target_link_libraries(edge_mapper "${HTS_LIB}")
target_link_libraries(extract_clips "${HTS_LIB}")
target_link_libraries(extract_regions "${HTS_LIB}")
target_link_libraries(filter "${HTS_LIB}")
