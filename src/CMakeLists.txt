set(CPUINFO_PATH "/proc/cpuinfo" CACHE STRING "Path to the CPU info file on the system")
set(TARGET_PLATFORM "" CACHE STRING "Target platform for the converter tool. Leave empty to autodetect CPU. Use converter --help to list possible platforms.")
set(TARGET_DEVICE "" CACHE STRING "Target device on the selected platform. Leave empty to use the default.")
set(SCENE_FILE "${PROJECT_SOURCE_DIR}/bench_scenes/bedroom/bedroom.obj" CACHE FILEPATH "Absolute path to an OBJ scene")
set(MEGAKERNEL_FUSION ON CACHE BOOLEAN "Set to true to enable simple shader fusion for megakernel mappings")
set(MAX_PATH_LEN "64" CACHE STRING "Maximum path length")
set(DISABLE_GUI TRUE CACHE BOOLEAN "Set to true to disable GUI")
set(SPP "4" CACHE STRING "Samples per pixel")
if (SCENE_FILE STREQUAL "")
    message(FATAL_ERROR "Please specify a valid OBJ scene in the SCENE_FILE variable")
endif()
set(CONVERTER_OPTIONS "")
set(CONVERTER_OPTIONS "--dev-num" "1")
if (NOT TARGET_PLATFORM STREQUAL "")
    set(CONVERTER_OPTIONS ${CONVERTER_OPTIONS} "--target" "${TARGET_PLATFORM}")
endif()
if (NOT TARGET_DEVICE STREQUAL "")
    set(CONVERTER_OPTIONS ${CONVERTER_OPTIONS} "--device" "${TARGET_DEVICE}")
    if (MEGAKERNEL_FUSION AND (TARGET_PLATFORM STREQUAL "nvvm-megakernel" OR TARGET_PLATFORM STREQUAL "amdgpu-megakernel"))
        set(CONVERTER_OPTIONS ${CONVERTER_OPTIONS} "--fusion")
    endif()
endif()
set(RODENT_SRCS
    core/color.impala
    core/common.impala
    core/cpu_common.impala
    core/matrix.impala
    core/random.impala
    core/sort.impala
    core/vector.impala
    core/file.impala
    render/chunk_hit.impala
    render/image.impala
    render/camera.impala
    render/geometry.impala
    render/light.impala
    render/material.impala
    render/sdtree.impala
    render/renderer.impala
    render/scene.impala
    render/driver.impala
    render/mapping_cpu.impala
    render/mapping_gpu.impala
    traversal/intersection.impala
    traversal/stack.impala
    traversal/mapping_cpu.impala
    traversal/mapping_gpu.impala
)

set(DRIVER_SRCS
    driver/interface.h
    driver/interface.cpp
    distributed/Node.h 
    distributed/ProcStatus.h
    driver/driver.cpp
    driver/obj.cpp
    driver/obj.h
    driver/image.cpp
    driver/image.h
    driver/bvh.h
    driver/float2.h
    driver/float3.h
    driver/float4.h
    driver/common.h
    driver/color.h
    driver/buffer.h
)

set(CONVERTER_SRCS
    driver/converter.cpp
    driver/obj.cpp
    driver/obj.h
    driver/interface.h
    driver/bvh.h)



find_package(MPI REQUIRED)
if (MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
else (MPI_FOUND)
    message(SEND_ERROR "This application cannot compile without MPI")
endif (MPI_FOUND)

anydsl_runtime_wrap(RODENT_OBJS
    NAME "rodent"
    CLANG_FLAGS ${CLANG_FLAGS}
    IMPALA_FLAGS --log-level info
    FILES ${RODENT_SRCS} ${CMAKE_CURRENT_BINARY_DIR}/main.impala)

anydsl_runtime_wrap(DISCARD_TMP_OBJS
    NAME "interface"
    FILES ${RODENT_SRCS} dummy_main.impala
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/driver/interface)


if (NOT DISABLE_GUI)
    find_package(SDL2 REQUIRED)
endif()
find_package(PNG REQUIRED)
find_package(JPEG REQUIRED)
find_package(TBB REQUIRED)
find_package(LZ4 REQUIRED)

add_executable(converter ${CONVERTER_SRCS})
target_include_directories(converter PUBLIC ${LZ4_INCLUDE_DIR} )
target_link_libraries(converter ${LZ4_LIBRARY} ${MPI_LIBRARIES} )
target_compile_definitions(converter PUBLIC -DCPUINFO_PATH="${CPUINFO_PATH}")
target_compile_features(converter PRIVATE cxx_std_14)
if (COLORIZE)
    target_compile_definitions(converter PUBLIC -DCOLORIZE)
endif()
if (EMBREE_FOUND)
    target_include_directories(converter PUBLIC ${EMBREE_ROOT_DIR}/include ${EMBREE_ROOT_DIR} ${EMBREE_LIBRARY_DIR})
    target_link_libraries(converter ${EMBREE_DEPENDENCIES})
    target_compile_definitions(converter PUBLIC ${EMBREE_DEFINITIONS} -DENABLE_EMBREE_BVH)
endif()

#add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/main.impala
#    COMMAND converter ${SCENE_FILE} ${CONVERTER_OPTIONS} --max-path-len ${MAX_PATH_LEN} --samples-per-pixel ${SPP}
#    COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_BINARY_DIR}/main.impala ${CMAKE_CURRENT_BINARY_DIR}/main.impala
#    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#    DEPENDS ${SCENE_FILE} converter)

set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/main.impala PROPERTIES GENERATED TRUE)

add_library(driver ${DRIVER_SRCS})
target_include_directories(driver PUBLIC ${LZ4_INCLUDE_DIR})
if (EMBREE_FOUND)
    target_include_directories(driver PUBLIC ${EMBREE_ROOT_DIR}/include ${EMBREE_ROOT_DIR} ${EMBREE_LIBRARY_DIR} ${PNG_INCLUDE_DIRS} ${JPEG_INCLUDE_DIRS} ${TBB_INCLUDE_DIRS} ${LZ4_INCLUDE_DIR})
    target_link_libraries(driver ${EMBREE_DEPENDENCIES})
    target_compile_definitions(driver PUBLIC -DENABLE_EMBREE_DEVICE)
endif()
if (NOT DISABLE_GUI)
    target_include_directories(driver PUBLIC ${SDL2_INCLUDE_DIR})
    target_link_libraries(driver ${SDL2_LIBRARY})
else()
    target_compile_definitions(driver PUBLIC -DDISABLE_GUI)
endif()
if (COLORIZE)
    target_compile_definitions(driver PUBLIC -DCOLORIZE)
endif()

add_executable(rodent ${RODENT_OBJS})
target_link_libraries(rodent driver ${AnyDSL_runtime_LIBRARIES} ${PNG_LIBRARIES} ${JPEG_LIBRARIES} ${TBB_LIBRARIES} ${LZ4_LIBRARY} ${MPI_LIBRARIES} -lpthread)
