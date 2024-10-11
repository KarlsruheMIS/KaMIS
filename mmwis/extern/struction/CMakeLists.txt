

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fno-omit-frame-pointer -march=native -DNDEBUG -fPIC -O3 ")

set(STRUCTION_INCLUDES
        app
        lib
        lib/tools
        lib/mis
        lib/mis/kernel
)

set(STRUCTION_SOURCES
        lib/tools/mis_log.cpp
	    $<TARGET_OBJECTS:libkaffpa> 
        $<TARGET_OBJECTS:lib_mmwis>
)

add_executable(branch_reduce ${LIB_MWIS} ${STRUCTION_SOURCES} "app/branch_reduce.cpp")
add_executable(branch_reduce_convergence ${LIB_MWIS} ${STRUCTION_SOURCES} "app/branch_reduce.cpp")
add_executable(cyclic_blow_up_convergency ${LIB_MWIS} ${STRUCTION_SOURCES} "app/cyclic_blow_up_convergency.cpp")
add_executable(kernelization ${LIB_MWIS} ${STRUCTION_SOURCES} "app/kernelization.cpp")
add_executable(struction_checker ${LIB_MWIS} ${STRUCTION_SOURCES} "app/struction_checker.cpp")

target_include_directories(branch_reduce PRIVATE ${MMWIS_INCLUDES})
target_include_directories(branch_reduce PRIVATE ${STRUCTION_INCLUDES})
target_include_directories(branch_reduce_convergence PRIVATE ${MMWIS_INCLUDES})
target_include_directories(branch_reduce_convergence PRIVATE ${STRUCTION_INCLUDES})
target_include_directories(cyclic_blow_up_convergency PRIVATE ${MMWIS_INCLUDES})
target_include_directories(cyclic_blow_up_convergency PRIVATE ${STRUCTION_INCLUDES})
target_include_directories(struction_checker PRIVATE ${MMWIS_INCLUDES})
target_include_directories(struction_checker PRIVATE ${STRUCTION_INCLUDES})
target_include_directories(kernelization PRIVATE ${MMWIS_INCLUDES})
target_include_directories(kernelization PRIVATE ${STRUCTION_INCLUDES})
target_compile_definitions(cyclic_blow_up_convergency PRIVATE -D OUTPUT_GRAPH_CONVERGENCE)
target_compile_definitions(branch_reduce_convergence PRIVATE -D OUTPUT_WEIGHT_CONVERGENCE)
target_link_libraries(branch_reduce LINK_PUBLIC ${LIBS})
target_link_libraries(branch_reduce_convergence LINK_PUBLIC ${LIBS})
target_link_libraries(cyclic_blow_up_convergency LINK_PUBLIC ${LIBS})
target_link_libraries(kernelization LINK_PUBLIC ${LIBS})
target_link_libraries(struction_checker LINK_PUBLIC ${LIBS})

