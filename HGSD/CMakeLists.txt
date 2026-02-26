cmake_minimum_required(VERSION 3.20)
project(HGSD CXX)
set(CMAKE_CXX_STANDARD 17)

message(STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(deps/geometry-central)
add_subdirectory(deps/polyscope)
include_directories( "/opt/homebrew/include/eigen3" )
include_directories("deps/fcpw-master/include")
add_executable(main
        ${CMAKE_SOURCE_DIR}/src/ZeroSignedDistance.cpp
#        ${CMAKE_SOURCE_DIR}/src/test.cpp
)
target_include_directories(main PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include)
target_link_libraries(main PRIVATE geometry-central polyscope)

add_executable(demo
        ${CMAKE_SOURCE_DIR}/src/demoGeoPath.cpp
        include/utils.h
)
target_include_directories(demo PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include)
target_link_libraries(demo PRIVATE geometry-central polyscope)

add_executable(demoapp
        ${CMAKE_SOURCE_DIR}/src/demo_app.cpp
        ${CMAKE_SOURCE_DIR}/src/utils.cpp
)
target_include_directories(demoapp PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include include)
target_link_libraries(demoapp PRIVATE geometry-central polyscope)


# draw some demo
add_executable(quicktest
        src/quicktest.cpp
        ${CMAKE_SOURCE_DIR}/src/utils.cpp
)
target_include_directories(quicktest PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include include)
target_link_libraries(quicktest PRIVATE geometry-central polyscope)

# draw some demo
add_executable(demoPatchTransfer
        src/demoPatchTransfer.cpp
        ${CMAKE_SOURCE_DIR}/src/utils.cpp
)
target_include_directories(demoPatchTransfer PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include include)
target_link_libraries(demoPatchTransfer PRIVATE geometry-central polyscope)

# demo QSLIM propcess
add_executable(demoQSLIM
        src/demoQSLIM.cpp
        ${CMAKE_SOURCE_DIR}/src/utils.cpp
)
target_include_directories(demoQSLIM PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include include)
target_link_libraries(demoQSLIM PRIVATE geometry-central polyscope)

# demo QSLIM RealFakeHole
add_executable(demoRealFakeHole
        src/demoRealFakeHole.cpp
        ${CMAKE_SOURCE_DIR}/src/utils.cpp
)
target_include_directories(demoRealFakeHole PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include include)
target_link_libraries(demoRealFakeHole PRIVATE geometry-central polyscope)

add_executable(isoTopoisson
        src/isoTopoisson.cpp
)
target_include_directories(isoTopoisson PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include)
target_link_libraries(isoTopoisson PRIVATE geometry-central polyscope)


add_executable(flipGeodesicsDemo
        src/flipGeodesicsDemo.cpp
)
target_include_directories(flipGeodesicsDemo PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include)
target_link_libraries(flipGeodesicsDemo PRIVATE geometry-central polyscope)

add_executable(remeshProj
        ${CMAKE_SOURCE_DIR}/src/remeshProj.cpp

)
target_include_directories(remeshProj PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include )
target_link_libraries(remeshProj PRIVATE geometry-central polyscope)


add_executable(MeshImplanting
        ${CMAKE_SOURCE_DIR}/src/main.cpp
        ${CMAKE_SOURCE_DIR}/src/utils.cpp
)
target_include_directories(MeshImplanting PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include include)
target_link_libraries(MeshImplanting PRIVATE geometry-central polyscope)

add_executable(cut
        ${CMAKE_SOURCE_DIR}/src/cut.cpp
        ${CMAKE_SOURCE_DIR}/src/utils.cpp
)
target_include_directories(cut PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include include)
target_link_libraries(cut PRIVATE geometry-central polyscope)

add_executable(cutAll
        ${CMAKE_SOURCE_DIR}/src/cutAll.cpp
        ${CMAKE_SOURCE_DIR}/src/utils.cpp
)
target_include_directories(cutAll PUBLIC ${CMAKE_SOURCE_DIR}/src deps/libigl/include include)
target_link_libraries(cutAll PRIVATE geometry-central polyscope)
