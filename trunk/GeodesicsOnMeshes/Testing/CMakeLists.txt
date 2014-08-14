project( vtkMeshGeodesicTests )

find_package( VTK REQUIRED )
include( ${VTK_USE_FILE} )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../Source
                     ${CMAKE_CURRENT_SOURCE_DIR}  )

# easier regression testing
include(${CMAKE_CURRENT_SOURCE_DIR}/vtkTestingObjectFactory2.cmake)

CREATE_TEST_SOURCELIST( vtkMeshGeodesicTests
  vtkMeshGeodesicTests.cxx
  TestGeodesicDistance.cxx
  TestGeodesicPath.cxx
  EXTRA_INCLUDE vtkMeshGeodesicTestingIncludes.h # needed when building with VTK58
  )

add_executable( vtkMeshGeodesicTests ${vtkMeshGeodesicTests} )
target_link_libraries( vtkMeshGeodesicTests vtkMeshGeodesics ${VTK_LIBRARIES} )

ADD_TEST( TestGeodesicDistance ${EXECUTABLE_OUTPUT_PATH}/vtkMeshGeodesicTests
  TestGeodesicDistance
  Bunny.vtp BunnyGeodesic.vtp
  -D ${CMAKE_CURRENT_SOURCE_DIR}/../Data
  -B ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines
  -T ${CMAKE_CURRENT_BINARY_DIR}
  -V ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines/TestGeodesicDistance.png )

ADD_TEST( TestGeodesicDistanceMaxDistTermCriteria ${EXECUTABLE_OUTPUT_PATH}/vtkMeshGeodesicTests
  TestGeodesicDistance
  Bunny.vtp BunnyGeodesicMaxDist.vtp --maxDist 0.7
  -D ${CMAKE_CURRENT_SOURCE_DIR}/../Data
  -B ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines
  -T ${CMAKE_CURRENT_BINARY_DIR}
  -V ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines/TestGeodesicDistanceMaxDistTermCriteria.png )

ADD_TEST( TestGeodesicDistanceExclusionRegion ${EXECUTABLE_OUTPUT_PATH}/vtkMeshGeodesicTests
  TestGeodesicDistance
  Bunny.vtp BunnyGeodesicExclusionContour.vtp --exclusionContour
  -D ${CMAKE_CURRENT_SOURCE_DIR}/../Data
  -B ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines
  -T ${CMAKE_CURRENT_BINARY_DIR}
  -V ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines/TestGeodesicDistanceExclusionRegion.png )

ADD_TEST( TestGeodesicDistancePropagationWt ${EXECUTABLE_OUTPUT_PATH}/vtkMeshGeodesicTests
  TestGeodesicDistance
  Ellipsoid.vtp EllipsoidCurvBasedPropagation.vtp --propagationWts CurvWt
  -D ${CMAKE_CURRENT_SOURCE_DIR}/../Data
  -B ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines
  -T ${CMAKE_CURRENT_BINARY_DIR}
  -V ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines/TestGeodesicDistancePropagationWt.png )

ADD_TEST( TestGeodesicDistanceNoPropagationWt ${EXECUTABLE_OUTPUT_PATH}/vtkMeshGeodesicTests
  TestGeodesicDistance
  Ellipsoid.vtp EllipsoidNoCurvBasedPropagation.vtp --propagationWts None
  -D ${CMAKE_CURRENT_SOURCE_DIR}/../Data
  -B ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines
  -T ${CMAKE_CURRENT_BINARY_DIR}
  -V ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines/TestGeodesicDistanceNoPropagationWt.png )

ADD_TEST( TestGeodesicPathFastMarchingLinear ${EXECUTABLE_OUTPUT_PATH}/vtkMeshGeodesicTests
  TestGeodesicPath
  Bunny.vtp 1 1 0
  -D ${CMAKE_CURRENT_SOURCE_DIR}/../Data
  -B ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines
  -T ${CMAKE_CURRENT_BINARY_DIR} 
  -V ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines/TestGeodesicPathFastMarchingLinear.png )

ADD_TEST( TestGeodesicPathFastMarchingNearest ${EXECUTABLE_OUTPUT_PATH}/vtkMeshGeodesicTests
  TestGeodesicPath
  Bunny.vtp 1 0 0
  -D ${CMAKE_CURRENT_SOURCE_DIR}/../Data
  -B ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines
  -T ${CMAKE_CURRENT_BINARY_DIR} 
  -V ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines/TestGeodesicPathFastMarchingNearest.png )

ADD_TEST( TestGeodesicPathDijkstraNearest ${EXECUTABLE_OUTPUT_PATH}/vtkMeshGeodesicTests
  TestGeodesicPath
  Bunny.vtp 0 0 0
  -D ${CMAKE_CURRENT_SOURCE_DIR}/../Data
  -B ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines
  -T ${CMAKE_CURRENT_BINARY_DIR} 
  -V ${CMAKE_CURRENT_SOURCE_DIR}/../TestBaselines/TestGeodesicPathDijkstraNearest.png )
