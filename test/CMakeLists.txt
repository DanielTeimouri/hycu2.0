if (BUILD_TEST)
  find_package(Boost COMPONENTS
    unit_test_framework
    REQUIRED
    )
  
  include_directories(.)
  aux_source_directory(. HyCu_SOURCES_TEST)

  add_executable(hycu-flint-fq-test
    test.cc
    flint_fq.cc
    flint_fq_polynomial.cc
    ${HyCu_SOURCE_DIR}/src/flint_fq.cc
    ${HyCu_SOURCE_DIR}/src/flint_fq_polynomial.cc
    )
  target_link_libraries(hycu-flint-fq-test
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    ${FLINT_LIBRARY}
    ${GMP_LIBRARY}
    )

  add_executable(hycu-test
    ${HyCu_SOURCES_TEST}
    )
  target_link_libraries(hycu-test
    hycu-test-lib
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    ${FLINT_LIBRARY}
    ${GMP_LIBRARY}
    ${CMAKE_THREAD_LIBS_INIT}
    )
  if (WITH_OPENCL)
    target_link_libraries(hycu-test ${OpenCL_LIBRARY})
  endif()

  add_custom_target(runtest ALL
    COMMAND hycu-test
    DEPENDS hycu-test
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif (BUILD_TEST)
