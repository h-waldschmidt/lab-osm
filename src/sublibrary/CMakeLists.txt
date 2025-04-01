set(HEADER_FILES sublib_component.h)
set(SOURCE_FILES sublib_component.cpp)

add_library(sublibrary SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(sublibrary PUBLIC fmt::fmt)

if(BUILD_TESTING)
    set(TEST_FILES sublib_component_test.cpp)
    add_library(sublibrary_test SHARED ${TEST_FILES})
    target_link_libraries(sublibrary_test PRIVATE sublibrary)
    target_link_libraries(sublibrary_test PRIVATE GTest::gtest GTest::gtest_main GTest::gmock GTest::gmock_main)

    append_coverage_compiler_flags()
endif()