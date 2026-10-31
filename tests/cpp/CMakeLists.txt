# C++ tests for NNS-core.
#
# Dependency-free test binaries (see check.hpp); each exits non-zero on
# failure so it can be driven by ctest.

set(NNSCORE_CPP_TESTS
    smoke_partial_moments
    test_prefix_backend
    test_degree_validation
    test_central_tendencies
)

foreach(test_name IN LISTS NNSCORE_CPP_TESTS)
    add_executable(${test_name} ${test_name}.cpp)
    target_link_libraries(${test_name} PRIVATE nnscore)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
