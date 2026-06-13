# C++ tests for NNS-core.
#
# Currently a dependency-free smoke test exercising the core partial moment
# functions. Catch2-based conformance tests can be added here later.

add_executable(smoke_partial_moments smoke_partial_moments.cpp)
target_link_libraries(smoke_partial_moments PRIVATE nnscore)

add_test(NAME smoke_partial_moments COMMAND smoke_partial_moments)
