# Live Rcpp fidelity helper.
#
# Builds the `nns_eval` CLI used by scripts/run_fidelity.py to compare the
# portable C++ core against truth generated from the installed R NNS package.
# This is intentionally NOT registered as a ctest: it requires a truth file
# produced from live R, which is generated in the sync workflow.

add_executable(nns_eval nns_eval.cpp)
target_link_libraries(nns_eval PRIVATE nnscore)
