// tests/main_tests.cpp
#define DOCTEST_CONFIG_IMPLEMENT
#include "../third_party/doctest.h"

// This file only defines the test runner.
// All tests in any test_*.cpp will be automatically registered.
int main(int argc, char** argv) {
    doctest::Context context;

    context.applyCommandLine(argc, argv);
    int res = context.run();   // run all registered tests

    // important — propagate test failures
    if (context.shouldExit())
        return res;

    return res;
}