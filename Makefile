CXX      := g++
CXXFLAGS := -std=c++20 -I third_party -I Multiview

TEST_BIN := unit_tests
TEST_SRC := tests/main_tests.cpp tests/dummy_test.cpp

# We should try as much as possible to separate Rcpp code from the rest,
# right now this is more of a "placeholder" solution.
MV_SRC := \
	Multiview/multiview_hyper.cpp \
	Multiview/multiview_state.cpp \
	#Multiview/multiview_gibbs.cpp \
	#Multiview/multiview_utils.cpp

.PHONY: all test clean

all: test

test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): $(TEST_SRC) $(MV_SRC)
	$(CXX) $(CXXFLAGS) $(TEST_SRC) $(MV_SRC) -o $(TEST_BIN)

clean:
	rm -f $(TEST_BIN)
