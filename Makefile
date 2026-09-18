.DEFAULT_GOAL := release

DIR_SRC = ./src
DIR_OBJ = ./obj
DIR_BIN = ./bin
DIR_TEST = ./tests
HTSLIB_DIR = ./deps/htslib
BINDIR=/usr/local/bin

SRC = $(wildcard ${DIR_SRC}/*.cpp)  
OBJ = $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC})) 
DEP = $(OBJ:.o=.d)
TEST_SRC = $(wildcard ${DIR_TEST}/*.cpp)
TEST_OBJ = $(patsubst %.cpp,${DIR_OBJ}/test_%.o,$(notdir ${TEST_SRC}))
TEST_DEP = $(TEST_OBJ:.o=.d)
TEST_LINK_OBJ = $(filter-out ${DIR_OBJ}/main.o,$(OBJ))
RELEASE_OBJ = $(patsubst %.cpp,${DIR_OBJ}/release/%.o,$(notdir ${SRC}))
RELEASE_DEP = $(RELEASE_OBJ:.o=.d)

TARGET = gencore

BIN_TARGET = ${DIR_BIN}/${TARGET}
TEST_TARGET = ${DIR_BIN}/${TARGET}_tests

CXX = g++
CPPFLAGS = -I${DIR_SRC} -I${HTSLIB_DIR}/include
CXXFLAGS = -std=c++20 -O3 -g -MMD -MP
RELEASE_CXXFLAGS = -std=c++20 -O3 -DNDEBUG -MMD -MP
HTSLIB_LIB = ${HTSLIB_DIR}/lib/libhts.a
LDLIBS = ${HTSLIB_LIB} -ldeflate -llzma -lbz2 -lz -lm -lpthread

# Statically link the complete application, including C/C++ runtimes and
# compression libraries. Keep test linking separate for sanitizer support.
${BIN_TARGET}:${OBJ} ${HTSLIB_LIB} Makefile | make_bin_dir
	$(CXX) $(LDFLAGS) -static $(OBJ) $(LDLIBS) -o $@
    
${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp | make_obj_dir
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

${DIR_OBJ}/test_%.o:${DIR_TEST}/%.cpp | make_obj_dir
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Separate objects prevent reusing a debug/sanitizer build for distribution.
${DIR_OBJ}/release/%.o:${DIR_SRC}/%.cpp Makefile | make_release_obj_dir
	$(CXX) $(CPPFLAGS) $(RELEASE_CXXFLAGS) -c $< -o $@

# -s also removes debug symbols embedded in third-party static archives.
release: $(RELEASE_OBJ) $(HTSLIB_LIB) | make_bin_dir
	$(CXX) $(LDFLAGS) -static -s $(RELEASE_OBJ) $(LDLIBS) -o $(BIN_TARGET)

${TEST_TARGET}:${TEST_LINK_OBJ} ${TEST_OBJ} | make_bin_dir
	$(CXX) $(LDFLAGS) $^ -lgtest_main -lgtest $(LDLIBS) -o $@

.PHONY: clean make_obj_dir make_release_obj_dir make_bin_dir install test release
clean:
	rm -f $(OBJ) $(DEP) $(TEST_OBJ) $(TEST_DEP) $(RELEASE_OBJ) $(RELEASE_DEP) $(BIN_TARGET) $(TEST_TARGET) $(TARGET)

make_obj_dir:
	mkdir -p $(DIR_OBJ)

make_release_obj_dir:
	mkdir -p $(DIR_OBJ)/release

make_bin_dir:
	mkdir -p $(DIR_BIN)

install: $(BIN_TARGET)
	install $(BIN_TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."

test: $(TEST_TARGET) $(BIN_TARGET)
	$(TEST_TARGET)
	python3 $(DIR_TEST)/test_umi_pipeline.py

-include $(DEP) $(TEST_DEP) $(RELEASE_DEP)
