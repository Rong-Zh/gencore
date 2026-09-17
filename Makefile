DIR_SRC = ./src
DIR_OBJ = ./obj
DIR_BIN = ./bin
HTSLIB_DIR = ./deps/htslib
BINDIR=/usr/local/bin

SRC = $(wildcard ${DIR_SRC}/*.cpp)  
OBJ = $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC})) 
DEP = $(OBJ:.o=.d)

TARGET = gencore

BIN_TARGET = ${DIR_BIN}/${TARGET}

CXX = g++
CPPFLAGS = -I${HTSLIB_DIR}/include
CXXFLAGS = -std=c++20 -O3 -g -MMD -MP
HTSLIB_LIB = ${HTSLIB_DIR}/lib/libhts.a
LDLIBS = ${HTSLIB_LIB} -ldeflate -llzma -lbz2 -lz -lm -lpthread

${BIN_TARGET}:${OBJ} | make_bin_dir
	$(CXX) $(LDFLAGS) $(OBJ) $(LDLIBS) -o $@
    
${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp | make_obj_dir
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
.PHONY: clean make_obj_dir make_bin_dir install
clean:
	rm -f $(OBJ) $(DEP) $(BIN_TARGET) $(TARGET)

make_obj_dir:
	mkdir -p $(DIR_OBJ)

make_bin_dir:
	mkdir -p $(DIR_BIN)

install: $(BIN_TARGET)
	install $(BIN_TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."

-include $(DEP)
