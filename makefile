WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

SRC_DIR = src
INC_DIR = include
OBJ_DIR = obj
BIN_DIR = bin

INC = -I$(INC_DIR)
CFLAGS = -Wall -fexceptions
RESINC = 
LIBDIR = 
LIB = 
LDFLAGS = 

INC_DEBUG = $(INC)
CFLAGS_DEBUG = $(CFLAGS) -g
RESINC_DEBUG = $(RESINC)
LIBDIR_DEBUG = $(LIBDIR)
LIB_DEBUG = $(LIB)
LDFLAGS_DEBUG = $(LDFLAGS)
OBJDIR_DEBUG = $(OBJ_DIR)/Debug
DEP_DEBUG = 
OUT_DEBUG = $(BIN_DIR)/Debug/GBLA

INC_RELEASE = $(INC)
CFLAGS_RELEASE = $(CFLAGS) -O2
RESINC_RELEASE = $(RESINC)
LIBDIR_RELEASE = $(LIBDIR)
LIB_RELEASE = $(LIB)
LDFLAGS_RELEASE = $(LDFLAGS) -s
OBJDIR_RELEASE = $(OBJ_DIR)/Release
DEP_RELEASE = 
OUT_RELEASE = $(BIN_DIR)/Release/GBLA

SRC_FILES = $(SRC_DIR)/GBLA_LC.cpp $(SRC_DIR)/Multinomial.cpp $(SRC_DIR)/Term.cpp $(SRC_DIR)/Code.cpp $(SRC_DIR)/GF.cpp $(SRC_DIR)/GFElem.cpp $(SRC_DIR)/Matrix.cpp $(SRC_DIR)/Word.cpp $(SRC_DIR)/main.cpp
OBJ_DEBUG = $(patsubst $(SRC_DIR)/%.cpp, $(OBJDIR_DEBUG)/%.o, $(SRC_FILES))
OBJ_RELEASE = $(patsubst $(SRC_DIR)/%.cpp, $(OBJDIR_RELEASE)/%.o, $(SRC_FILES))

all: debug release

clean: clean_debug clean_release

before_debug: 
	test -d $(BIN_DIR)/Debug || mkdir -p $(BIN_DIR)/Debug
	test -d $(OBJDIR_DEBUG) || mkdir -p $(OBJDIR_DEBUG)

after_debug: 
	cp codes.txt $(BIN_DIR)/Debug/

debug: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG) $(DEP_DEBUG)
	$(LD) $(LIBDIR_DEBUG) -o $(OUT_DEBUG) $(OBJ_DEBUG) $(LDFLAGS_DEBUG) $(LIB_DEBUG)

$(OBJDIR_DEBUG)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c $< -o $@

clean_debug: 
	rm -f $(OBJ_DEBUG) $(OUT_DEBUG)
	rm -rf $(BIN_DIR)/Debug
	rm -rf $(OBJDIR_DEBUG)

before_release: 
	test -d $(BIN_DIR)/Release || mkdir -p $(BIN_DIR)/Release
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)

after_release: 
	cp codes.txt $(BIN_DIR)/Release/

release: before_release out_release after_release

out_release: before_release $(OBJ_RELEASE) $(DEP_RELEASE)
	$(LD) $(LIBDIR_RELEASE) -o $(OUT_RELEASE) $(OBJ_RELEASE) $(LDFLAGS_RELEASE) $(LIB_RELEASE)

$(OBJDIR_RELEASE)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c $< -o $@

clean_release: 
	rm -f $(OBJ_RELEASE) $(OUT_RELEASE)
	rm -rf $(BIN_DIR)/Release
	rm -rf $(OBJDIR_RELEASE)

.PHONY: before_debug after_debug clean_debug before_release after_release clean_release
