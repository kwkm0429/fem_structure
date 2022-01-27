CXX     := g++            # compiler
CFLAGS  := -g -O3 -MMD -MP -Wall -std=c++11 #-fopenmp  # -g : generate debug information, -MMD : output *.d , -MP : output dependent, -O2 : optimization
LDFLAGS := -g #-fopenmp             # math library -> add -lm
LIBS    :=                # declare .a, .so
INCLUDE := -I ./include   # -I : search .h in directory
SRC_DIR := src          # source directory
OBJ_DIR := obj          # .o , .d 
EX_DIR  := data-output     # .vtk
SOURCES := $(wildcard src/*.cpp) # all file, if you want except a file, use shell command
OBJS    := $(addprefix obj/,$(notdir $(SOURCES:.cpp=.o)))
TARGET  := run
DEPENDS := $(OBJS:.o=.d) # convert .o to .d

all: $(TARGET) data-output

$(TARGET): $(OBJS) $(LIBS)
	$(CXX) -o  $@ $(OBJS) $(LDFLAGS)

obj/%.o : src/%.cpp
	@if [ ! -e OBJ_DIR ]; \
		then mkdir -p obj; \
		fi
	$(CXX) $(CFLAGS) $(INCLUDE) -o $@ -c $<  

# clean files
clean: 
	-rm -f $(OBJS) $(TARGET) $(DEPENDS) 
	-rm -r $(EX_DIR) $(OBJ_DIR)

-include $(DEPENDS)

# create data-output 
data-output: # make data-output directory
	mkdir -p data-output

# clean files
.PHONY: all clean
