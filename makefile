CXX     := g++

CFLAGS  := -g -O3 -std=c++11
LDFLAGS := -g -lm

LIBS    := -L/usr/local/lib -lm -lstdc++
INCLUDE := -I./include -I/usr/include/eigen3 -I/usr/local/include
OBJ_DIR := obj
EX_DIR  := data-output
SRC     := $(wildcard src/*.cpp)

ifdef FLAGS
	comma:= ,
	empty:=
	space:= $(empty) $(empty)
	DFLAGS = $(subst $(comma), $(space), $(FLAGS))

	ifeq ($(findstring OPENMP, $(DFLAGS)), OPENMP)
		CFLAGS = -g -O3 -std=c++11 -fopenmp
		LDFLAGS = -g -lm -fopenmp
	endif
endif

OBJ     := $(addprefix obj/,$(notdir $(SRC:.cpp=.o)))
TARGET  := bin/fem
DEPENDS := $(OBJ:.o=.d)

all: $(TARGET) $(EX_DIR)

$(TARGET): $(OBJ)
	$(CXX) $(INCLUDE) -o $@ $(OBJ) $(CFLAGS) $(LIBS)

obj/%.o : src*/%.cpp
	@if [ ! -e OBJ_DIR ]; \
		then mkdir -p obj; \
		fi
	$(CXX) $(CFLAGS) $(INCLUDE) -o $@ -c $<  

# clean files
clean: 
	-rm -f $(OBJ) $(TARGET) $(DEPENDS) 
	-rm -r $(EX_DIR) $(OBJ_DIR)

-include $(DEPENDS)

$(EX_DIR):
	mkdir -p data-output

# clean files
.PHONY: all clean
