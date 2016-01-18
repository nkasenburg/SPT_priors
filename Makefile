MAYBE_USE_THESE_FLAGS = -Wshadow

# Controls which vectorisation extennsion to use.
VECT_EXT = #-DENABLE_AVX -mavx
#VECT_EXT = -DENABLE_SSE2 -msse2

# The C compiler is needed for compiling sleef.
CC = gcc
CC_FLAGS = -pedantic -Wall -Wextra -Werror -O3 $(VECT_EXT)


CPP = g++
CPP_FLAGS = -pedantic -Wall -Wextra -Werror -std=c++11 -pthread -O3 $(VECT_EXT) -g #-ftree-vectorizer-verbose=2  #-g # 

SRC_DIR = src
LOCAL_LIB_DIR = lib
LOCAL_INC_DIR = include

APP_DIR = apps

PYTHON_DIR = python

# Library for vectorisation
SLEEF_DIR = sleef-2.80
SLEEF_LIB_DIR = $(SLEEF_DIR)/simd
SLEEF_INC = $(SLEEF_DIR)


INCLUDE = -I$(LOCAL_INC_DIR)\
	  -I$(PARSER_INCLUDE)\
	  -I$(SLEEF_INC)

LIB_DIRS = -L$(LOCAL_LIB_DIR)\
	   -L$(SLEEF_LIB_DIR)

OBJECTS = $(LOCAL_LIB_DIR)/Graph.o\
	  $(LOCAL_LIB_DIR)/Dijkstra.o

RCSP = $(LOCAL_LIB_DIR)/RCSP.so

.PHONY: all cpp_binding python_binding apps clean objects

all: cpp_binding python_binding

cpp_binding: $(RCSP)

python_binding: $(OBJECTS)
	$(MAKE) -C $(PYTHON_DIR)

$(RCSP): $(OBJECTS)
	$(CPP) -shared $^ -o $@

$(LOCAL_LIB_DIR)/%.o: $(SRC_DIR)/%.cxx
	$(CPP) -c -fPIC $(CPP_FLAGS) $< $(INCLUDE) -o $@

$(SLEEF_LIB_DIR)/%.o: $(SLEEF_LIB_DIR)/%.c
	$(CC) -c -fPIC $(CC_FLAGS) $< $(INCLUDE) -o $@

$(APP_DIR)/%.o: $(APP_DIR)/%.cxx
	$(CPP) -c $(CPP_FLAGS) $< $(INCLUDE) -o $@

zip:
	zip $(PYTHON_DIR)/_RCSP.so $(PYTHON_DIR)/RCSP.py \
	shortest_path_bg.py brain_graph_sample.py utils.py

clean:
	rm -f $(OBJECTS) $(RCSP)
	$(MAKE) clean -C $(PYTHON_DIR)
