SWIG = swig
SWIG_FLAGS = -python -c++

PREPROCESSOR_SYMBOLS = -DSWIG

CPP = g++
CPP_FLAGS = -pedantic -Wall -Wextra -std=c++11 -O3 -g $(PREPROCESSOR_SYMBOLS) #-Werror is too much hassle with the generated interface

LOCAL_INC_DIR = ../include
LOCAL_LIB_DIR = ../lib
PYTHON_INC = /usr/include/python2.7
INCLUDE = -I$(LOCAL_INC_DIR)\
	  -I$(PYTHON_INC)

OBJECTS = $(LOCAL_LIB_DIR)/Dijkstra.o\
	  $(LOCAL_LIB_DIR)/Graph.o\
	  RCSP_wrap.o

.PHONY: python_binding

python_binding: _RCSP.so

_RCSP.so: $(OBJECTS)
	$(CPP) -shared $^ -o $@

RCSP_wrap.o: RCSP_wrap.cxx
	$(CPP) -c -fPIC $(CPP_FLAGS) $< $(INCLUDE) -o $@

RCSP_wrap.cxx : RCSP.i
	$(SWIG) $(SWIG_FLAGS) $<

clean:
	rm -f *_wrap.*
	rm -f _RCSP.so RCSP.py RCSP.pyc
