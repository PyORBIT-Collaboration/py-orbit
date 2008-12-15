include ../../conf/make_root_config

SRCS_WRAP = $(wildcard wrap_*.cc)
SRCS = $(wildcard *.cc)

OBJS = $(patsubst %.cc,./obj/%.o,$(SRCS))
OBJS_WRAP = $(patsubst wrap_%.cc,./obj/wrap_%.o,$(SRCS_WRAP))

INC  = $(wildcard *.hh)
INC += $(wildcard *.h)

#include files could be everywhere
INC_DIRS =$(filter %/,$(shell ls -F ../../src))
INCLUDES_LOCAL = $(patsubst %/, -I../../src/%, $(INC_DIRS))

#wrappers CC FLAGS
WRAPPER_FLAGS = -fno-strict-aliasing

#CXXFLAGS
CXXFLAGS += -fPIC

#shared library flags
SHARED_LIB = -shared

#tracker shared library
spacecharge_lib = spacecharge.so

#========rules=========================
compile: $(OBJS_WRAP) $(OBJS) $(INC)
	$(CXX) -fPIC $(SHARED_LIB) -o ../../modules/$(spacecharge_lib) $(OBJS) 

./obj/wrap_%.o : wrap_%.cc $(INC)
	$(CXX) $(CXXFLAGS) $(WRAPPER_FLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

./obj/%.o : %.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

clean:
	rm -rf ./obj/*.o
	rm -rf ../../modules/$(spacecharge_lib)
