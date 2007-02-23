include ../../conf/make_root_config

SRCS = $(wildcard *.cc)

OBJS = $(patsubst %.cc,./obj/%.o,$(SRCS))

INC  = $(wildcard *.hh)
INC += $(wildcard *.h)

#include files could be everywhere
INC_DIRS =$(filter %/,$(shell ls -F ..))
INCLUDES_LOCAL = $(patsubst %/, -I../%, $(INC_DIRS))

#========rules=========================
compile: $(OBJS) $(INC)

./obj/%.o : %.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

clean:
	rm -rf ./obj/*.o
