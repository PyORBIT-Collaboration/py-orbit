include ../conf/make_root_config

DIRS = $(patsubst %/, %, $(filter %/,$(shell ls -F)))

#objects files list for ORBIT executable
OBJS := $(foreach dir,$(DIRS),$(patsubst $(dir)/%.cc,$(dir)/obj/%.o,$(wildcard $(dir)/*.cc)))

#========rules=========================

all: compile

compile:
	@for dir in $(DIRS); do \
		(cd $$dir; if [ -f ./Makefile ]; then $(MAKE) compile; fi;); \
	done
	$(CXX) -o $(orbit_main) $(OBJS) $(LIBS) $(LINKFLAGS)

clean:
	@for dir in $(DIRS); do \
		(cd $$dir; if [ -f ./Makefile ]; then $(MAKE) clean; fi;); \
	done
	rm -rf $(orbit_main)
