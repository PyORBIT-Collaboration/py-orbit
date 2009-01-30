include ../conf/make_root_config

DIRS = $(filter-out test%,$(patsubst %/, %, $(filter %/,$(shell ls -F))))
#========rules=========================

all: compile

compile:
	@for dir in $(DIRS); do \
		(cd $$dir; if [ -f ./Makefile ]; then $(MAKE) compile; fi;); \
	done
	$(MAKE) --makefile=MakefileMain compile

clean:
	@for dir in $(DIRS); do \
		(cd $$dir; if [ -f ./Makefile ]; then $(MAKE) clean; fi;); \
	done
	rm -rf $(orbit_main)
