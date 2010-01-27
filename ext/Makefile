DIRS  = 
DIRS += ./ecloud
DIRS += ./laserstripping
DIRS += ./H0starkeffect

#========rules=========================

all: compile

compile:  
	@for dir in $(DIRS); do \
		(cd $$dir; if [ -f ./Makefile ]; then $(MAKE) compile; fi;); \
	done	

clean:
	@for dir in $(DIRS); do \
		(cd $$dir; if [ -f ./Makefile ]; then $(MAKE) clean; fi;); \
	done

