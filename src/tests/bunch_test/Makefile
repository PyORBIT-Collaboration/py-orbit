include ../../../conf/make_root_config

TEST_MAIN = test_bunch

SRCS  = $(wildcard *.cc)
OBJS  = $(patsubst %.cc,./obj/%.o,$(SRCS))

#==========PUT all object files HERE==============
OBJS += ../../mpi/obj/orbit_mpi.o
OBJS += ../../orbit/obj/ParticleAttributes.o
OBJS += ../../orbit/obj/OrbitConst.o
OBJS += ../../orbit/obj/Bunch.o
OBJS += ../../orbit/obj/SyncPart.o
OBJS += ../../orbit/obj/ParticleMacroSize.o
OBJS += ../../orbit/obj/ParticleAttributesFactory.o

OBJS += ../../utils/obj/StringUtils.o
OBJS += ../../utils/obj/BufferStore.o
OBJS += ../../utils/obj/AttributesBucket.o

INC  = $(wildcard *.hh)
INC += $(wildcard *.h)

#include files could be everywhere
INC_DIRS =$(filter %/,$(shell ls -F ../..))
INCLUDES_LOCAL = $(patsubst %/, -I../../%, $(INC_DIRS))

#========rules=========================
compile: $(OBJS) $(INC)
	$(CXX) -o $(TEST_MAIN) $(OBJS) $(LIBS) $(LINKFLAGS)

./obj/%.o : %.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

clean:
	rm -rf $(TEST_MAIN)
	rm -rf ./obj/*.o
