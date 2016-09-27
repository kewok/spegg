GTEST_DIR = ./../googletest/googletest
GMOCK_DIR = ./../googletest/googlemock

CFLAGS  = -O1
CFLAGS += -arch=sm_30

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include -isystem $(GMOCK_DIR)/include

# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wall -Wextra -pthread

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# All Google Mock headers. Note that all Google Test headers are
# included here too, as they are #included by Google Mock headers.
# Usually you shouldn't change this definition.	
GMOCK_HEADERS = $(GMOCK_DIR)/include/gmock/*.h \
                $(GMOCK_DIR)/include/gmock/internal/*.h \
                $(GTEST_HEADERS)

SPEGG_HEADERS = -I./../../header
OBJDIR_SPEGG = ./../../objdir

# Directory for individual tests
IND_TESTS = Individual_Tests

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = TestSPEGGspp # $(IND_TESTS)/WhichEqualToIntegerSuccess  $(IND_TESTS)/WhichEqualToFloatSuccess $(IND_TESTS)/WhichEqualToSubsetIntegerSuccess $(IND_TESTS)/AmplifyIntVectorSuccess  $(IND_TESTS)/WhichGreaterThanIntegerSuccess $(IND_TESTS)/WhichGreaterThanFloatSuccess $(IND_TESTS)/WhichGreaterThanSubsetIntegerSuccess

# Object directory for unit tests
OBJDIR = unitTests_objdir

#UTIL_CUDA_FILES := $(wildcard *.cu)
#UTIL_OBJ_FILES := $(addprefix $(OBJDIR)/,$(notdir $(UTIL_CUDA_FILES:.cu=.o)))

# Note all the objects to be created:

OBJS := $(OBJDIR)/Test_Core.o $(OBJDIR)/RemoveDeadTest.o 

# House-keeping build targets.
all : $(OBJDIR) $(OBJDIR_SPEGG) $(TESTS) $(IND_TESTS) TestSPEGGspp

clean : 
	rm -rf $(TESTS) gmock.a gmock_main.a $(OBJDIR) $(IND_TESTS)

# Make sure you have access to cuda lib
export LD_LIBRARY_PATH=/usr/local/lib /usr/local/cuda/lib

# Paths to sPEGG specific files
SPEGG_HEADERS = -I./../../header
SPEGG_SRC = ./../../src
CUDA = ${SPEGG_SRC}/cuda
UTIL = ${CUDA}/util

####################
#
#
# Priming Google test's own stuff
#
#
###################

# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on implementation details of
# Google Mock and Google Test, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Mock and
# Google Test compile fast and for ordinary users their source rarely
# changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) -I$(GMOCK_DIR) $(CXXFLAGS) \
            -c $(GTEST_DIR)/src/gtest-all.cc

gmock-all.o : $(GMOCK_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) -I$(GMOCK_DIR) $(CXXFLAGS) \
            -c $(GMOCK_DIR)/src/gmock-all.cc

gmock_main.o : $(GMOCK_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) -I$(GMOCK_DIR) $(CXXFLAGS) \
            -c $(GMOCK_DIR)/src/gmock_main.cc

gmock.a : gmock-all.o gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gmock_main.a : gmock-all.o gtest-all.o gmock_main.o
	$(AR) $(ARFLAGS) $@ $^


####################
#
#
# Building the local tests for sPEGG core code
#
#
###################

$(OBJDIR):
	mkdir $(OBJDIR) $(IND_TESTS)

# Compile the vanilla sPEGG functions:
#$(OBJDIR)/inds.o : ${SPECIES}/inds.cu
#	nvcc -c $(CFLAGS) ${HEADERS}  ${SPECIES}/inds.cu -o $(OBJDIR)/inds.o

# Compile the unit tests
$(OBJDIR)/%.o : %.cu
	nvcc -c $(CFLAGS) $(SPEGG_HEADERS) -I./../googletest/googletest/include/ -I./../googletest/googlemock/include/ -o $@ $<

$(OBJDIR)/Test_Core.o : Test_Core.cu
	nvcc -c $(CFLAGS) $(SPEGG_HEADERS) -I./../googletest/googletest/include/ -I./../googletest/googlemock/include/ Test_Core.cu -o $(OBJDIR)/Test_Core.o

####################
#
#
# Link the local tests to Google test
#
#
###################

# TO DEBUG and find the value of variable named var
# $(info $$var is [${var}])

TestSPEGGspp : $(OBJDIR_SPEGG)/*.o $(OBJS) gmock_main.a
	nvcc -O1 -lcurand -lrt -lcuda -lconfig++ -Xcompiler -I./../googletest/googlemock/include/ -I./../googletest/googletest/include/ --compiler-options -Wall --compiler-options -Wextra --compiler-options -pthread -lpthread $^ -o $@

