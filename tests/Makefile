# Simple Makefile for AudioTransportN tests

CXX = c++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -g
INCLUDES = -I../deps -I../deps/memory -I../deps/flucoma-core -I../deps/hisstools_library -I../deps/Eigen

# FluCoMa requires these defines
DEFINES = -DEIGEN_MPL2_ONLY -DEIGEN_HAS_CXX11=1

# Platform specific libraries
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    LDFLAGS = -framework Accelerate -L../libs -lfoonathan_memory-0.7.3
else
    LDFLAGS = -lfftw3 -L../libs -lfoonathan_memory-0.7.3
endif

# Test targets
TARGET1 = test_audiotransportn
TARGET2 = test_phase_interpolation

SOURCES1 = TestAudioTransportN.cpp ../src/Explorer/AudioTransportN.cpp
SOURCES2 = TestPhaseInterpolation.cpp ../src/Explorer/AudioTransportN.cpp

OBJECTS1 = $(SOURCES1:.cpp=.o)
OBJECTS2 = $(SOURCES2:.cpp=.o)

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJECTS1)
	$(CXX) $(CXXFLAGS) $(OBJECTS1) $(LDFLAGS) -o $(TARGET1)

$(TARGET2): $(OBJECTS2)
	$(CXX) $(CXXFLAGS) $(OBJECTS2) $(LDFLAGS) -o $(TARGET2)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c $< -o $@

run: $(TARGET1)
	./$(TARGET1)

run_phase: $(TARGET2)
	./$(TARGET2)

run_all: $(TARGET1) $(TARGET2)
	@echo "=== Running AudioTransportN Tests ==="
	./$(TARGET1)
	@echo "\n=== Running Phase Interpolation Tests ==="
	./$(TARGET2)

clean:
	rm -f $(OBJECTS1) $(OBJECTS2) $(TARGET1) $(TARGET2)

.PHONY: all run run_phase run_all clean