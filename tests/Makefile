# Simple Makefile for AudioTransportN tests

CXX = c++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -g
INCLUDES = -I../deps/flucoma-core/include -I../deps/eigen -I../deps/hisstools_library -I../deps

# FluCoMa requires these defines
DEFINES = -DEIGEN_MPL2_ONLY -DEIGEN_HAS_CXX11=1

TARGET = test_audiotransportn
SOURCES = TestAudioTransportN.cpp ../src/Explorer/AudioTransportN.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all run clean