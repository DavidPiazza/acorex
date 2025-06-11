# Simple Makefile for AudioTransportN tests

CXX = c++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -g
# Project dependencies
INCLUDES = -I../deps -I../deps/memory -I../deps/flucoma-core -I../deps/hisstools_library -I../deps/Eigen

# Core openFrameworks modules
INCLUDES += -I../../../../libs/openFrameworks
INCLUDES += -I../../../../libs/openFrameworks/3d
INCLUDES += -I../../../../libs/openFrameworks/app
INCLUDES += -I../../../../libs/openFrameworks/communication
INCLUDES += -I../../../../libs/openFrameworks/events
INCLUDES += -I../../../../libs/openFrameworks/gl
INCLUDES += -I../../../../libs/openFrameworks/graphics
INCLUDES += -I../../../../libs/openFrameworks/math
INCLUDES += -I../../../../libs/openFrameworks/sound
INCLUDES += -I../../../../libs/openFrameworks/types
INCLUDES += -I../../../../libs/openFrameworks/utils
INCLUDES += -I../../../../libs/openFrameworks/video

# Third-party libraries
INCLUDES += -I../../../../libs/boost/include
INCLUDES += -I../../../../libs/cairo/include/cairo
INCLUDES += -I../../../../libs/curl/include
INCLUDES += -I../../../../libs/fmod/include
INCLUDES += -I../../../../libs/FreeImage/include
INCLUDES += -I../../../../libs/freetype/include
INCLUDES += -I../../../../libs/glew/include
INCLUDES += -I../../../../libs/glfw/include
INCLUDES += -I../../../../libs/glm/include
INCLUDES += -I../../../../libs/json/include
INCLUDES += -I../../../../libs/pugixml/include
INCLUDES += -I../../../../libs/rtAudio/include
INCLUDES += -I../../../../libs/tess2/include
INCLUDES += -I../../../../libs/uriparser/include
INCLUDES += -I../../../../libs/utf8/include

# FluCoMa requires these defines
DEFINES = -DEIGEN_MPL2_ONLY -DEIGEN_HAS_CXX11=1

# Platform specific libraries
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    # Base libraries
    LDFLAGS = -L../libs -lfoonathan_memory-0.7.3
    # macOS frameworks
    LDFLAGS += -framework Accelerate
    LDFLAGS += -framework CoreFoundation
    LDFLAGS += -framework CoreGraphics
    LDFLAGS += -framework CoreServices
    LDFLAGS += -framework OpenGL
    LDFLAGS += -framework IOKit
    LDFLAGS += -framework Cocoa
    LDFLAGS += -framework AVFoundation
    LDFLAGS += -framework CoreMedia
    LDFLAGS += -framework CoreVideo
    LDFLAGS += -framework QTKit
    LDFLAGS += -framework CoreAudio
    LDFLAGS += -framework AudioToolbox
    LDFLAGS += -framework Security
    # OpenFrameworks path
    OF_LIBS_PATH = ../../../../libs/openFrameworksCompiled/lib/osx
    # Additional library paths
    LDFLAGS += -L$(OF_LIBS_PATH)
    LDFLAGS += -L../../../../libs/boost/lib/osx
    LDFLAGS += -L../../../../libs/cairo/lib/osx
    LDFLAGS += -L../../../../libs/curl/lib/osx
    LDFLAGS += -L../../../../libs/FreeImage/lib/osx
    LDFLAGS += -L../../../../libs/freetype/lib/osx
    LDFLAGS += -L../../../../libs/glew/lib/osx
    LDFLAGS += -L../../../../libs/glfw/lib/osx
    LDFLAGS += -L../../../../libs/rtAudio/lib/osx
    LDFLAGS += -L../../../../libs/tess2/lib/osx
    LDFLAGS += -L../../../../libs/uriparser/lib/osx
else
    LDFLAGS = -lfftw3 -L../libs -lfoonathan_memory-0.7.3
    OF_LIBS_PATH = ../../../../libs/openFrameworksCompiled/lib/linux64
endif

# Test targets
TARGET1 = test_audiotransportn
TARGET2 = test_phase_interpolation
TARGET3 = test_transport_frame_cache
TARGET4 = test_audio_morph_engine

SOURCES1 = TestAudioTransportN.cpp ../src/Explorer/AudioTransportN.cpp ../src/Explorer/SIMDUtils.cpp
SOURCES2 = TestPhaseInterpolation.cpp ../src/Explorer/AudioTransportN.cpp ../src/Explorer/SIMDUtils.cpp
SOURCES3 = TestTransportFrameCache.cpp
SOURCES4 = TestAudioMorphEngine.cpp ../src/Explorer/AudioMorphEngine.cpp ../src/Explorer/AudioTransportN.cpp ../src/Explorer/SIMDUtils.cpp

OBJECTS1 = $(SOURCES1:.cpp=.o)
OBJECTS2 = $(SOURCES2:.cpp=.o)
OBJECTS3 = $(SOURCES3:.cpp=.o)
OBJECTS4 = $(SOURCES4:.cpp=.o)

all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

$(TARGET1): $(OBJECTS1)
	$(CXX) $(CXXFLAGS) $(OBJECTS1) $(LDFLAGS) -o $(TARGET1)

$(TARGET2): $(OBJECTS2)
	$(CXX) $(CXXFLAGS) $(OBJECTS2) $(LDFLAGS) -o $(TARGET2)

$(TARGET3): $(OBJECTS3)
	$(CXX) $(CXXFLAGS) $(OBJECTS3) $(LDFLAGS) -pthread -o $(TARGET3)

$(TARGET4): $(OBJECTS4)
	$(CXX) $(CXXFLAGS) $(OBJECTS4) $(LDFLAGS) -pthread \
		-lopenFrameworks \
		../../../../libs/fmod/lib/osx/libfmod.dylib \
		../../../../libs/freetype/lib/osx/freetype.a \
		../../../../libs/FreeImage/lib/osx/freeimage.a \
		../../../../libs/cairo/lib/osx/cairo.a \
		../../../../libs/cairo/lib/osx/pixman-1.a \
		../../../../libs/cairo/lib/osx/png.a \
		../../../../libs/boost/lib/osx/boost_system.a \
		../../../../libs/boost/lib/osx/boost_filesystem.a \
		../../../../libs/curl/lib/osx/curl.a \
		../../../../libs/rtAudio/lib/osx/rtaudio.a \
		../../../../libs/tess2/lib/osx/tess2.a \
		../../../../libs/glfw/lib/osx/glfw3.a \
		../../../../libs/glew/lib/osx/glew.a \
		../../../../libs/uriparser/lib/osx/uriparser.a \
		-lz -lbz2 -liconv \
		-o $(TARGET4)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEFINES) -c $< -o $@

# Convenience targets
run: $(TARGET1)
	./$(TARGET1)

run_phase: $(TARGET2)
	./$(TARGET2)

run_cache: $(TARGET3)
	./$(TARGET3)

run_morph: $(TARGET4)
	./$(TARGET4)

run_all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)
	@echo "=== Running AudioTransportN Tests ==="
	./$(TARGET1)
	@echo "\n=== Running Phase Interpolation Tests ==="
	./$(TARGET2)
	@echo "\n=== Running TransportFrameCache Tests ==="
	./$(TARGET3)
	@echo "\n=== Running Audio Morph Engine Tests ==="
	./$(TARGET4)

clean:
	rm -f $(OBJECTS1) $(OBJECTS2) $(OBJECTS3) $(OBJECTS4) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

.PHONY: all run run_phase run_cache run_morph run_all clean