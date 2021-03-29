CXX = g++ -std=c++11
CXXFLAGS = -O3 -Wall -Wno-unused-result -mavx2 -pthread -g #-lmetis

SOURCES = $(shell find . | grep -e ".hpp")
OBJECTS = $(SOURCES:%.hpp=%.o)
OBJECTS := $(OBJECTS)

all: mc reorder

mc: $(OBJECTS) mc.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@

reorder: $(OBJECTS) reorder.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@


%.o: %.cpp %.hpp util.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -rf mc reorder $(OBJECTS)