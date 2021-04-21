CXX = g++ -std=c++11
CXXFLAGS = -O3 -Wall -Wno-unused-result -mavx2 -pthread #-g #-lmetis

SOURCES = $(shell find . | grep -e ".hpp")
OBJECTS = $(SOURCES:%.hpp=%.o)
OBJECTS := $(OBJECTS)

all: mc reorder loi org bp

mc: $(OBJECTS) mc.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@

loi: $(OBJECTS) loi_mc.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@

bp: $(OBJECTS) bp_mc.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@

org: $(OBJECTS) org_mc.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@

reorder: $(OBJECTS) reorder.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@

%.o: %.cpp %.hpp util.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf mc bp org loi reorder $(OBJECTS)