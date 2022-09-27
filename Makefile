vpath %.cpp src/

PROGRAM=mini-pic

CXXFLAGS = -std=c++20 -O2 -g -fopenmp

all: $(PROGRAM)

$(PROGRAM):

%: %.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf *.o  $(PROGRAM)
