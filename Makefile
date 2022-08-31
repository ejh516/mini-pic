vpath %.cpp src/
vpath %.o		obj/

PROGRAM=mini-pic


all: $(PROGRAM)

$(PROGRAM):

%: %.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o obj/$@ -Jobj

clean:
	rm -rf obj/*  $(PROGRAM)
