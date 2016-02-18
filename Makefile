LFLAGS = -lm  -O3 -stdc++0x

main: main.cpp 
	c++ $(LFLAGS) -o $@ $(MOBLIB) main.cpp $(LIBS)
