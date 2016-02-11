LFLAGS = -lm  -O3 -std=c++11

main: main.cpp 
	g++ $(LFLAGS) -o $@ $(MOBLIB) main.cpp $(LIBS)
