#include <stdio.h>     //Include the standard input/output libraries
#include <iostream>  //Cout and Cin etc.
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdlib.h>    //Include standard function libraries
#include <cmath>      //Use the math function libraries
#include <ctime>
#include <chrono>
#include "./functions/smemory.h"

using namespace std;
using namespace chrono;


#define Nr 200
#define Ds 100
#define ChainType 6
#define Pi 3.14159


double kappaC,kappa_ABA;
double r_0;
double phi_bulk;
double ds,dr;
int initial;
int Coord;
int poly;

//define my update parameters
double gamma_up = 0.01;
double epsilon_up = 0.01;
//set to 0.05 for chi<30
//set to 0.005 for 30<chi<40



