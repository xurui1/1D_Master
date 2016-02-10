#include <stdio.h>     //Include the standard input/output libraries
#include <iostream>  //Cout and Cin etc.
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdlib.h>    //Include standard function libraries
#include <cmath>      //Use the math function libraries
#include <time.h>      //Call system time libraries to define the integer seed for random numbers

using namespace std;

#define Nr 200
#define Ds 200
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
double gamma_up = 0.05;
double epsilon_up = 0.05;
//set to 0.05 for chi<30
//set to 0.005 for 30<chi<40



