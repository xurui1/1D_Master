/*
 In this file we analyze the command line input arguments.
 - The arguments read in are in order:
 - 1-> P=plane, C = Cylinder, S=Sphere
 - 2->
 - 3-> Mod number, which determines what type of calculatin is being done
 */
void input_Arguments(int numb_of_args, char* arg_input[]){
    
   
    
    if(strcmp( arg_input[1], "P") == 0){
        Coord=1;
    }else if(strcmp( arg_input[1], "C") == 0){
        Coord=2;
    }else if(strcmp( arg_input[1], "S") == 0){
        Coord=3;
    }else{
        std::cout<<"The phase you have chosen does not exists! check: inputArguments.h"<<std::endl;
        exit(1);
    }
    
    if(strcmp( arg_input[2], "T") == 0){
        poly=0;
    }
    else if(strcmp( arg_input[2], "D") == 0){
        poly=1;
    }
    else if (strcmp( arg_input[2], "B") == 0){
        poly=2;
    }else{
        std::cout<<"The polymer you have chosen does not exists! check: inputArguments.h"<<std::endl;
        exit(1);
    }
    
   
};
