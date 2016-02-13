
#include "global.h"

typedef vector<double> Row;
typedef vector<Row> Matrix;

#include "./functions/inputarguments.h"
#include "./functions/parameters.h"
#include "./functions/filename.h"
#include "./functions/vol.h"
#include "./functions/integrate.h"
#include "./functions/omega.h"
#include "./functions/TDMA.h"
#include "./functions/solvediffeq.h"
#include "./functions/phi.h"
#include "./functions/Q_partition.h"
#include "./functions/polymers.h"
#include "./functions/loop.h"
#include "./functions/conc.h"
#include "./functions/Incomp.h"
#include "./functions/output.h"
#include "./functions/fE.h"
#include "./functions/homogfE.h"
#include "./functions/secant.h"
#include "./functions/radius.h"
#include "./functions/calcexcess.h"
#include "./functions/FreeEnergy.h"
#include "./functions/curvefitting.h"
#include "./mods/mod_0.h"
#include "./mods/mod_width.h"
#include "./mods/mod_radius.h"
#include "./mods/mod_phi.h"
#include "./mods/mod_phif50.h"
#include "./mods/mod_loop.h"
#include "./mods/mod_main.h"


int main( int argc, char* argv[] ){
    
    int nradii=15,nfa=21;                               //number of radius & fa measurements
    double *A=create_1d_double_array(1,"A");
    double *B=create_1d_double_array(1,"B");
    double *C=create_1d_double_array(1,"C");
       //Fitting parameters

    double **w=create_2d_double_array(ChainType,Nr,"w");
    double **chiMatrix=create_2d_double_array(ChainType,ChainType,"chiMatrix");
    double **phi=create_2d_double_array(ChainType,Nr,"phi");
    
    
    double *eta=create_1d_double_array(Nr,"eta");
    double *chi=create_1d_double_array(ChainType,"chi");
    double *f=create_1d_double_array(ChainType,"f");
    double *mu=create_1d_double_array(ChainType,"mu");
    double *mu_vector=create_1d_double_array(nfa,"mu_vector");
    double *dFE=create_1d_double_array(nradii,"dFE");
    int *Ns=create_1d_integer_array(ChainType,"Ns");

    input_Arguments(argc, argv);
    
    //Set parameters & interaction matrix
    parameters(chi,f,Ns,mu);
    Xmatrix(chiMatrix,chi);
    
    if(atoi(argv[3])==0){
        mod_0(f,mu,chiMatrix,w,phi,eta,Ns,chi);
    }
    else if(atoi(argv[3])==1){
        mod_radius(f,mu,chiMatrix,w,phi,eta,Ns,chi,A,B,C,nfa,mu_vector); //find rad to centre mmb
        parameters(chi,f,Ns,mu); //reset parameters
        mod_phi(f,mu,chiMatrix,w,phi,eta,Ns,chi,nfa,A,B,C,nradii,dFE,mu_vector);
    }
    else if(atoi(argv[3])==2){
        mod_radius(f,mu,chiMatrix,w,phi,eta,Ns,chi,A,B,C,nfa,mu_vector); //find rad to centre mmb
        parameters(chi,f,Ns,mu); //reset parameters
        mod_main(f,mu,chiMatrix,w,phi,eta,Ns,chi,nfa,A,B,C,nradii,dFE,mu_vector);
    }
    else if(atoi(argv[3])==3){
        mod_phif50(f,mu,chiMatrix,w,phi,eta,Ns,chi);
    }
    else if (atoi(argv[3])==4){
        mod_width(f,mu,chiMatrix,w,phi,eta,Ns,chi,nfa);
    }
    else if (atoi(argv[3])==5){
        mod_loop(f,mu,chiMatrix,w,phi,eta,Ns,chi,nfa);
    }
    else{
        cout<<"The mod does not exist, try again"<<endl;
    }

    
    
    return 0;
}
