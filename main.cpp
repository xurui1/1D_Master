
#include "global.h"

typedef vector<double> Row;
typedef vector<Row> Matrix;

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
#include "./mods/mod_width.h"
#include "./mods/mod_radius.h"
#include "./mods/mod_phi.h"
#include "./mods/mod_phif50.h"
#include "./mods/mod_main.h"


int main( ){
    
    int nradii=15,nfa=21;                               //number of radius & fa measurements
    vector <double>A;
    vector <double>B;
    vector <double>C;                                       //Fitting parameters

    Matrix w(ChainType,Row(Nr));
    Matrix chiMatrix(ChainType,Row(ChainType));
    Matrix phi(ChainType,Row(Nr));
    
    vector<double>eta(Nr);
    vector<double>chi(ChainType);
    vector<double>f(ChainType);
    vector<double>mu(ChainType);
    vector<double>mu_vector(nfa);
    vector<double>dFE(nradii);
    vector<int>Ns(ChainType);
    
    //Set parameters & interaction matrix
    parameters(chi,f,Ns,mu);
    Xmatrix(chiMatrix,chi);

    //mod_width(f,mu,chiMatrix,w,phi,eta,Ns,chi,nfa);
    
    //calculate radius of membrane center
    mod_radius(f,mu,chiMatrix,w,phi,eta,Ns,chi,A,B,C,nfa,mu_vector);
    
    //reset parameters
    //parameters(chi,f,Ns,mu);
    
    //calculate concentration profiles
    //mod_phi(f,mu,chiMatrix,w,phi,eta,Ns,chi,nfa,A,B,C,nradii,dFE,mu_vector);
    mod_phif50(f,mu,chiMatrix,w,phi,eta,Ns,chi);
    

    
    ofstream outputrad;
    outputrad.open("./results/radius_fit.dat");
    outputrad<<A[0]<<" "<<B[0]<<" "<<C[0]<<endl;
    
    outputrad.close();
    
    //reset parameters
    parameters(chi,f,Ns,mu);
    
    //main function for finding bending moduli
    mod_main(f,mu,chiMatrix,w,phi,eta,Ns,chi,nfa,A,B,C,nradii,dFE,mu_vector);
    
    
    
    return 0;
}
