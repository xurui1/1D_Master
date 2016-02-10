/*********************A simple test function to run SCFT*************************************************/
void mod_0 (vector <double> &f,vector <double> &mu,Matrix &chiMatrix,Matrix &w,Matrix &phi,vector <double> &eta,vector <int> &Ns,vector <double> &chi){

    //fA = fB = 0.5
    for (int i=0;i<ChainType-1;i++){
        f[i]=0.5;
        Ns[i]=f[i]*Ds;
    }
    Ns[3] = Ds; //B block in triblock is longer
    
    //Initiate omega field
    omega(w);
    
    //calculate volume
    double volume = vol();
    
    //calculate homogeneous free energy
    double fE_hom=homogfE(mu,chiMatrix,f);
    
    //calculate volume
    volume=vol();
    
    //calculate free energy minus homogeneneous free energy
    double dFE=FreeEnergy(w,phi,eta,Ns,chi,chiMatrix,mu,volume,f,2*Nr/5,0);
    
    
}