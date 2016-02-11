/*********************A simple test function to run SCFT*************************************************/
void mod_0 (vector <double> &f,vector <double> &mu,Matrix &chiMatrix,Matrix &w,Matrix &phi,vector <double> &eta,vector <int> &Ns,vector <double> &chi){

    //fA = fB = 0.5
    for (int i=0;i<ChainType-1;i++){
        f[i]=0.5;
        Ns[i]=f[i]*Ds;
    }
    Ns[3] = Ds; //B block in triblock is longer
    
    time_point<high_resolution_clock> beginning,end;
    beginning=high_resolution_clock::now();
    
    //Initiate omega field
    omega(w);
    
    //calculate homogeneous free energy
    double fE_hom=homogfE(mu,chiMatrix,f);
    
    //calculate volume
    double volume=vol();
    
    //calculate free energy minus homogeneneous free energy
    double dFE=FreeEnergy(w,phi,eta,Ns,chi,chiMatrix,mu,volume,f,2*Nr/5,0);
    
    end = high_resolution_clock::now();
    duration<double> measured_time = end-beginning;
    
    
    cout<<"Run Time: "<<measured_time.count()<<endl;
}