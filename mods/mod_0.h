/*********************A simple test function to run SCFT*************************************************/
void mod_0 (double *f,double *mu,double **chiMatrix,double **w,double **phi,double *eta,int *Ns,double *chi){

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
    double dFE=FreeEnergy(w,phi,eta,Ns,chi,chiMatrix,mu,f,2*Nr/5,0);
    
    end = high_resolution_clock::now();
    duration<double> measured_time = end-beginning;
    
    
    cout<<"Run Time: "<<measured_time.count()<<endl;
}