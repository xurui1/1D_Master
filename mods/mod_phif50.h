/**************My function for calculating the concentration profile for an FA = 0.5 membrane***********************************************/

void mod_phif50(vector <double> &f,vector <double> &mu,Matrix &chiMatrix,Matrix &w,Matrix &phi,vector <double> &eta,vector <int> &Ns,vector <double> &chi){
    
    for (int i=0;i<ChainType-1;i++){
        f[i]=0.5;
        Ns[i]=f[i]*Ds;
    }
    Ns[3] = Ds; //B block in triblock is longer
    
    double fE_hom;
    
    //Initiate omega field
    omega(w);
    
    //calculate volume
    double volume = vol();
    
    //Find tensionless mmb
    secant(w,phi,eta,Ns,chi,chiMatrix,mu,f,2*Nr/5);

    //calculate homogeneous free energy
    fE_hom=homogfE(mu,chiMatrix,f);
    
    //define pinning location
    int pin = 2*Nr/5;
    
    //calculate volume
    volume=vol();
    
    //reset radius
    r_0=0.5;
    
    //reset omega
    omega(w);

    
    //calculate free energy minus homogeneneous free energy
    double dFE=FreeEnergy(w,phi,eta,Ns,chi,chiMatrix,mu,volume,f,pin,1);
    
    
        
    //output concentration profile
    ofstream outfile;
    outfile.open("./results/phi_FA50.dat");
    for (int i=0;i<Nr;i++){
        outfile<< (double)i*dr<<" "<<phi[0][i]<<" "<<phi[1][i]<<" "<<phi[2][i]<<" "<<phi[3][i]<<" "<<phi[4][i]<<" "<<phi[5][i]<<endl;
    }
    
        
        
    
    
}