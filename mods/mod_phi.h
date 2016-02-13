void mod_phi(double *f,double *mu,double **chiMatrix,double **w,double **phi,double *eta,int *Ns,double *chi, int nfa,double *A,double *B,double *C, int nradii,double *dFE,double *mu_vec){
    
    double *r_0vector= create_1d_double_array(nradii+1,"r_0vector");
    
    int counter=0;
    double volume,diameter;
    double OP;
    double fE_hom;
    
    int ds_increment = 0.4*(double)Ds/((double)nfa-1.0);
    
    for (int dds=0; f[0]<0.7 ;dds+= ds_increment){
        counter+=1;
        //Set parameters s
        updateparameters(f,Ns,dds);
        mu[5] = mu_vec[counter-1];                        //don't want to calc mu again
        fE_hom=homogfE(mu,chiMatrix,f);                 //calculate homog. fE
        omega(w);                                       //Initiate omega field
        
        double pin_location=10.8-A[0]-B[0]*f[0]-C[0]*f[0]*f[0];
        int pin = pin_location/dr;
        
        volume=vol();                                 //calculate volume
        OP = calcOP(phi,volume);                     //calculate order parameter
        //Set radius vector
        set_radius(r_0vector,nradii);
        
        r_0=0.5;                                        //reset radius
 
        volume=vol();
            
        //initialize omega field
        omega(w);
            
        //calculate free energy minus homogeneneous free energy
        dFE[0]=FreeEnergy(w,phi,eta,Ns,chi,chiMatrix,mu,f,pin,1);
        OP = calcOP(phi,volume);                    //calculate order parameter
        diameter = calc_excess(phi,volume); //calculate copolymer excess
            
        
        //output concentration profile
        outputphi_fa(phi,f[0],nfa);
            

        
    }
    
    
}