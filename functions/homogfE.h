//Here I calculate the homogeneous free energy using SCFT

double homogfE(vector <double> mu, Matrix chimatrix, vector <double> f){
    
    //homogeneous concentrations
    vector<double>p_ave(ChainType);
    
    //homogeneous chemical potentials
    vector<double>w_ave(ChainType);
    
    //difference in hom chemical potential
    double dpp_ave=0.0;
    double eta_ave=0.0;
    vector<double>dw_ave(ChainType);
    
    //homogeneous energies
    double f_int=0, f_omeg=0,fE_hom=0;
    
    //set intial concentrations
    p_ave[ChainType-1]=1.0;             //C concentration
    for (int i=0;i<ChainType-1;i++){
        p_ave[i]=0.002;
        p_ave[ChainType-1]-=p_ave[i];
    }
    
    //set initial chemical potentials
    for (int i=0;i<ChainType;i++){
        w_ave[i]=0.0;
        for (int j=0;j<ChainType;j++) {
            if (i != j){
                w_ave[i] += chimatrix[i][j]*p_ave[j];
            }
        }
        w_ave[i]+=eta_ave;
    }
    
    
    for (int i=0;i<10000000;i++){
        
        eta_ave=eta_ave-0.05*(1.0-(p_ave[0]+p_ave[1]+p_ave[2]+p_ave[3]+p_ave[4]+p_ave[5]));
        
        //diblock concentration
        p_ave[0]=exp(mu[0]-w_ave[0]*f[0]-w_ave[1]*f[1])*f[0];
        p_ave[1]=exp(mu[1]-w_ave[0]*f[0]-w_ave[1]*f[1])*f[1];
        
        //triblock concentration
        p_ave[2]=exp(kappa_ABA*mu[2]-w_ave[2]*f[2]/2.0-w_ave[3]*f[3]-w_ave[4]*f[4]/2.0)*f[2]/2.0;
        p_ave[3]=exp(kappa_ABA*mu[3]-w_ave[2]*f[2]/2.0-w_ave[3]*f[3]-w_ave[4]*f[4]/2.0)*(2.0*f[3])/2.0;
        p_ave[4]=exp(kappa_ABA*mu[4]-w_ave[2]*f[2]/2.0-w_ave[3]*f[3]-w_ave[4]*f[4]/2.0)*f[4]/2.0;
        
        //Homopolymer concentration
        p_ave[5]=exp(kappaC*(mu[5]-w_ave[5]));
        
        //change in chemical potential
        for (int i=0;i<ChainType;i++){
            dw_ave[i]=0.0;
            for (int j=0;j<ChainType;j++) {
                if (i != j){
                    dw_ave[i] += chimatrix[i][j]*p_ave[j];
                }
            }
            dw_ave[i]+=eta_ave;
            dw_ave[i]-=w_ave[i];
            
        }
       
        //change in total concentration
        dpp_ave=1.0-(p_ave[0]+p_ave[1]+p_ave[2]+p_ave[3]+p_ave[4]+p_ave[5]);
        
        //update chemical potentials
        for (int i=0;i<ChainType;i++){
            w_ave[i]+=0.005*dw_ave[i];
        }
       
    }
   
    phi_bulk=0.0;                       //concentration of non-homopolymer
    for (int i=0;i<ChainType-1;i++){
        phi_bulk +=p_ave[i];
    }
    
    //calculate interaction and potential energies
    for (int i=0;i<ChainType;i++){
        for (int j=0;j<ChainType;j++){
            f_int+=p_ave[i]*p_ave[j]*chimatrix[i][j];
        }
        f_omeg+=p_ave[i]*w_ave[i];
    }
    
    //change below to be consistent
    
    //combine all energies
    fE_hom=f_int/2.0-f_omeg-(exp(mu[0]-w_ave[0]*f[0]-w_ave[1]*f[1]));         //diblock
    fE_hom-=(exp(2.0*mu[2]-w_ave[2]*f[2]-w_ave[4]*f[4]-2.0*w_ave[3]*f[3])/kappa_ABA);//triblock
    fE_hom-=(exp(kappaC*(mu[5]-w_ave[5]))/kappaC);                              //homopolymer
    
    return fE_hom;
}


