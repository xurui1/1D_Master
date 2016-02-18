//Here I calculate the homogeneous free energy using SCFT

double homogfE(double *mu, double **chimatrix,double *f){
    //chain partition functions
    double QT_ave,QD_ave,QC_ave;
    
    //activities
    double act_h=exp(kappaC*(mu[5]));
    double act_t=exp(kappa_ABA*(mu[2]));
    double act_d=exp(mu[0]);
    
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
        
        QT_ave = exp(-kappa_ABA*(w_ave[2]*f[2]+w_ave[3]*f[3]+w_ave[4]*f[4]));
        QD_ave = exp(-(w_ave[0]*f[0]+w_ave[1]*f[1]));
        QC_ave = exp(-kappaC*w_ave[5]*f[5]);
        
        //diblock concentration
        p_ave[0]=act_d*QD_ave*f[0];
        p_ave[1]=act_d*QD_ave*f[1];
        
        //triblock concentration
        p_ave[2]=act_t*QT_ave*f[2];
        p_ave[3]=act_t*QT_ave*f[3];
        p_ave[4]=act_t*QT_ave*f[4];
        
        //Homopolymer concentration
        p_ave[5]=act_h*QC_ave;
        
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
    fE_hom=f_int/2.0-f_omeg;
    fE_hom-=act_t*QT_ave/kappa_ABA;//triblock
    fE_hom-=act_d*QD_ave;          //diblock
    fE_hom-=act_h*QC_ave/kappaC;
    
    return fE_hom;
}


//Here I calculate the homogeneous free energy using SCFT

double homogfE_new(double *mu, double **chimatrix, double *f){
    
    //homogeneous concentrations
    double pA1_ave,pA2_ave,pA3_ave;
    double pB1_ave,pB2_ave;
    double pC_ave;
    
    //homogeneous chemical potentials
    double wA1_ave,wA2_ave,wA3_ave;
    double wB1_ave,wB2_ave;
    double wC_ave;
    
    //difference in hom chemical potential
    double dwA1_ave,dwA2_ave,dwA3_ave;
    double dwB1_ave,dwB2_ave;
    double dwC_ave,dpp_ave;
    double eta_ave;
    
    //homogeneous energies
    double f_int, f_omeg,fE_hom;
    double *p_vect;
    double *w_vect;
    
    p_vect=create_1d_double_array(6,"p_vect");
    w_vect=create_1d_double_array(6,"w_vect");
    
    f_int=0.0;
    f_omeg=0.0;
    
    //set initial change in chemical potentials
    dwA1_ave=0.0;
    dwA2_ave=0.0;
    dwA2_ave=0.0;
    dwB1_ave=0.0;
    dwB2_ave=0.0;
    dwC_ave=0.0;
    
    eta_ave=0.0;
    
    //set intial concentrations
    pA1_ave=0.002;
    pB1_ave=pA1_ave;
    pA2_ave=0.002;
    pB2_ave=pA2_ave;
    pA3_ave=pA2_ave;
    pC_ave=1.0-(pA1_ave+pA2_ave+pA3_ave+pB1_ave+pB2_ave);
    
    //set initial chemical potentials
    wA1_ave=chimatrix[0][1]*pB1_ave+chimatrix[0][2]*pA2_ave+chimatrix[0][3]*pB2_ave+chimatrix[0][4]*pA3_ave+chimatrix[0][5]*pC_ave+eta_ave;
    
    wB1_ave=chimatrix[1][0]*pA1_ave+chimatrix[1][2]*pA2_ave+chimatrix[1][3]*pB2_ave+chimatrix[1][4]*pA3_ave+chimatrix[1][5]*pC_ave+eta_ave;
    
    wA2_ave=chimatrix[2][0]*pA1_ave+chimatrix[2][1]*pB1_ave+chimatrix[2][3]*pB2_ave+chimatrix[2][4]*pA3_ave+chimatrix[2][5]*pC_ave+eta_ave;
    
    wB2_ave=chimatrix[3][0]*pA1_ave+chimatrix[3][1]*pB1_ave+chimatrix[3][2]*pA2_ave+chimatrix[3][4]*pA3_ave+chimatrix[3][5]*pC_ave+eta_ave;
    
    wA3_ave=chimatrix[4][0]*pA1_ave+chimatrix[4][1]*pB1_ave+chimatrix[4][2]*pA2_ave+chimatrix[4][3]*pB2_ave+chimatrix[4][5]*pC_ave+eta_ave;
    
    wC_ave=chimatrix[5][0]*pA1_ave+chimatrix[5][1]*pB1_ave+chimatrix[5][2]*pA2_ave+chimatrix[5][3]*pB2_ave+chimatrix[5][4]*pA3_ave+eta_ave;
    
    
    for (int i=0;i<10000000;i++){
        
        eta_ave=eta_ave-0.05*(1.0-(pA1_ave+pA2_ave+pA3_ave+pB1_ave+pB2_ave+pC_ave));
        
        //diblock concentration
        pA1_ave=exp(mu[0]-wA1_ave*f[0]-wB1_ave*f[1])*f[0];
        pB1_ave=exp(mu[0]-wA1_ave*f[0]-wB1_ave*f[1])*f[1];
        
        //triblock concentration
        pA2_ave=exp(2.0*mu[1]-wA2_ave*f[0]/2.0-wB2_ave*f[1]-wA3_ave*f[0]/2.0)*f[0]/2.0;
        pB2_ave=exp(2.0*mu[1]-wA2_ave*f[0]/2.0-wB2_ave*f[1]-wA3_ave*f[0]/2.0)*(2.0*f[1])/2.0;
        pA3_ave=exp(2.0*mu[1]-wA2_ave*f[0]/2.0-wB2_ave*f[1]-wA3_ave*f[0]/2.0)*f[0]/2.0;
        
        //Homopolymer concentration
        pC_ave=exp(kappaC*(mu[2]-wC_ave));
        
        //change in chemical potential for diblock
        dwA1_ave=(chimatrix[0][1]*pB1_ave+chimatrix[0][2]*pA2_ave+chimatrix[0][3]*pB2_ave+chimatrix[0][4]*pA3_ave+chimatrix[0][5]*pC_ave+eta_ave)-wA1_ave;
        dwB1_ave=(chimatrix[1][0]*pA1_ave+chimatrix[1][2]*pA2_ave+chimatrix[1][3]*pB2_ave+chimatrix[1][4]*pA3_ave+chimatrix[1][5]*pC_ave+eta_ave)-wB1_ave;
        
        //change in chemical potential for triblock
        dwA2_ave=(chimatrix[2][0]*pA1_ave+chimatrix[2][1]*pB1_ave+chimatrix[2][3]*pB2_ave+chimatrix[2][4]*pA3_ave+chimatrix[2][5]*pC_ave+eta_ave)-wA2_ave;
        dwB2_ave=(chimatrix[3][0]*pA1_ave+chimatrix[3][1]*pB1_ave+chimatrix[3][2]*pA2_ave+chimatrix[3][4]*pA3_ave+chimatrix[3][5]*pC_ave+eta_ave)-wB2_ave;
        dwA3_ave=(chimatrix[4][0]*pA1_ave+chimatrix[4][1]*pB1_ave+chimatrix[4][2]*pA2_ave+chimatrix[4][3]*pB2_ave+chimatrix[4][5]*pC_ave+eta_ave)-wA3_ave;
        
        //change in chemical potential for homopolymer
        dwC_ave=(chimatrix[5][0]*pA1_ave+chimatrix[5][1]*pB1_ave+chimatrix[5][2]*pA2_ave+chimatrix[5][3]*pB2_ave+chimatrix[5][4]*pA3_ave+eta_ave)-wC_ave;
        
        //change in total concentration
        dpp_ave=1.0-(pA1_ave+pA2_ave+pA3_ave+pB1_ave+pB2_ave+pC_ave);
        
        //update chemical potentials
        wA1_ave=wA1_ave+0.005*dwA1_ave;
        wB1_ave=wB1_ave+0.005*dwB1_ave;
        wA2_ave=wA2_ave+0.005*dwA2_ave;
        wB2_ave=wB2_ave+0.005*dwB2_ave;
        wA3_ave=wA3_ave+0.005*dwA3_ave;
        wC_ave=wC_ave+0.005*dwC_ave;
        
    }
    //build concentration vector
    p_vect[0]=pA1_ave;
    p_vect[1]=pB1_ave;
    p_vect[2]=pA2_ave;
    p_vect[3]=pB2_ave;
    p_vect[4]=pA3_ave;
    p_vect[5]=pC_ave;
    
    phi_bulk=0.0;
    for (int i=0;i<5;i++){
        phi_bulk +=p_vect[i];
    }
    
    //build chemical potential vector
    w_vect[0]=wA1_ave;
    w_vect[1]=wB1_ave;
    w_vect[2]=wA2_ave;
    w_vect[3]=wB2_ave;
    w_vect[4]=wA3_ave;
    w_vect[5]=wC_ave;
    
    
    //calculate interaction and potential energies
    for (int i=0;i<6;i++){
        for (int j=0;j<6;j++){
            f_int+=p_vect[i]*p_vect[j]*chimatrix[i][j];
        }
        f_omeg+=p_vect[i]*w_vect[i];
    }
    
    //combine all energies
    fE_hom=f_int/2.0-f_omeg-(exp(mu[0]-wA1_ave*f[0]-wB1_ave*f[1]));         //diblock
    fE_hom-=(exp(2.0*mu[1]-wA2_ave*f[0]-wA3_ave*f[0]-2.0*wB2_ave*f[1])/2.0);//triblock
    fE_hom-=(exp(kappaC*(mu[2]-wC_ave))/kappaC);                              //homopolymer
    
    return fE_hom;
}



