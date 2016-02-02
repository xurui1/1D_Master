void phi_total(Matrix phi, double volume){
    
    //Here I am calculating the total concentration of each species using a trapezoidal (?) rule.
    //This is an ugly function, and I'd like to rewrite it
    
    vector <double> phi_tot(Nr);
   
    
    //integrate concentrations and normalize by volume
    for (int j=0;j<ChainType;j++){
        phi_tot[j] = integratedV(phi[j],0,Nr);
        phi_tot[j] /= volume;
    }
    
    
    
    double total=0.0;
    
    for (int j=0;j<ChainType;j++){
        total += phi_tot[j];
    }

    //output total average concentration, if ya want!
    //cout<<total<<endl;
    

}

/**************************************Multiply qdag vectors *********************************************/


/***************Here I calculate the various concentration profiles from the propagators******************/
void phi_calc(Matrix &phi,Matrix qA1,Matrix qdagA1,Matrix qB1,Matrix qdagB1,Matrix qA2,Matrix qB2,Matrix qA3,Matrix qC,vector <int> Ns,vector <double> mu){

    
    for(int i=0;i<Nr;i++){
        
        //Empty array elements
        phi[0][i]=0.0;
        phi[1][i]=0.0;
        phi[2][i]=0.0;
        phi[3][i]=0.0;
        phi[4][i]=0.0;
        phi[5][i]=0.0;
        
        
        //phiA1 integration
        for(int s=0;s<(int)Ns[0]+1;s++){
            if(s==0 || s==(int)Ns[0]){
                phi[0][i]+=0.5*qA1[i][s]*qdagA1[i][Ns[0]-s]*ds;
            }
            else{
                phi[0][i]+=qA1[i][s]*qdagA1[i][Ns[0]-s]*ds;
            }
        }
            
        //phiB1 integration
        for(int s=0;s<(int)Ns[1]+1;s++){
            if(s==0 || s==(int)Ns[1]){
                phi[1][i]+=0.5*qB1[i][s]*qdagB1[i][Ns[1]-s]*ds;
            }
            else{
                phi[1][i]+=qB1[i][s]*qdagB1[i][Ns[1]-s]*ds;
            }
        }
        
        //phiA2 integration
        for(int s=0;s<(int)Ns[2]+1;s++){
            if(s==0 || s==(int)Ns[2]){
                phi[2][i]+=0.5*qA2[i][s]*qA3[i][Ns[2]-s]*ds;
            }
            else{
                phi[2][i]+=qA2[i][s]*qA3[i][Ns[2]-s]*ds;
            }
        }
        
        //phiB2 integration
        for(int s=0;s<(int)Ns[3]+1;s++){
            if(s==0 || s==(int)Ns[3]){
                phi[3][i]+=0.5*qB2[i][s]*qB2[i][Ns[3]-s]*ds;
            }
            else{
                phi[3][i]+=qB2[i][s]*qB2[i][Ns[3]-s]*ds;
            }
        }
        
        //phiA3 integration
        for(int s=0;s<(int)Ns[4]+1;s++){
            if(s==0 || s==(int)Ns[4]){
                phi[4][i]+=0.5*qA3[i][s]*qA2[i][Ns[4]-s]*ds;
            }
            else{
                phi[4][i]+=qA3[i][s]*qA2[i][Ns[4]-s]*ds;
            }
        }
        
        
            //phiC integration
            for(int s=0;s<(int)Ns[5]+1;s++){
                if(s==0 || s==(int)Ns[5]){
                    phi[5][i]+=0.5*qC[i][s]*qC[i][Ns[5]-s]*ds;
                }
                else{
                    phi[5][i]+=qC[i][s]*qC[i][Ns[5]-s]*ds;
                }
            }
            
            //Grand canonical relation
            phi[0][i]=exp(mu[0])*phi[0][i];
            phi[1][i]=exp(mu[1])*phi[1][i];
            phi[2][i]=exp(kappa_ABA*mu[2])*phi[2][i]/kappa_ABA;
            phi[3][i]=exp(kappa_ABA*mu[3])*phi[3][i]/kappa_ABA;
            phi[4][i]=exp(kappa_ABA*mu[4])*phi[4][i]/kappa_ABA;
            phi[5][i]=exp((mu[5])*kappaC)*phi[5][i]*(1.0/kappaC);
    }
    
    
}

/************Here I calculate the diblock/triblock order parameter************/
double calcOP(Matrix phi, double volume){
    
    //define average concentrations
    double phi_ABA;
    double phi_AB;
    double OP;
    
    vector <double> phi_tot(ChainType);

    //integrate concentrations and normalize by volume
    for (int j=0;j<ChainType;j++){
        phi_tot[j] = integratedV(phi[j],0,Nr);
        phi_tot[j] /= volume;
    }
    
    phi_AB = phi_tot[0]+phi_tot[1];
    phi_ABA = phi_tot[0]+phi_tot[1]+phi_tot[2];
    OP = (phi_AB-phi_ABA)/(phi_AB+phi_ABA);
    
    
    return OP;
    
}

//calculate centre hydrophobic maximum
int mmbcentre(Matrix phi){
    int imax;
    double phiB1B2,phiB1B2new;
    imax=0;
    phiB1B2=phi[1][0]+phi[3][0];
    
    for (int i=0;i<Nr;i++){
        phiB1B2new=phi[1][i]+phi[3][i];
        
        if (phiB1B2new>phiB1B2){
            imax=i;
            phiB1B2=phiB1B2new;
        }
    }
    
    return imax;
    
}

//calculate right hydrophilic maximum
int mmbright(Matrix phi,int imax){
    int iright=imax;
    double phiA1A2A3,phiA1A2A3new;
    
    phiA1A2A3 = phi[0][imax]+phi[2][imax]+phi[4][imax];
    
    for (int i=imax;i<Nr;i++){
        phiA1A2A3new = phi[0][i]+phi[2][i]+phi[4][i];
        
        if (phiA1A2A3new>phiA1A2A3){
            iright = i;
            phiA1A2A3=phiA1A2A3new;
        }
    }
    
    return iright;
    
}

//calculate left hydrophobic maximum
int mmbleft(Matrix phi,int imax){
    int ileft=imax;
    double phiA1A2A3,phiA1A2A3new;
    
    phiA1A2A3 = phi[0][imax]+phi[2][imax]+phi[4][imax];
    
    for (int i=imax;i>=0;i--){
        phiA1A2A3new = phi[0][i]+phi[2][i]+phi[4][i];
        
        if (phiA1A2A3new>phiA1A2A3){
            ileft = i;
            phiA1A2A3=phiA1A2A3new;
        }
    }
    
    return ileft;
    
}

//calculate the midpoint between the two locations where phiA = phi B
int mmb_half(Matrix phi, int imax, int pin){
    
    //int pin is the pinning location, which is predetermined
    
    int outer_intersection = 0;
    double del_phi=1.0;
    double del_phi_new = 1.0;
    
    for (int i=imax;i<Nr;i++){
        
        del_phi_new = phi[0][i] + phi[2][i]+ phi[4][i] - phi[1][i] - phi[3][i];
        
        if (del_phi_new<del_phi){
            del_phi = del_phi_new;
            outer_intersection = i;
        }
            
    }
    
    double result = ((double)outer_intersection+(double)pin)/2.0;
    
    int output = result;
    
    return output;
    
}