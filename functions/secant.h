/***************Here I find the chemical potential for tensionless membranes****************************/

void secant(double **w,double **phi,double *eta,int *Ns,double *chi,double **chiMatrix,double *mu,double *f, int pin_location){
    
    
    double  currentfE, oldfE, deltafE;
    int     maxIter=1e6;
    double precision=1.0e-5;          //convergence condition
    double mu1,mu2,mu3;
    double fE1,fE2,fE3;
    int     mmb;
    int coordinate;
    double  Q;                      //total partition function (3 chains= AB, ABA, C)
    double  fE_int, fES;            //interaction free energy and chain partition function fE
    double  deltaW;
    double fE_hom=0.0;
    double volume;
    
    //Arrays for updating the omega fields
    double **delW=create_2d_double_array(ChainType,Nr,"delW");
    double **newW=create_2d_double_array(ChainType,Nr,"newW");
    
    //allocate normal propagators
    double **qA1=create_2d_double_array(Nr,Ns[0]+1,"qA1");
    double **qB1=create_2d_double_array(Nr,Ns[1]+1,"qB1");
    double **qA2=create_2d_double_array(Nr,Ns[2]+1,"qA1");
    double **qB2=create_2d_double_array(Nr,Ns[3]+1,"qB1");
    double **qA3=create_2d_double_array(Nr,Ns[4]+1,"qB1");
    double **qC=create_2d_double_array(Nr,Ns[5]+1,"qB1");
    
    //allocate complementary propagators for diblock
    double **qdagA1=create_2d_double_array(Nr,Ns[0]+1,"qdagA1");
    double **qdagB1=create_2d_double_array(Nr,Ns[1]+1,"qdagB1");
    
    
    //allocate looping propagators
    double **qB2LoopLeft=create_2d_double_array(Nr,Ns[3]+1,"qB1LoopLeft");
    double **qB2LoopRight=create_2d_double_array(Nr,Ns[3]+1,"qB1LoopRight");
    
    
    double *delphi=create_1d_double_array(Nr,"delphi");
    double *sigma = create_1d_double_array(Nr, "sigma");
    double *loop = create_1d_double_array(2,"loop");
    
    currentfE=0.0;
    deltafE=0.0;
    
    
    //we are adjusting mu_C to get tensionless membrane
    mu1=mu[5];
    mu2=mu[5]+0.01;
    
    
    mmb=1;
    coordinate=Coord;
    Coord=1;
    
    for (int xx=0;xx<3;xx++){
        fE_hom=homogfE(mu,chiMatrix,f);
        volume=vol();

    for (int iter=0;iter<maxIter;iter++){
        
        fE_int=0.0;
        fES=0.0;
        deltaW=0.0;
        
        
        Q=Conc(phi,w,Ns,mu,volume,loop,0,iter, pin_location,qA1,qB1,qA2,qB2,qA3,qC,qdagA1,qdagB1,qB2LoopLeft,qB2LoopRight);       //Calculate Chain partition functions
        
        
        Incomp(eta,phi,delphi);              //Enforce incompressibility condition
        
        if (iter%100==0) {
            output(phi);                   //Output some data to file
        }
        
        if (mmb==1){Pin(sigma, phi,pin_location);}
        
        
        
        //Calculate components for new field and interaction free energies
        for(int i=0;i<Nr;i++){
            for(int ii=0;ii<ChainType;ii++){
                newW[ii][i]=0.0;            //set field update to zero
                for(int jj=0;jj<ChainType;jj++){
                    newW[ii][i]+=(chiMatrix[ii][jj]*phi[jj][i]);
                }
                newW[ii][i]+=eta[i];
                
                if (mmb==1){
                    if (ii==0 || ii == 2 || ii==4){
                        newW[ii][i]-=sigma[i];
                    }
                    else if (ii==1 || ii==3){
                        newW[ii][i]+=sigma[i];
                    }
                }
                delW[ii][i]=newW[ii][i]-w[ii][i];
                w[ii][i]+=(gamma_up*delW[ii][i]-epsilon_up*delphi[i]);     //update omega field
                deltaW+=fabs(delW[ii][i]);
            }
        }
        fE_int=fE(newW,phi,chiMatrix,volume);
        
        //Normalize by box size
        deltaW/=volume;
        
        
        //Update free energy
        fES=Q;
        oldfE=currentfE;
        currentfE=-fES+fE_int;
        deltafE=fabs(currentfE-oldfE);
        
        //Print free energy, difference in free energy, change in omega field to screen
        if (iter%100==0){cout<<xx<<" "<<iter<<" fE:"<<currentfE<< " dfE:"<<currentfE-fE_hom<<" " << deltaW<<" "<<fE_hom<<endl;}
        
        
        if (deltafE<precision && deltaW<precision){break;} //Convergence condition
        
    }
    
        
        //Secant method for finding chemical potential
        if (xx==0){
            fE1=currentfE-fE_hom;
            mu[5]=mu2;
        }
        else if(xx==1){
            fE2=currentfE-fE_hom;
            mu3=mu2-(fE2*((mu2-mu1)/(fE2-fE1)));
            mu[5]=mu3;
        }
        else if (xx==2){
            fE3=currentfE-fE_hom;
            if (abs(fE3)>precision){
                mu1=mu2;
                mu2=mu3;
                mu[5]=mu1;
                xx=-1;
            }
        }
        cout<<"muC: "<<mu[5]<<endl;
        if (xx>2){cout<<"You messed up secant mod"<<endl;}
    }

 //   mu[2]=mu3;
    cout<<"muC: "<<mu[5]<<" done secant mod"<<endl;
    Coord=coordinate;
    
    ofstream outfile;
    outfile.open("./results/mu.dat");
    outfile << mu[5];
    outfile.close();

    destroy_1d_double_array(loop);
    destroy_2d_double_array(newW);
    destroy_1d_double_array(sigma);
    destroy_1d_double_array(delphi);
    destroy_2d_double_array(qA1);
    destroy_2d_double_array(qA2);
    destroy_2d_double_array(qA3);
    destroy_2d_double_array(qB1);
    destroy_2d_double_array(qB2);
    destroy_2d_double_array(qC);
    destroy_2d_double_array(qB2LoopLeft);
    destroy_2d_double_array(qB2LoopRight);
    
}
