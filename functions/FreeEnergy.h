/*************************This is my main function for calculating free energies*************************/

double FreeEnergy(Matrix &w, Matrix &phi, vector <double> &eta, vector <int> &Ns,vector <double> &chi,  Matrix chiMatrix, vector <double> &mu, double volume, vector <double> &f, int pin_location, int out_loop){
    
    
    double  currentfE, oldfE, deltafE;
    int     maxIter=1e6;
    double precision=1e-5;          //convergence condition
    int     mmb;                    //Turns on pinning condition
    double  Q;                      //Chain partition functions
    double  fE_int, fES;            //interaction free energy and chain partition function fE
    double  deltaW;
    double fE_hom = homogfE(mu,chiMatrix,f);

    
    //Arrays for updating the omega fields
    Matrix delW(ChainType,Row(Nr));
    Matrix newW(ChainType,Row(Nr));

    vector <double> delphi(Nr);
    vector <double>sigma(Nr);
    vector <double>loop(2);
    
    //set energies to zero
    currentfE=0.0;
    deltafE=0.0;
    
    
    //Turn pinning condition on for membrane
    mmb=1;

        
    for (int iter=0;iter<maxIter;iter++){
        
        fE_int=0.0;
        fES=0.0;
        deltaW=0.0;

        
        Q=Conc(phi,w,Ns,mu,volume,loop,out_loop,iter, pin_location);      //Calculate Chain partition functions
        
        
        Incomp(eta,phi,delphi);           //Enforce incompressibility condition
        
        if (iter%1000==0){
            output(phi);                   //Output concentration data to file
        }
        
        if (mmb==1){
            Pin(sigma, phi, pin_location);
        }
        
        //Calculate components for new field and interaction free energies
        for(int i=0;i<Nr;i++){
            for(int ii=0;ii<ChainType;ii++){
                
                //reset chem potential update
                newW[ii][i]=0.0;
                
                for(int jj=0;jj<ChainType;jj++){
                    newW[ii][i]+=(chiMatrix[ii][jj]*phi[jj][i]);
                }
                
                //add incompressibility
                newW[ii][i]+=eta[i];
                
                //apply pinning condition
                if (mmb==1){
                    if (ii==0 || ii == 2 || ii==4){
                        newW[ii][i]-=sigma[i];
                    }
                    else if (ii==1 || ii==3){
                        newW[ii][i]+=sigma[i];
                    }
                }
                delW[ii][i]=newW[ii][i]-w[ii][i];                    //change in omega field
                w[ii][i]+=(gamma_up*delW[ii][i]-epsilon_up*delphi[i]);     //update omega field
                deltaW+=fabs(delW[ii][i])*dV(i);                  //total change
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
        if (iter%100==0 && out_loop==1){std::cout<<iter<<" fE:"<<currentfE<< " dfE:"<<currentfE-fE_hom<<" " << deltaW<<" "<<fE_hom<<" "<<loop[0]<<" "<<loop[1]<<std::endl;}
        else if  (iter%100==0 && out_loop==0){std::cout<<iter<<" fE:"<<currentfE<< " dfE:"<<currentfE-fE_hom<<" " << deltaW<<" "<<fE_hom<<std::endl;}

        

        if (deltafE<precision && deltaW<precision){break;} //Convergence condition
        
    }

    //output loop data if certain conditions are met
    if ((f[0] == 0.3 || f[0] == 0.4 || f[0] == 0.5 || f[0] == 0.6 || f[0] == 0.7) && out_loop==1){
        string filename="./results/loop/loopr"+IntToStr((int)10*f[0])+".dat";
        std::ofstream outputloop;
        int imax=mmbcentre(phi);
    
        outputloop.open(filename.c_str(), std::ofstream::app);
        outputloop <<4.3/(r_0+imax*dr)<<" "<<loop[0]<<" "<<loop[1]<<" "<<0.5*(loop[0]+loop[1])<< endl;
        outputloop.close();
    }
    
 
    
    return currentfE-fE_hom;
    
}
