/*********Here I include an incompressibility function and a pinning function*******/

void Incomp(vector <double> &eta, Matrix phi,vector  <double> &delphi){
    
    double  ptot=0.0; //total concentration
    
    //Calculate incompressibility at each point
    for(int i=0;i<Nr;i++){
        
        ptot=0.0;
        delphi[i]=0.0;
                
        for(int chain=0;chain<ChainType;chain++){
            ptot+=phi[chain][i];
        }
                            
        delphi[i]=1.0-ptot;
        eta[i]-=delphi[i];
        
        //if concentration is too high, something has gone wrong
        if (fabs(delphi[i])>1e3){
            cout<<i<<" incomp: "<<delphi[i]<<endl;
            exit(EXIT_FAILURE);
        }
    }
}


void Pin(vector <double> &sigma, Matrix phi, int pin_location){
    
    //Implement pinning condition
    for (int i=0;i<Nr;i++){
        if (i==pin_location){
            sigma[i]=sigma[i]-10.0*(phi[0][i]+phi[2][i]+phi[4][i]-phi[1][i]-phi[3][i]);
        }
        else{
            sigma[i]=0.0;
        }
    }
}