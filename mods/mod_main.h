void mod_main(vector <double> &f,vector <double>&mu,Matrix &chiMatrix,Matrix &w,Matrix &phi,vector <double> &eta,vector <int> &Ns, vector <double> &chi,int nfa,vector <double> &A,vector <double> &B,vector <double> &C, int nradii, vector <double> &dFE, vector <double> &mu_vec){
    
    vector <double> r_0vector(nradii+1);    //Box radius
    vector <double> Rad(nradii);            //Radius for fitting
    vector <double> diameter(nradii);       //membrane diameter
    vector <double> Curv(nradii);           //curvature
    vector <double> Curvsq(nradii);         //squared curvature
    
    Matrix fitting(6,Row(nfa));             //matrix of fitting constants
    
    //open main output file
    ofstream outFile2;
    string filename2;
    filename2="./results/fA_test.dat";
    outFile2.open(filename2.c_str());
    
    //radius ouput
    ofstream radiout;
    radiout.open("./results/main_radius.dat");
    
    
    //output quadratic bending modulus
    ofstream bendingout;
    bendingout.open("./results/bending_mod.dat");
    
    //output average diameter
    ofstream diamout;
    diamout.open("./results/diameter_output.dat");
    
    
    int counter=0;
    double volume;
    double OP;
    double fE_hom;
    
    int ds_increment = 4;
    
    for (int dds=0; f[0]<0.7 ;dds+= ds_increment){
        counter+=1;
        //Set parameters s
        updateparameters(f,Ns,dds);
        mu[5] = mu_vec[counter-1];                        //don't want to calc mu again
        fE_hom=homogfE(mu,chiMatrix,f);                 //calculate homog. fE
        omega(w);                                       //Initiate omega field
        
        double pin_location=10.8-A[0]-(B[0]*f[0])-(C[0]*f[0]*f[0]);
        int pin = pin_location/dr;
        
        volume=vol();                                 //calculate volume
        OP = calcOP(phi,volume);                     //calculate order parameter
        //Set radius vector
        set_radius(r_0vector,nradii);
        
        r_0=r_0vector[0];                                        //reset radius
        double avgradius=0.0;                                  //reset avgradius
        double avgmiddle=0.0;
        double avgdiameter=0.0;
        
        for (int radius=0;radius<nradii;radius++){
            volume=vol();
            
            //initialize omega field
            omega(w);
            
            //calculate free energy minus homogeneneous free energy
            dFE[radius]=FreeEnergy(w,phi,eta,Ns,chi,chiMatrix,mu,volume,f,pin,1);
            OP = calcOP(phi,volume);                    //calculate order parameter
            diameter[radius] = calc_excess(phi,volume); //calculate copolymer excess
            avgdiameter+=diameter[radius];
            
            int imax=mmbcentre(phi);                       //membrane center (max phib)
            int ihalf=mmb_half(phi,imax,pin);                      //membrane middle (1/2 phib = phiA)
            
            avgradius+=(double)imax*dr;                    //avg membrane center
            avgmiddle+=(double)ihalf*dr;
            
            Rad[radius]=r_0;                               //set radius vectors
            
            //output free energy data
            outFile2 <<f[0]<<" "<< r_0 << " "<<r_0+(double)imax*dr<<" "<<dFE[radius]<<std::endl;
            cout<<f[0]<<" "<< r_0 << " "<<r_0+(double)imax*dr<<" "<<dFE[radius]<<std::endl;
            
            //output concentration profile
            outputphi(phi);
            
            //set new radius
            r_0=r_0vector[radius+1];
            
        }
        //output average radius
        avgradius/= (double)nradii;
        avgmiddle/= (double)nradii;
        radiout<<f[0]<<" "<<avgradius<<" "<<avgmiddle<<endl;
        
        //calculate average diameter
        avgdiameter/=nradii;
        diamout<<f[0]<<" "<<avgdiameter<<endl;
        
        
        //build curvature and curvature squared vectors
        for (int radius=0;radius<nradii;radius++){
    
            //membrane should be centered
            Rad[radius]+=6.0;

            //mmb thickness approx 4.3 for consistency
            Curv[radius] =(4.3/Rad[radius]);
            Curvsq[radius] = pow(Curv[radius],2.0);
    
        }
        
        //output files of free energy as a function of radius
        outputfE_FA(f[0],Curv,dFE,nradii);
        
        //quadratic curve fit with max hydrophobic
        curvefit(Curv,dFE,nradii,counter,fitting[0],fitting[1],fitting[2]);
        
        //quartic curve fit with max hydrophobic
        curvefit(Curvsq,dFE,nradii,counter,fitting[3],fitting[4],fitting[5]);
        
        
        if (Coord==2){
            bendingout<<f[0]<<" "<<fitting[4][counter-1]<<" "<<2.0*fitting[4][counter-1]*(double)Ds/(sqrt(chi[0])*pow(4.3,2.0))<<endl;
        }
        else if (Coord==3){
            bendingout<<f[0]<<" "<<fitting[4][counter-1]<<" "<<fitting[4][counter-1]*(double)Ds/(sqrt(chi[0])*pow(4.3,2.0))<<endl;
        }
        
    }
    outFile2.close();
    radiout.close();
    bendingout.close();
    diamout.close();
    
    //output all curvefit results
    outputkappa(fitting[0],fitting[1],fitting[2],fitting[3],fitting[4],fitting[5],nfa,chi,0);
    
    
}