void mod_loop(vector <double> &f,vector <double> &mu,Matrix &chiMatrix,Matrix &w,Matrix &phi,vector <double> &eta,vector <int> &Ns, vector <double> &chi, int nradii){
    
    vector <double> r_0vector(nradii+1);    //Box radius
    vector <double> mu_calc(5);
    vector <double>  rad_test(5);
    vector <double> fA(5);
    vector <double> A(1);
    vector <double> B(1);
    vector <double> C(1);

    
    double fE_hom;
    double volume;
    double displacer;
    
    int counter=0;
    int ds_increment = 20;
    
    
    for (int dds=0;f[0]<=0.7;dds+=ds_increment){
        double avgmiddle=0.0,avghalf=0.0;

        updateparameters(f,Ns,dds);
        fA[counter] = f[0];
        fE_hom=homogfE(mu,chiMatrix,f);                 //calculate homog. fE
        omega(w);                                       //Initiate omega field

        secant(w,phi,eta,Ns,chi,chiMatrix,mu,f,2*Nr/5);  //Find tensionless mmb
        mu_calc[counter]=mu[5];
        r_0 = 1.0;
        for (int radius=0;radius<4;radius++){
            volume=vol();
            omega(w);
            int imax=mmbcentre(phi);
            int ihalf=mmb_half(phi,imax,2*Nr/5);
            avgmiddle+=(double)imax*dr;
            avghalf+=(double)ihalf*dr;
            
            r_0*=3.0;
        }
        avgmiddle/=4.0;
        avghalf/=4.0;
        rad_test[counter] = avghalf;
        counter++;
        
    }
    curvefit(fA,rad_test,5,1,A,B,C);
    
    counter =0;
    for (int dds=0; f[0]<=0.7 ;dds+= ds_increment){
        //Set parameters
        updateparameters(f,Ns,dds);
        fE_hom=homogfE(mu,chiMatrix,f);                 //calculate homog. fE
        omega(w);                                       //Initiate omega field
        mu[5] = mu_calc[counter];
        
        volume=vol();                                 //calculate volume
        
        set_radius(r_0vector,nradii);
        
        for (int radius=0;radius<nradii;radius++){
            radius = r_0vector[radius];
            volume=vol();
            omega(w);
            double pin_location=10.8-A[0]-(B[0]*f[0])-(C[0]*f[0]*f[0]);

            displacer=FreeEnergy(w,phi,eta,Ns,chi,chiMatrix,mu,volume,f,pin_location,1);
            
           
        }
       
        
        counter++;
        
    }
    
    
    
    
}
