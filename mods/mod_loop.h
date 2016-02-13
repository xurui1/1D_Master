void mod_loop(double *f,double *mu,double **chiMatrix,double **w,double **phi,double *eta,int *Ns, double *chi, int nradii){
    
    double *r_0vector=create_1d_double_array(nradii+1,"rO_vector");    //Box radius
    double *mu_calc=create_1d_double_array(5,"mu_calc");
    double *rad_test=create_1d_double_array(5,"rad_test");
    double *fA=create_1d_double_array(5,"fA");
    double *A=create_1d_double_array(1,"A");
    double *B=create_1d_double_array(1,"B");
    double *C=create_1d_double_array(1,"C");

    
    double fE_hom;
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
        
        double volume=vol();                                 //calculate volume
        
        set_radius(r_0vector,nradii);
        
        for (int radius=0;radius<nradii;radius++){
            radius = r_0vector[radius];
            double volume=vol();
            omega(w);
            double pin_location=10.8-A[0]-(B[0]*f[0])-(C[0]*f[0]*f[0]);

            displacer=FreeEnergy(w,phi,eta,Ns,chi,chiMatrix,mu,f,pin_location,1);
            
           
        }
       
        
        counter++;
        
    }
    
    
    
    
}
