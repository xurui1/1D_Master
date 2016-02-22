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
    int ds_increment = 0.1*Ds;
    
    
    for (int dds=0;f[0]<=0.7;dds+=ds_increment){
        double avgmiddle=0.0,avghalf=0.0;

        updateparameters(f,Ns,dds);
        fA[counter] = f[0];
        fE_hom=homogfE(mu,chiMatrix,f);                 //calculate homog. fE
        omega(w);                                       //Initiate omega field

        secant(w,phi,eta,Ns,chi,chiMatrix,mu,f,2*Nr/5);  //Find tensionless mmb
        mu_calc[counter]=mu[5];
        r_0 = 1.0;
        int r_max = 5;
        for (int radius=0;radius<r_max;radius++){
            omega(w);
            displacer=FreeEnergy(w,phi,eta,Ns,chi,chiMatrix,mu,f,2*Nr/5,1);

            int imax=mmbcentre(phi);
            int ihalf=mmb_half(phi,imax,2*Nr/5);
            avgmiddle+=(double)imax*dr;
            avghalf+=(double)ihalf*dr;
            
            r_0*=3.0;
        }
        avgmiddle/=r_max;
        avghalf/=r_max;
        rad_test[counter] = avghalf;
        counter++;
        
    }
    curvefit(fA,rad_test,5,1,A,B,C);
    
    
    //reset parameters
    parameters(chi,f,Ns,mu);
    counter =0;

    
    for (int dds=0; f[0]<=0.7 ;dds+= ds_increment){
        //Set parameters
        updateparameters(f,Ns,dds);
        fE_hom=homogfE(mu,chiMatrix,f);                 //calculate homog. fE
        omega(w);                                       //Initiate omega field
        mu[5] = mu_calc[counter];
        
        
        set_radius(r_0vector,nradii);
        
        for (int radius=0;radius<nradii;radius++){
            radius = r_0vector[radius];
            omega(w);
            double pin_location=10.8-A[0]-(B[0]*f[0])-(C[0]*f[0]*f[0]);

            displacer=FreeEnergy(w,phi,eta,Ns,chi,chiMatrix,mu,f,pin_location,1);
            
           
        }
       
        
        counter++;
        
    }
    
    
    
    
}
