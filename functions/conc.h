double Conc(double **phi,double **w,int *Ns,double *mu, double volume, double *loop,  int out_loop, int iter,int pin_location, double **qA1,double **qB1,double **qA2,double **qB2,double **qA3,double **qC,double **qdagA1,double **qdagB1,double **qB2LoopLeft,double **qB2LoopRight){
    
    
    //solve diffusion equations for propagators
    diblock(qA1,qdagA1,qB1,qdagB1,w,Ns);
    triblock(qA2,qB2,qA3,w,Ns);
    homopolymer(qC,w,Ns);
    
    // Here we get the single chain partition functions Q_AB+Q_C
    double Q=q_partition(qB1,qA3,qC,Ns,mu,volume);
        
    //cout<<"Q: "<< Q<<endl;
    
    // Here we do the concentration calculation by integration over chain
    phi_calc(phi,qA1,qdagA1,qB1,qdagB1,qA2,qB2,qA3,qC,Ns,mu);

    //calculation of average concentrations over entire computation box
    if (iter%100==0){
        phi_total(phi,volume);
    }
    //find max conc
    int imax = mmbcentre(phi);
    int ihalf = mmb_half(phi,imax,pin_location);
    
    if (out_loop==1 && iter%100==0 ){
        //calculate looping fraction
        calcloop(qA2,qB2LoopLeft,qB2LoopRight,qA3,Ns,w,mu,imax,loop);
    }
    

    return Q;
    
}