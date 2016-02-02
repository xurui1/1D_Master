double Conc(Matrix &phi,Matrix &w,vector <int> &Ns,vector <double> &mu, double volume, vector <double> &loop,  int out_loop, int iter,int pin_location){
    
    double      Q;
    
    //Forwards propagators
    Matrix qA1(Nr, Row(Ns[0]+1));
    Matrix qB1(Nr, Row(Ns[1]+1));
    Matrix qA2(Nr, Row(Ns[2]+1));
    Matrix qB2(Nr, Row(Ns[3]+1));
    Matrix qA3(Nr, Row(Ns[4]+1));
    Matrix qC(Nr, Row(Ns[5]+1));
    
    //Complementary propagators
    Matrix qdagA1(Nr, Row(Ns[0]+1));
    Matrix qdagB1(Nr, Row(Ns[1]+1));
    
    //Looping propagators
    Matrix qA2LoopLeft(Nr,Row(Ns[2]+1));
    Matrix qB2LoopLeft(Nr,Row(Ns[3]+1));
    Matrix qA2LoopRight(Nr,Row(Ns[2]+1));
    Matrix qB2LoopRight(Nr,Row(Ns[3]+1));

    //solve diffusion equations
    diblock(qA1,qdagA1,qB1,qdagB1,w,Ns);
    triblock(qA2,qB2,qA3,w,Ns);
    homopolymer(qC,w,Ns);
    
    // Here we get the single chain partition functions Q_AB+Q_C
    Q=q_partition(qB1,qA3,qC,Ns,mu,volume);
        
    //cout<<"Q: "<< Q<<endl;
    
    // Here we do the concentration calculation by integration over chain
    phi_calc(phi,qA1,qdagA1,qB1,qdagB1,qA2,qB2,qA3,qC,Ns,mu);

    //calculation of average concentrations over entire computation box
    phi_total(phi,volume);
    
    //find max conc
    int imax = mmbcentre(phi);
    int ihalf = mmb_half(phi,imax,pin_location);
    
    if (out_loop==1 && iter%100==0 ){
        //calculate looping fraction
        calcloop(qA2,qA2LoopLeft,qB2LoopLeft,qA2LoopRight,qB2LoopRight,qA3,Ns,w,mu,ihalf,loop);
    }
    

    return Q;
    
}