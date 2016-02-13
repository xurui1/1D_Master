/*************This is my function for calulating the looping fraction for triblock************************/

void calcloop(double **qA2,double **qA2LL,double **qB2LL,double **qA2LR,double **qB2LR,double **qA3,int *Ns,double **w,double *mu,int imax,double *loop){
    
    double Q_ABA;

    //reset looping fractions
    loop[0]=0.0;
    loop[1]=0.0;

    imax = Nr/2;
    
    //Generate constrained case
    for (int i=0;i<Nr;i++){
        if (i<imax){
            qA2LL[i][Ns[0]] = qA2[i][Ns[2]];
            qA2LR[i][Ns[0]] = 0.0;
        }
        else{
            qA2LL[i][Ns[0]] = 0.0;
            qA2LR[i][Ns[0]] = qA2[i][Ns[2]];
        }
    }
    
    
    //solve diffusion equation
    for (int i=0;i<Nr;i++){
        qB2LL[i][0]=qA2LL[i][Ns[0]];
        qB2LR[i][0]=qA2LR[i][Ns[0]];
    }
    solvediffyQ(qB2LL,w[3],Ns[3]);
    solvediffyQ(qB2LR,w[3],Ns[3]);
    
    //Calculate ABA chain partition function
    Q_ABA = integrate2d(qA3,0,Nr,Ns[4]);       //integrate
    Q_ABA = exp(mu[2]*kappa_ABA)*Q_ABA/kappa_ABA;         //G-C relation
    
    
    
    //Calculate probability of looping
    for (int i=0;i<Nr/2;i++){
        loop[0]+=(qB2LL[i][Ns[3]])*(qA2[i][Ns[2]])/(Q_ABA);
    }
    for (int i=Nr/2;i<Nr;i++){
        loop[1]+=(qB2LR[i][Ns[3]])*(qA2[i][Ns[2]])/(Q_ABA);
    }

}

