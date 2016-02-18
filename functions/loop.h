/*************This is my function for calulating the looping fraction for triblock************************/

void calcloop(double **qA2,double **qB2LL,double **qB2LR,double **qA3,int *Ns,double **w,double *mu,int imax,double *loop){
    
    double Q_ABA;

    //reset looping fractions
    loop[0]=0.0;
    loop[1]=0.0;

    //imax = Nr/2;
    
    //Generate constrained case
    for (int i=0;i<Nr;i++){
        if (i<imax){
            qB2LL[i][0] = qA2[i][Ns[2]];
            qB2LR[i][0] = 0.0;
        }
        else{
            qB2LL[i][0] = 0.0;
            qB2LR[i][0] = qA2[i][Ns[2]];
        }
    }
    
    solvediffyQ(qB2LL,w[3],Ns[3]);
    solvediffyQ(qB2LR,w[3],Ns[3]);
    
    //Calculate ABA chain partition function
    Q_ABA=integrate2d_dV(qA3,0,Nr,Ns[4]);       //integrate
    Q_ABA = exp(mu[2]*kappa_ABA)*Q_ABA/(kappa_ABA);         //G-C relation
    cout<<Q_ABA<<endl;
    

    //Calculate probability of looping
    for (int i=0;i<=imax;i++){
        double term =((qB2LL[i][Ns[3]])*(qA2[i][Ns[2]])*dV(i))/Q_ABA;
        loop[0]+=term;
    }
    for (int i=imax;i<Nr;i++){
        double term = ((qB2LR[i][Ns[3]])*(qA2[i][Ns[2]])*dV(i))/Q_ABA;
        loop[1]+=term;
    }
    double unconfined=0.0;
    for (int i=0;i<Nr;i++){
        unconfined += qA3[i][0]*qA2[i][Ns[2]]*dV(i)/Q_ABA;
    }
    
    loop[0]/=unconfined;
    loop[1]/=unconfined;

}

