/**********Here I calculate three chain partition functions by integrating over space****************/

double q_partition(double **qB1,double **qA3,double **qC,int *Ns,double *mu, double volume){
    
    double Q,Q_AB,Q_C,Q_ABA;
 
    Q_AB=integrate2d_dV(qB1,0,Nr,Ns[1]);
    Q_C=integrate2d_dV(qC,0,Nr,Ns[5]);
    Q_ABA=integrate2d_dV(qA3,0,Nr,Ns[4]);
    
    Q_AB=exp(mu[0])*Q_AB;
    Q_ABA=exp(mu[2]*kappa_ABA)*Q_ABA/kappa_ABA;
    Q_C=(exp(mu[5]*kappaC)*Q_C)/kappaC;
    
    //I'm adding the three single chain partition functions together for the return function
    Q=Q_AB+Q_C+Q_ABA;
    
    // Normalizing with respect to box volume
    Q/=volume;
    
    
    return Q;
}