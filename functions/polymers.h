//Here I solve the propagators for the diblock and the triblock and the homopolymer

void diblock(double **qA1,double **qdagA1,double **qB1,double **qdagB1,double **w,int *Ns){
    
    //empty propagators
    for (int i=0;i<Nr;i++){
        for (int j=0;j<=Ns[0];j++){
            qA1[i][j]=0.0;
            qdagA1[i][j]=0.0;
        }
        for (int j=0;j<=Ns[1];j++){
            qB1[i][j]=0.0;
            qdagB1[i][j]=0.0;
        }
    }
    
    // Here is the for loop for setting the A propagator initial conditions to 1.0
    for(int i=0;i<Nr;i++){
        qA1[i][0]=1.0;
    }
    
    // Here we solve the diffusion equation for the A forwards propagator
    solvediffyQ(qA1,w[0],Ns[0]);
    
    //Initial condition for B forwards propagator from A propagator
    for(int i=0;i<Nr;i++){
        qB1[i][0]=qA1[i][Ns[0]];
    }
    
    //Here we solve the diffusion equation for the B forwards propagator
    solvediffyQ(qB1,w[1],Ns[1]);
    
    //Set initial condition for the B Complementary propagator
    for(int i=0;i<Nr;i++){
        qdagB1[i][0]=1.0;
    }
    
    //Here we solve the diffusion for the B complementary propagator
    solvediffyQ(qdagB1,w[1],Ns[1]);
    
    //Initial condition for A complementary propagator from B comp
    for(int i=0;i<Nr;i++){
        qdagA1[i][0]=qdagB1[i][Ns[1]];
    }
    
    //Here we solve the diffusion equation for the A complementary propagator
    solvediffyQ(qdagA1,w[0],Ns[0]);
    
}

void triblock(double **qA2,double **qB2,double **qA3,double **w,int *Ns){
    //Triblock is symmetric so we don't need complementary propagator
    
    //empty propagators
    for (int i=0;i<Nr;i++){
        for (int j=0;j<=Ns[2];j++){
            qA2[i][j]=0.0;
            qA3[i][j]=0.0;
        }
        for (int j=0;j<=Ns[3];j++){
            qB2[i][j]=0.0;
        }
    }
    //Set initial condition for A2 forward propagator
    for (int i=0;i<Nr;i++){
        qA2[i][0]=1.0;
    }
    
    //Solve diffusion equation for A2 forward propagator
    solvediffyQ(qA2,w[2],Ns[2]);
    
    //Initial condition for B2 forward propagator
    for (int i=0;i<Nr;i++){
        qB2[i][0]=qA2[i][Ns[2]];
    }
    
    //Solve diffusion equation for B2 forward propagator
    solvediffyQ(qB2,w[3],Ns[3]);
    
    //Set initial condition for A3 forward propagator
    for (int i=0;i<Nr;i++){
        qA3[i][0]=qB2[i][Ns[3]];
    }
    
    //Solve diffusion equation for A3 forward propagator
    solvediffyQ(qA3,w[4],Ns[4]);
        
}

void homopolymer(double **qC,double **w,int *Ns){
    
    for (int i=0;i<Nr;i++){
        for (int j=0;j<=Ns[5];j++){
            qC[i][j]=0.0;
        }
    }
    
    //Set initial condition for homopolymer propagator
    for (int i=0;i<Nr;i++){
        qC[i][0]=1.0;
    }
    
    //Solve diffusion equation for homopolymer
    solvediffyQ(qC,w[5],Ns[5]);
    
}

