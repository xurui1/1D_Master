/****************************Here I set some parameters *****************************/

void parameters(vector <double> &chi,vector <double> &f, vector <int> &Ns, vector <double> &mu){
    
    //inner box radius
    r_0=1.0;
    
    //initial settings
    initial=3;
    
    double chi_input;
    
    ifstream inputchi;
    inputchi.open("chi.dat");
    inputchi >> chi_input;
    inputchi.close();
    
    //chi_input=30.0;
        
    //Interaction parameters
    chi[0]=chi_input;   //Chi_AB
    chi[1]=chi_input;   //Chi_BC
    chi[2]=0.0;         //Chi_AC
    
    cout<<"Chi input:"<<chi_input<<endl;
    
    //Type of system
    //poly=0;
    
    //Chemical potential array
    if (poly==0){
        mu[0]=-50.0;    //AB -A
        mu[1]=mu[0];    //AB -B
        mu[2]=0.0;      //ABA -A
        mu[3]=mu[2];    //ABA -B
        mu[4]=mu[2];    //ABA -A
        mu[5]=-5.0;    //C
    }
    else if (poly==1){
        mu[0]=0.0;    //AB -A
        mu[1]=mu[0];    //AB -B
        mu[2]=-50.0;      //ABA -A
        mu[3]=mu[2];    //ABA -B
        mu[4]=mu[2];    //ABA -A
        mu[5]=-5.0;    //C
    }
    else if (poly==2){
        mu[0]=0.0;      //Other
        mu[1]=0.0;
        mu[2]=0.0;
        mu[3]=0.0;
        mu[4]=0.0;
        mu[5]=0.0;
    }
    
    //Chain length array
    Ns[0]=60;            //A1 blocks
    Ns[1]=140;            //B1 blocks
    Ns[2]=Ns[0];        //A2 block
    Ns[3]=2*Ns[1];      //B2 blocks
    Ns[4]=Ns[0];        //A3 blocks
    Ns[5]=Ds;            //C blocks
    
    double diblock = (double)Ns[0]+(double)Ns[1];
    double triblock= (double)Ns[2]+(double)Ns[3]+(double)Ns[4];
    double homopoly = (double)Ns[5];
    
    //Chain fraction array
    f[0]=(double)Ns[0]/diblock;      //A1
    f[1]=(double)Ns[1]/diblock;      //B1
    f[2]=(double)Ns[2]/triblock;
    f[3]=(double)Ns[3]/triblock;
    f[4]=(double)Ns[4]/triblock;
    f[5]=(double)Ns[5]/homopoly;
    
    //Length ratio of c homopolymer to diblock copolymer
    kappaC=homopoly/diblock;
    kappa_ABA= triblock/diblock;
    
    
    //cout<<Ns[0]<<" "<<Ns[1]<<" "<<Ns[2]<<endl;
    
    //Step size in r,z direction
    dr=12.0/(double)Nr;
    
    //Step length along polymer
    ds=1.0/(double)Ds;
    
}


/****************Here I build the interaction Matrix************/

void Xmatrix(Matrix &chiMatrix, vector <double> &chi){
    //Interaction Matrix
    chiMatrix[0][0]=0.0;    //ChiA1,A1
    chiMatrix[0][1]=chi[0]; //ChiA1,B1
    chiMatrix[0][2]=0.0;    //ChiA1,A2
    chiMatrix[0][3]=chi[0]; //ChiA1,B2
    chiMatrix[0][4]=0.0;    //ChiA1,A3
    chiMatrix[0][5]=chi[2]; //ChiA1,C
    
    chiMatrix[1][0]=chi[0]; //ChiB1,A1
    chiMatrix[1][1]=0.0;    //ChiB1,B1
    chiMatrix[1][2]=chi[0]; //ChiB1,A2
    chiMatrix[1][3]=0.0;    //ChiB1,B2
    chiMatrix[1][4]=chi[0]; //ChiB1,A3
    chiMatrix[1][5]=chi[1]; //ChiB1,C
    
    chiMatrix[2][0]=0.0;    //ChiA2,A1
    chiMatrix[2][1]=chi[0]; //ChiA2,B1
    chiMatrix[2][2]=0.0;    //ChiA2,A2
    chiMatrix[2][3]=chi[0]; //ChiA2,B2
    chiMatrix[2][4]=0.0;    //ChiA2,A3
    chiMatrix[2][5]=chi[2]; //ChiA2,C
    
    chiMatrix[3][0]=chi[0]; //ChiB2,A1
    chiMatrix[3][1]=0.0;    //ChiB2,B1
    chiMatrix[3][2]=chi[0]; //ChiB2,A2
    chiMatrix[3][3]=0.0;    //ChiB2,B2
    chiMatrix[3][4]=chi[0]; //ChiB2,A3
    chiMatrix[3][5]=chi[1]; //ChiB2,C
    
    chiMatrix[4][0]=0.0;    //ChiA3,A1
    chiMatrix[4][1]=chi[0]; //ChiA3,B1
    chiMatrix[4][2]=0.0;    //ChiA3,A2
    chiMatrix[4][3]=chi[0]; //ChiA3,B2
    chiMatrix[4][4]=0.0;    //ChiA3,A3
    chiMatrix[4][5]=chi[2]; //ChiA3,C
    
    
    chiMatrix[5][0]=chi[2]; //ChiC,A1
    chiMatrix[5][1]=chi[1]; //ChiC,B1
    chiMatrix[5][2]=chi[2]; //ChiC,A2
    chiMatrix[5][3]=chi[1]; //ChiC,B2
    chiMatrix[5][4]=chi[2]; //ChiC,A3
    chiMatrix[5][2]=0.0;    //ChiC,C
    
}

/******Here I update parameters for looping through fA*******/

void updateparameters(vector <double> &f, vector <int> &Ns, int dds){
    
    //Chain length array
    Ns[0]=60+dds;
    Ns[1]=140-dds;
    Ns[2]=60+dds;
    Ns[3]=2*(140-dds);
    Ns[4]=60+dds;
    Ns[5]=Ds;
    
    //Chain fraction array
    f[0]=(double)Ns[0]/(double)Ds;    //A
    f[1]=(double)Ns[1]/(double)Ds;  //B
    f[2]=f[0];
    f[3]=f[1];
    f[4]=f[0];

    
}
