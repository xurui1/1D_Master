double rterm(int i){
    
    if (Coord==1){
        return 0.0;
    }
    else{
        return ((double)Coord-1.0)*(ds/((dr*((double)i)+r_0)*4.0*dr));
    }
}




//Build LHS matrix
void Matrix_r(double *w,double *rmid,double *rupper,double *rlower){
    int i;
    
    for (i=0;i<Nr;i++){
        rmid[i]=1.0+(ds/(pow(dr,2.0)))+(ds/2.0)*w[i];
        rupper[i]=-(ds/(2.0*pow(dr,2.0)))-rterm(i);
        rlower[i]=-(ds/(2.0*pow(dr,2.0)))+rterm(i);
    }
    
    rupper[0]=rupper[0]+rlower[0];
    rlower[Nr-1]=rupper[Nr-1]+rlower[Nr-1];
    rupper[Nr-1]=0.0;
    rlower[0]=0.0;

}



//Apply Crank-Nicolson method to solve the modified diffusion equation
void solvediffyQ(double **q,double *w, int Ns){

    double gamma, betaL, betaU;
    
    //Allocate memory for `matrix' operations
    
    double *bvecr=create_1d_double_array(Nr,"bvecr");
    double *rupper=create_1d_double_array(Nr,"rupper");
    double *rmid=create_1d_double_array(Nr,"rmid");
    double *rlower=create_1d_double_array(Nr,"rmid");

    for (int s=1;s<(int)Ns+1;s++){
    
        Matrix_r(w,rmid,rupper,rlower);
            
        for (int i=0;i<Nr;i++){
            gamma=1.0-(ds/(pow(dr,2.0)))-((ds/2.0)*w[i]);
            betaL=ds/(2.0*pow(dr,2.0))-rterm(i);
            betaU=ds/(2.0*pow(dr,2.0))+rterm(i);
            //cout<<"i: "<<i<<" gamma: "<<gamma<<" betaL: "<<betaL<<" betaU: "<<betaU<<endl;
            
            if(i==0){
                bvecr[i]=gamma*q[0][s-1]+(betaL+betaU)*q[1][s-1];
            }
            else if(i==(int)Nr-1){
                bvecr[i]=gamma*q[Nr-1][s-1]+(betaL+betaU)*q[Nr-2][s-1];
            }
            else{
                bvecr[i]=gamma*q[i][s-1]+betaL*q[i-1][s-1]+betaU*q[i+1][s-1];
            }
        }
        
        //Apply TDMA
        TDMA(bvecr,Nr,rlower,rmid,rupper);
            
        //Now we have our solution for all i for s, from s-1. Full step completed.
        for (int i=0;i<Nr;i++){
            q[i][s]=bvecr[i];
            bvecr[i]=0.0;
            if (fabs(q[i][s])>1e5){
                cout<<i<<" "<<s<<" propagator problem: "<<q[i][s]<<endl;
                exit(EXIT_FAILURE);
            }
            
        }
        

    }

    destroy_1d_double_array(bvecr);
    destroy_1d_double_array(rupper);
    destroy_1d_double_array(rmid);
    destroy_1d_double_array(rlower);
    
}


