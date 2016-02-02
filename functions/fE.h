/**************This is my function for calculating free energy from SCFT******/

double fE(Matrix newW, Matrix phi, Matrix chiMatrix, double volume){
    
    double fEW,fEchi;
    double fE_int;
    fEW=0.0;
    fEchi=0.0;
    fE_int=0.0;
    
    //inner end
    for(int ii=0;ii<ChainType;ii++){
        for(int jj=0;jj<ChainType;jj++){
            fEchi+=0.5*phi[ii][0]*chiMatrix[ii][jj]*phi[jj][0]*dV(0);
        }
        fEW+=0.5*(newW[ii][0]*phi[ii][0]*dV(0));
    }
    
    //outer end
    for(int ii=0;ii<ChainType;ii++){
        for(int jj=0;jj<ChainType;jj++){
            fEchi+=0.5*phi[ii][(int)Nr-1]*chiMatrix[ii][jj]*phi[jj][(int)Nr-1]*dV((int)Nr-1);
        }
        fEW+=0.5*(newW[ii][(int)Nr-1]*phi[ii][(int)Nr-1]*dV((int)Nr-1));
    }

    
    //middle
    for (int i=1;i<(int)Nr-1;i++){
        for(int ii=0;ii<ChainType;ii++){
            for(int jj=0;jj<ChainType;jj++){
                fEchi+=phi[ii][i]*chiMatrix[ii][jj]*phi[jj][i]*dV(i);
            }
            fEW+=(newW[ii][i]*phi[ii][i]*dV(i));
        }
    }

    //normalize by box size
    fEchi/=2.0*volume;
    fEW/=volume;
    
    fE_int=fEchi-fEW;

    return fE_int;
    
}