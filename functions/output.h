/****************Here I output the concentration profile for plotting***********************/
void output(double **phi){
    
    int i;
    ofstream outputFile1("./results/phi.dat");
    
    for(i=0;i<Nr;i++){
        outputFile1 <<i*dr<<" "<<phi[0][i]<<" "<<phi[1][i]<<" "<<phi[2][i]<<" "<<phi[3][i]<<" "<<phi[4][i]<<" "<<phi[5][i]<<std::endl;
    }
    
    outputFile1.close();
    
}

/*******************Here I output concentration profile for various radii***********************/
void outputphi(double **phi){
    
    ofstream outphi;
    string filename;
    
    if (Coord==2){
        if (poly==0){
            filename="./results/phi/ABAcyl"+DoubleToStr(r_0)+".dat";
            outphi.open(filename.c_str());
        }
        else if (poly==1){
            filename="./results/phi/ABcyl"+DoubleToStr(r_0)+".dat";
            outphi.open(filename.c_str());
        }
    }
    else if (Coord==3){
        if (poly==0){
            filename="./results/phi/ABAsph"+DoubleToStr(r_0)+".dat";
            outphi.open(filename.c_str());
        }
        else if (poly==1){
            filename="./results/phi/ABsph"+DoubleToStr(r_0)+".dat";
            outphi.open(filename.c_str());
        }
    }
    
    for (int i=0;i<Nr;i++){
        outphi<<i*dr<<" "<<phi[0][i]<<" "<<phi[1][i]<<" "<<phi[2][i]<<" "<<phi[3][i]<<" "<<phi[4][i]<<" "<<phi[5][i]<<std::endl;
    }
    
    outphi.close();
    
}

/*************************Here I output quadratic and quartic fit parameters****************************************/
void outputkappa(double *a1,double *a2,double *a3,double *a4,double *a5,double *a6, int nfa,double *chi, int rad_type){
    
    ofstream outkappa;
    string filename;
    
        if (Coord==2){
            if (poly==0){
                filename="./results/ABAcyl"+DoubleToStr(chi[0])+".dat";
                outkappa.open(filename.c_str());
            }
            else if (poly==1){
                filename="./results/ABcyl"+DoubleToStr(chi[0])+".dat";
                outkappa.open(filename.c_str());
            }
        }
        else if (Coord==3){
            if (poly==0){
                filename="./results/ABAsph"+DoubleToStr(chi[0])+".dat";
                outkappa.open(filename.c_str());
            }
            else if (poly==1){
                filename="./results/ABsph"+DoubleToStr(chi[0])+".dat";
                outkappa.open(filename.c_str());
            }
        }
    
    
    for (int i=0;i<nfa;i++){
        outkappa<<0.3+(double)i*0.02<<"  "<<a1[i]<<" "<<a2[i]<<" "<<a3[i]<<" "<<a4[i]<<" "<<a5[i]<<" "<<a6[i]<<endl;
    }
    
}

/***************Here I output a concentration profile**********************/
void outputphi_fa(double **phi, double fA, int nfa){
    
    ofstream outphifa;
    string filename;
    
    if (Coord==2){
        if (poly==0){
            filename="./results/phi_FA/ABAcyl"+DoubleToStr(fA)+".dat";
            outphifa.open(filename.c_str());
        }
        else if (poly==1){
            filename="./results/phi_FA/ABcyl"+DoubleToStr(fA)+".dat";
            outphifa.open(filename.c_str());
        }
    }
    else if (Coord==3){
        if (poly==0){
            filename="./results/phi_FA/ABAsph"+DoubleToStr(fA)+".dat";
            outphifa.open(filename.c_str());
        }
        else if (poly==1){
            filename="./results/phi_FA/ABsph"+DoubleToStr(fA)+".dat";
            outphifa.open(filename.c_str());
        }
    }
    else if (Coord==1){
        if (poly==0){
            filename="./results/phi_FA/ABAplane"+DoubleToStr(fA)+".dat";
            outphifa.open(filename.c_str());
        }
        else if (poly==1){
            filename="./results/phi_FA/ABplane"+DoubleToStr(fA)+".dat";
            outphifa.open(filename.c_str());
        }
    }
    
    for (int i=0;i<Nr;i++){
        outphifa<<(double)i*dr<<"  "<<phi[0][i]<<" "<<phi[1][i]<<" "<<phi[2][i]<<" "<<phi[3][i]<<" "<<phi[4][i]<<" "<<phi[5][i]<<std::endl;
    }
    
    
}

void outputfE_FA(double fA,double *Curv,double *dFE,int nradii){
    
    fstream outfE_FA;
    string filename;
    
    if (Coord==2){
        if (poly==0){
            filename="./results/fE/ABAcyl"+DoubleToStr(fA)+".dat";
            outfE_FA.open(filename.c_str());
        }
        else if (poly==1){
            filename="./results/fE/ABcyl"+DoubleToStr(fA)+".dat";
            outfE_FA.open(filename.c_str());
        }
    }
    else if (Coord==3){
        if (poly==0){
            filename="./results/fE/ABAsph"+DoubleToStr(fA)+".dat";
            outfE_FA.open(filename.c_str());
        }
        else if (poly==1){
            filename="./results/fE/ABsph"+DoubleToStr(fA)+".dat";
            outfE_FA.open(filename.c_str());
        }
    }
    
    for (int i=0;i<nradii;i++){
        outfE_FA<<fA<<Curv[i]<<dFE[i]<<endl;
    }
    
    
    
    
}



