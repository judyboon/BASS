/*
    Bayesian group factor Analysis with Structured Sparsity (BASS)
    Copyright (C) 2015 Shiwen Zhao
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include"BASS.h"
#include"funWrapper.h"
#include"gig_par.h"

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>

#include <fstream>
#include <sstream>
#include <time.h> // used for setting seed
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

using namespace std;

// constructor
BASS::BASS(int * Parray, int v, int n, int k):
    v(v),n(n),k(k)
{
    P = gsl_vector_int_alloc(v);
    p = 0;
    for(int w=0; w<v; w++){
        gsl_vector_int_set(P,w,Parray[w]);
        p += Parray[w];
    }

    Y = gsl_matrix_alloc(p,n);

    Eta = gsl_matrix_alloc(k,n);
    Lambda = gsl_matrix_alloc(p,k);

    // sufficient statistics
    EtaTEta = gsl_matrix_alloc(k,k);
    EtaTY = gsl_matrix_alloc(k,p);

    // parameters
    Theta = gsl_matrix_alloc(p,k);
    Delta = gsl_matrix_alloc(p,k);
    Tau = gsl_matrix_alloc(v,k);
    Phi = gsl_matrix_alloc(v,k);
    Z = gsl_matrix_alloc(v,k);

    Gamma = gsl_vector_alloc(v);
    H = gsl_vector_alloc(v);
    Pi = gsl_vector_alloc(v);

    Sigma2inv = gsl_vector_alloc(p);

    a = b = 0.5;
    c = d = 0.5;
    e = f = 0.5;
    nu = 1.0;
    asig = 1.0;
    bsig = 0.25;

    // -----------------------------  random generator
    r=gsl_rng_alloc(gsl_rng_taus);

    // -----------------------------  for rgig
    mygig = new gig();
    myrs = new RngStream();

    // ----------------------------  for EM
    SigmaEta = gsl_matrix_alloc(k,k);
}

// destructor
BASS::~BASS()
{
    gsl_vector_int_free(P);
    gsl_matrix_free(Y);
    gsl_matrix_free(Eta);
    gsl_matrix_free(Lambda);

    gsl_matrix_free(EtaTEta);
    gsl_matrix_free(EtaTY);

    gsl_vector_free(Sigma2inv);

    gsl_matrix_free(Delta);
    gsl_matrix_free(Theta);
    gsl_matrix_free(Tau);
    gsl_matrix_free(Phi);

    gsl_matrix_free(Z);

    gsl_vector_free(Gamma);
    gsl_vector_free(H);
    gsl_vector_free(Pi);

    gsl_rng_free(r);
    delete mygig;
    delete myrs;

    gsl_matrix_free(SigmaEta);
}

void BASS::LoadData(string fileY, string sep){

    int p_ = 0, n_ = 0;
    fstream fs;
    string line, field;

    // -------------------------------------
    fs.open(fileY.c_str());
    if (! fs.is_open())
    {
        cout<<"Failed to open joint input matrix file"<<endl;
        exit(EXIT_FAILURE);
    }
    getline(fs,line);
    p_++;
    istringstream iss(line);
    if(sep.compare("space")==0){
        while(getline(iss,field,' ')){n_++;}
    }else if(sep.compare("tab")==0){
        while(getline(iss,field,'\t')){n_++;}
    }else{
        cout << "Please specify a valid separator." << endl << endl;
    }
    while(getline(fs,line)){p_++;}
    cout<<"The dimension of the joint input matrix is: "<< p_ << " " <<n_ <<endl;
    fs.close();

    if(p_ != p || n_ != n){
        cout << "Input joint matrix dimension does not match"<<endl;
        exit(EXIT_FAILURE);
    }


    // -------------------------------------- Actual read matrix
    FILE *fp = fopen(fileY.c_str(),"r");
    gsl_matrix_fscanf(fp,Y);
    fclose(fp);

    // --------------------------------------- Standardize data
    double rmean = 0.0;
    int j = 0;
    gsl_vector_view Vview;
    for(j = 0; j<p; j++){
        Vview = gsl_matrix_row(Y,j);
        rmean = gsl_stats_mean(Vview.vector.data,1,n);
        gsl_vector_add_constant(&Vview.vector,-1.0*rmean);
        double stdj = gsl_stats_sd(Vview.vector.data,1,n);
        gsl_vector_scale(&Vview.vector,1.0/stdj);
    }

}

void BASS::InitializePara()
{

    seed = time (NULL) * getpid();
    //seed = 8851296622221;
    gsl_rng_set (r, seed);    // set seed

    unsigned long seeds6 [6];
    for(int i = 0; i<6; i++){
        seeds6[i] = seed % (gsl_rng_uniform_int(r,500)+1);
    }
    myrs->SetSeed(seeds6);
    mygig->setRngStream(myrs);

    double gammaR,betaR;
    int pw;   // pw means the dimension of current data block
    int pcul = 0; // pcul means the sum of dimensions in previous data block
    for(int w=0; w<v; w++){
        gammaR = gsl_ran_gamma(r,f,1.0/nu);
        gsl_vector_set(Gamma,w,gammaR);
        gammaR = gsl_ran_gamma(r,e,1.0/gammaR);
        gsl_vector_set(H,w,gammaR);

        randomGamma(Tau->data+w*k,k,&d,1,&gammaR,1,r);
        randomGamma(Phi->data+w*k,k,&c,1,Tau->data+w*k,k,r);

        pw = gsl_vector_int_get(P,w);
        randomGamma(Delta->data+pcul*k,pw*k,&b,1,Phi->data+w*k,k,r);
        randomGamma(Theta->data+pcul*k,pw*k,&a,1,Delta->data+pcul*k,pw*k,r);
        pcul += pw;

        betaR = gsl_ran_beta(r,1,1);
        gsl_vector_set(Pi,w,betaR);
        randomBernoulli(Z->data+w*k,k,&betaR,1,r);
    }
    gsl_matrix_set_all(Z,0.0);

    double zero = 0.0;
    double one = 1.0;
    randomUniNormal(Lambda->data,p*k,&zero,1,&one,1,r);
    randomUniNormal(Eta->data,k*n,&zero,1,&one,1,r);
    randomGamma(Sigma2inv->data,p,&asig,1,&bsig,1,r);

    gsl_blas_dsyrk(CblasLower,CblasNoTrans,1.0,Eta,0.0,EtaTEta); // EtaTEta is lower triangular
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,Eta,Y,0.0,EtaTY);
}


void BASS::Write2File(string outputDir){

    string filename;
    FILE *fp;
    filename = outputDir + "/ETA";
    fp = fopen(filename.c_str(),"a"); gsl_matrix_fprintf(fp,Eta,"%e"); fclose(fp);

    filename = outputDir + "/LAMBDA";
    fp = fopen(filename.c_str(),"a"); gsl_matrix_fprintf(fp,Lambda,"%e"); fclose(fp);

    filename = outputDir + "/ZS";
    fp = fopen(filename.c_str(),"a"); gsl_matrix_fprintf(fp,Z,"%e"); fclose(fp);

    filename = outputDir + "/SIGMA2INV";
    fp = fopen(filename.c_str(),"a"); gsl_vector_fprintf(fp,Sigma2inv,"%e"); fclose(fp);

    filename = outputDir + "/THETA";
    fp = fopen(filename.c_str(),"a"); gsl_matrix_fprintf(fp,Theta,"%e"); fclose(fp);
    filename = outputDir + "/DELTA";
    fp = fopen(filename.c_str(),"a"); gsl_matrix_fprintf(fp,Delta,"%e"); fclose(fp);
    filename = outputDir + "/PHI";
    fp = fopen(filename.c_str(),"a"); gsl_matrix_fprintf(fp,Phi,"%e"); fclose(fp);
    filename = outputDir + "/TAU";
    fp = fopen(filename.c_str(),"a"); gsl_matrix_fprintf(fp,Tau,"%e"); fclose(fp);

    filename = outputDir + "/PI";
    fp = fopen(filename.c_str(),"a"); gsl_vector_fprintf(fp,Pi,"%e"); fclose(fp);
    filename = outputDir + "/H";
    fp = fopen(filename.c_str(),"a"); gsl_vector_fprintf(fp,H,"%e"); fclose(fp);
    filename = outputDir + "/GAMMA";
    fp = fopen(filename.c_str(),"a"); gsl_vector_fprintf(fp,Gamma,"%e"); fclose(fp);

}



// ---------------------------------------------- MCMC steps


void BASS::LambdaUpdateMCMC()
{
    using namespace MCMCFuns;

    // calculate sufficient statistics
    gsl_blas_dsyrk(CblasLower,CblasNoTrans,1.0,Eta,0.0,EtaTEta); // EtaTEta is lower triangular
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,Eta,Y,0.0,EtaTY);

    int pw;
    int pcul = 0;
    for(int w = 0; w < v; w++){
        pw = gsl_vector_int_get(P,w);
        gsl_matrix_view Lambdaw = gsl_matrix_submatrix(Lambda,pcul,0,pw,k);
        gsl_matrix_view Thetaw = gsl_matrix_submatrix(Theta,pcul,0,pw,k);
        gsl_vector_view Phiw = gsl_matrix_row(Phi,w);
        gsl_vector_view Zw = gsl_matrix_row(Z,w);
        gsl_vector_view Sigma2invw = gsl_vector_subvector(Sigma2inv,pcul,pw);
        gsl_matrix_view EtaTYw = gsl_matrix_submatrix(EtaTY,0,pcul,k,pw);

        updateLambda(&Lambdaw.matrix,EtaTEta,&EtaTYw.matrix,
                     &Thetaw.matrix,&Phiw.vector,&Zw.vector,&Sigma2invw.vector,r);

        pcul += pw;
    }

}


void BASS::EtaUpdateMCMC()
{
    using namespace MCMCFuns;

    updateEta(Eta,SigmaEta,Lambda,Y,Sigma2inv,r);
}

void BASS::HyperParamsUpdateMCMC(int collapse)
{
    using namespace MCMCFuns;

    int pw;
    int pcul = 0;
    for(int w = 0; w<v; w++){
        pw = gsl_vector_int_get(P,w);

        gsl_matrix_view Lambdaw = gsl_matrix_submatrix(Lambda,pcul,0,pw,k);
        gsl_matrix_view Thetaw = gsl_matrix_submatrix(Theta,pcul,0,pw,k);
        gsl_matrix_view Deltaw = gsl_matrix_submatrix(Delta,pcul,0,pw,k);
        gsl_vector_view Phiw = gsl_matrix_row(Phi,w);
        gsl_vector_view Tauw = gsl_matrix_row(Tau,w);
        gsl_vector_view Zw = gsl_matrix_row(Z,w);

        updateThetaDeltaPhi(&Thetaw.matrix,&Deltaw.matrix,&Phiw.vector,
                            &Lambdaw.matrix,&Tauw.vector,&Zw.vector,a,b,c,r,mygig);
        updateTauHGamma(&Tauw.vector,H->data+w,Gamma->data+w,&Phiw.vector,c,d,e,f,nu,r);


        if(collapse == 1)
            updatePiZCollapse(Pi->data+w,&Zw.vector,
                              &Lambdaw.matrix,&Thetaw.matrix,
                              &Deltaw.matrix,&Phiw.vector,a,b,r);
        else
            updatePiZnonCollapse(Pi->data+w,&Zw.vector,
                                 &Lambdaw.matrix,&Thetaw.matrix,
                                 &Deltaw.matrix,&Phiw.vector,a,b,r);

        pcul += pw;
    }

}

void BASS::Sigma2invUpdateMCMC()
{
    using namespace MCMCFuns;

    updateSigma2inv(Sigma2inv,Lambda,Eta,Y,asig,bsig,r);
}


// ============================================================ EM steps


void BASS::LambdaUpdateEM()
{

    using namespace EMFuns;

    // calculate sufficient statistics
    gsl_blas_dsyrk(CblasLower,CblasNoTrans,1.0,Eta,0.0,EtaTEta); // EtaTEta is lower triangular
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,Eta,Y,0.0,EtaTY);

    // Use SigmaEta to store EEtaTEta = EtaTEta + n*SigmaEta;
    // cholesky decomposition of EEtaTEta in updateEtaPX needs the lower part
    gsl_blas_dsyrk(CblasLower,CblasNoTrans,1.0,Eta,double(n),SigmaEta);

    updateLambda(Lambda,P,Theta,Phi,Z,Eta,Y,SigmaEta,Sigma2inv);
    //gsl_blas_dsyrk(CblasLower,CblasNoTrans,1.0,Eta,double(n),SigmaEta);
}

void BASS::EtaUpdateEM(int PX)
{
    using namespace EMFuns;

    if(PX == 1){
        updateEtaPX(Eta,SigmaEta,Lambda,Y,Sigma2inv);
    }
    else{
        updateEta(Eta,SigmaEta,Lambda,Y,Sigma2inv);
    }
}

void BASS::HyperParamsUpdateEM(int upz)
{
    using namespace EMFuns;

    int pw;
    int pcul = 0;
    for(int w = 0; w<v; w++){
        pw = gsl_vector_int_get(P,w);

        gsl_matrix_view Lambdaw = gsl_matrix_submatrix(Lambda,pcul,0,pw,k);
        gsl_matrix_view Thetaw = gsl_matrix_submatrix(Theta,pcul,0,pw,k);
        gsl_matrix_view Deltaw = gsl_matrix_submatrix(Delta,pcul,0,pw,k);
        gsl_vector_view Phiw = gsl_matrix_row(Phi,w);
        gsl_vector_view Tauw = gsl_matrix_row(Tau,w);
        gsl_vector_view Zw = gsl_matrix_row(Z,w);

        updateThetaDeltaPhi(&Thetaw.matrix,&Deltaw.matrix,&Phiw.vector,
                            &Lambdaw.matrix,&Tauw.vector,&Zw.vector,a,b,c);
        updateTauHGamma(&Tauw.vector,H->data+w,Gamma->data+w,&Phiw.vector,c,d,e,f,nu);

        if (upz == 0){
            gsl_vector_set_all(&Zw.vector,1e-1);
        }
        else{
            updatePiZnonCollapse(Pi->data+w,&Zw.vector,
                                 &Lambdaw.matrix,&Thetaw.matrix,
                                 &Deltaw.matrix,&Phiw.vector,a,b);
        }

        pcul += pw;
    }
}

void BASS::Sigima2invUpdateEM()
{
    using namespace EMFuns;

    updateSigma2inv(Sigma2inv,Lambda, Eta,Y,asig,bsig);

}


// ======================================================== Other functions


void BASS::DeleteEmptyFactor()
{

    deleteNullFactors(Lambda,Eta,SigmaEta,EtaTEta,EtaTY,
                      Theta,Delta,Phi,Tau,Z);
    k = Lambda->size2;
}

void BASS::PruneFactors()
{
    pruneFactors(Lambda,Eta,SigmaEta,EtaTEta,EtaTY,
                 Theta,Delta,Phi,Tau,Z,Sigma2inv);
    k = Lambda->size2;
}

void BASS::Rotate()
{
    // because updateLambda is called before updateEta,
    // so recalculate SigmaEta = 1/n*EtaTEta + SigmaEta;

    //gsl_matrix_set_identity(SigmaEta);
    gsl_blas_dsyrk(CblasLower,CblasNoTrans,1.0/n,Eta,1.0,SigmaEta);
    gsl_linalg_cholesky_decomp(SigmaEta);
    gsl_blas_dtrmm(CblasRight,CblasUpper,CblasTrans,CblasNonUnit,1.0, SigmaEta,Lambda);



}


void BASS::WriteRotation(string outputDir)
{
    // because updateLambda is called before updateEta,
    // so recalculate SigmaEta = 1/n*EtaTEta + SigmaEta;

    gsl_blas_dsyrk(CblasUpper,CblasNoTrans,1.0/n,Eta,1.0,SigmaEta);
    string filename;
    FILE *fp;
    filename = outputDir + "/A";
    fp = fopen(filename.c_str(),"a"); gsl_matrix_fprintf(fp,Eta,"%e"); fclose(fp);


}

int BASS::GetN0num()
{
    int num = countNZeroLoadingElement(Lambda);
    return num;
}

double BASS::GetFNorm()
{
    double fnorm = calculateMatrixNorm(Lambda);
    return fnorm;
}

double BASS::GetLogLikelihood()
{
    double logl = calculateLogLikelihood(Y,Lambda,Eta,Sigma2inv);
    return logl;
}


// ---------------------------------------------- VI steps
