#ifndef BASS_H
#define BASS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <iostream>

#include <gig_par.h>

class BASS
{
public:
    // members elements
    gsl_vector_int * P;
    int p;
    int v;  // the number of matrices
    int k;
    int n;

    // input data
    gsl_matrix * Y;

    // model parameters
    gsl_matrix * Lambda;  // when using PXL-EM, this is Lambda*
    gsl_matrix * Eta;
    gsl_matrix * SigmaEta;
    gsl_vector * Sigma2inv;

    // expected sufficient statistics
    gsl_matrix * EtaTEta;  // k * k, is a lower triangular matrix
    gsl_matrix * EtaTY;  // k * p

    // TPB prior, hyperparameters
    gsl_matrix * Theta;
    gsl_matrix * Delta;
    gsl_matrix * Phi;  // Phi here is a matrix, v * k
    gsl_matrix * Tau;  // v * k
    gsl_matrix * Z;  // v * k

    gsl_vector * Gamma; // 1 * v
    gsl_vector * H; // 1 * v
    gsl_vector * Pi;  // 1 * v

    double a,b,c,d,e,f,nu;
    double asig,bsig;

    // random seed
    unsigned long seed;
    gsl_rng * r;
    gig * mygig; // to sample from Generalized Inverse Gaussian
    RngStream * myrs; // a random U[0,1] generator for GIG


public:
    // member methods
    BASS(int *Parray, int v, int n, int k);
    ~BASS();

    void LoadData(std::string fileY, std::string sep);
    void InitializePara();

    // MCMC methods
    void LambdaUpdateMCMC();
    void EtaUpdateMCMC();
    void HyperParamsUpdateMCMC(int collapse);
    void Sigma2invUpdateMCMC();

    // EM methods
    void LambdaUpdateEM();
    void EtaUpdateEM(int);
    void HyperParamsUpdateEM(int);
    void Sigima2invUpdateEM();

    // VI methods
    // void BlockAUpdateVI();
    // void BlockB1UpdateVI();
    // void BlockB2UpdateVI();
    // void Sigma2invUpdateVI();

public:
    int GetN0num();
    double GetFNorm();
    double GetLogLikelihood();

    void DeleteEmptyFactor();
    void PruneFactors();
    void Rotate();
    // due to identifiability problem, specific latent factors may occur in shared factor
    // PruneFactors() delete those cases;
    void WriteRotation(std::string outputDir);
    void Write2File(std::string outputDir);

};

#endif // BASS_H
