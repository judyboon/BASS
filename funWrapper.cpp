#include "funWrapper.h"

#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>

using namespace std;

/* ------------------------------------------------------------------- */
int randomGamma(double * rv, int lrv, const double * a, int la,
                const double * b, int lb, gsl_rng * r)
{
    gsl_vector * A = gsl_vector_alloc(lrv);
    gsl_vector * B = gsl_vector_alloc(lrv);
    /* deal with the parameter a */
    if(la == 1){
        gsl_vector_set_all(A, *a);
    }
    else{
        int dev = lrv/la;
        for(int i = 0; i<dev; i++){
            for(int j = 0; j<la; j++){
                gsl_vector_set(A,i*la+j,*(a+j));
            }
        }
        for(int j = 0; j<lrv%la; j++){
            gsl_vector_set(A,dev*la+j,*(a+j));
        }
    }
    /* deal with parameter b */
    if(lb == 1){
        gsl_vector_set_all(B, *b);
    }
    else{
        int dev = lrv/lb;
        for(int i = 0; i<dev; i++){
            for(int j = 0; j<lb; j++){
                gsl_vector_set(B,i*lb+j,*(b+j));
            }
        }
        for(int j = 0; j<lrv%lb; j++){
            gsl_vector_set(B,dev*lb+j,*(b+j));
        }
    }
    /* sample random variable */
    for(int i = 0; i<lrv; i++){
        *(rv + i) = gsl_ran_gamma(r,gsl_vector_get(A,i),
                                  1.0/gsl_vector_get(B,i));
    }
    gsl_vector_free(A);
    gsl_vector_free(B);

    return 1;
}

int randomUniNormal(double * rv, int lrv, const double * mu, int lm,
                    const double * sigma, int ls, gsl_rng *r){
    gsl_vector * Mu = gsl_vector_alloc(lrv);
    gsl_vector * Sigma = gsl_vector_alloc(lrv);
    /* deal with the parameter mu */
    if(lm == 1){
        gsl_vector_set_all(Mu, *mu);
    }
    else{
        int dev = lrv/lm;
        for(int i = 0; i<dev; i++){
            for(int j = 0; j<lm; j++){
                gsl_vector_set(Mu,i*lm+j,*(mu+j));
            }
        }
        for(int j = 0; j<lrv%lm; j++){
            gsl_vector_set(Mu,dev*lm+j,*(mu+j));
        }
    }
    /* deal with parameter sigma */
    if(ls == 1){
        gsl_vector_set_all(Sigma, *sigma);
    }
    else{
        int dev = lrv/ls;
        for(int i = 0; i<dev; i++){
            for(int j = 0; j<ls; j++){
                gsl_vector_set(Sigma,i*ls+j,*(sigma+j));
            }
        }
        for(int j = 0; j<lrv%ls; j++){
            gsl_vector_set(Sigma,dev*ls+j,*(sigma+j));
        }
    }
    /* sample random variable */
    for(int i = 0; i<lrv; i++){
        *(rv + i) = gsl_vector_get(Mu,i) +
                gsl_ran_gaussian(r,gsl_vector_get(Sigma,i));
    }
    gsl_vector_free(Mu);
    gsl_vector_free(Sigma);

    return 1;
}

int randomGIG(double * rv, int lrv, const double * p, int lp,
              const double * a, int la, const double * b, int lb, gig * mygig)
{
    gsl_vector * P = gsl_vector_alloc(lrv);
    gsl_vector * A = gsl_vector_alloc(lrv);
    gsl_vector * B = gsl_vector_alloc(lrv);
    /* deal with the parameter p */
    if(lp == 1){
        gsl_vector_set_all(P, *p);
    }
    else{
        int dev = lrv/lp;
        for(int i = 0; i<dev; i++){
            for(int j = 0; j<lp; j++){
                gsl_vector_set(P,i*lp+j,*(p+j));
            }
        }
        for(int j = 0; j<lrv%lp; j++){
            gsl_vector_set(P,dev*lp+j,*(p+j));
        }
    }
    /* deal with parameter a */
    if(la == 1){
        gsl_vector_set_all(A, *a);
    }
    else{
        int dev = lrv/la;
        for(int i = 0; i<dev; i++){
            for(int j = 0; j<la; j++){
                gsl_vector_set(A,i*la+j,*(a+j));
            }
        }
        for(int j = 0; j<lrv%la; j++){
            gsl_vector_set(A,dev*la+j,*(a+j));
        }
    }
    /* deal with parameter b */
    if(lb == 1){
        gsl_vector_set_all(B, *b);
    }
    else{
        int dev = lrv/lb;
        for(int i = 0; i<dev; i++){
            for(int j = 0; j<lb; j++){
                gsl_vector_set(B,i*lb+j,*(b+j));
            }
        }
        for(int j = 0; j<lrv%lb; j++){
            gsl_vector_set(B,dev*lb+j,*(b+j));
        }
    }

    for(int i = 0; i<lrv; i++){
        mygig->rgig(*(rv + i), &(P->data[i]), &(A->data[i]), &(B->data[i]));
        //cout<<*(rv + i)<<" ";
    }
    //cout<<endl;
    gsl_vector_free(P);
    gsl_vector_free(A);
    gsl_vector_free(B);

    return 1;
}

int randomBernoulli(double * rv, int lrv, const double * p, int lp, gsl_rng * r)
{
    gsl_vector * P = gsl_vector_alloc(lrv);
    if(lp == 1){
        gsl_vector_set_all(P,*p);
    }
    else{
        int dev = lrv/lp;
        for(int i = 0; i<dev; i++){
            for(int j = 0; j<lp; j++){
                gsl_vector_set(P,i*lp+j,*(p+j));
            }
        }
        for(int j = 0; j<lrv%lp; j++){
            gsl_vector_set(P,dev*lp+j,*(p+j));
        }
    }
    for(int i = 0; i<lrv; i++){
        *(rv + i) = double(gsl_ran_bernoulli(r,gsl_vector_get(P,i)));
    }
    gsl_vector_free(P);

    return 1;
}

/* ------------------------------------------------------------------- */
int mulM1M2(gsl_matrix * Dest, const gsl_matrix * M1, const gsl_matrix * M2)
{
    if(M1->size2 != M2->size1 || Dest->size1 != M1->size1 ||
            Dest->size2 != M2->size2){
        cout<<"matrix matrix multiplication dimension error"<<endl;
        exit(EXIT_FAILURE);
    }
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,M1,M2,0.0,Dest);

    /*
    int m = M1->size1;
    int k = M1->size2;
    int n = M2->size2;
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                m,n,k,1.0,
                M1->data,k,M2->data,n,0.0,Dest->data,n);
    /*
    double v = 0.0;
    for(int i = 0; i < Dest->size1; i ++){
        for(int j = 0; j < Dest->size2; j ++){
            v = 0.0;
            for(int k = 0; k<M1->size2; k++){
                v += gsl_matrix_get(M1,i,k) * gsl_matrix_get(M2,k,j);
            }
            gsl_matrix_set(Dest,i,j,v);
        }
    }*/

    return 1;
}


int mulM1V1byrow(gsl_matrix * M1, const gsl_vector * V1)
{
    if(M1->size2 != V1->size){
        cout<<"matrix vector by row multiplication dimension error"<<endl;
        exit(EXIT_FAILURE);
    }
    gsl_vector_view Vview;
    for(int i = 0; i<M1->size1; i++){
        Vview = gsl_matrix_row(M1,i);
        gsl_vector_mul(&Vview.vector,V1);
    }

    return 1;
}

int divM1V1byrow(gsl_matrix * M1, const gsl_vector * V1)
{
    if(M1->size2 != V1->size){
        cout<<"matrix vector by row division dimension error"<<endl;
        exit(EXIT_FAILURE);
    }
    gsl_vector_view Vview;
    for(int i = 0; i<M1->size1; i++){
        Vview = gsl_matrix_row(M1,i);
        gsl_vector_div(&Vview.vector,V1);
    }

    return 1;
}

int mulM1V1bycol(gsl_matrix * M1, const gsl_vector * V1)
{
    if(M1->size1 != V1->size){
        cout<<"matrix vector by column multiplication dimension error"<<endl;
        exit(EXIT_FAILURE);
    }
    gsl_vector_view Vview;
    for(int j = 0; j<M1->size2; j++){
        Vview = gsl_matrix_column(M1,j);
        gsl_vector_mul(&Vview.vector,V1);
    }

    return 1;
}

int divM1V1bycol(gsl_matrix * M1, const gsl_vector * V1)
{
    if(M1->size1 != V1->size){
        cout<<"matrix vector by column division dimension error"<<endl;
        exit(EXIT_FAILURE);
    }
    gsl_vector_view Vview;
    for(int j = 0; j<M1->size2; j++){
        Vview = gsl_matrix_column(M1,j);
        gsl_vector_div(&Vview.vector,V1);
    }

    return 1;
}

int mulM1diagV1M2(gsl_matrix * Dest, const gsl_matrix * M1,
                  const gsl_vector * V1, const gsl_matrix * M2)
{
    if(M1->size2 != V1->size || M2->size1 != V1->size){
        cout<<"matrix diagonal vector matrix multiplication error"<<endl;
        exit(EXIT_FAILURE);
    }
    int n1 = M1->size1;
    int n2 = M2->size2;
    int p = M1->size2;

    double v = 0.0;
    for(int i = 0; i<n1; i++){
        for(int j = 0; j<n2; j++){
            v = 0.0;
            for(int t = 0; t<p; t++){
                v += gsl_matrix_get(M1,i,t)*gsl_vector_get(V1,t)*gsl_matrix_get(M2,t,j);
            }
            gsl_matrix_set(Dest,i,j,v);
        }
    }

    return 1;
}

int addM1V1diag(gsl_matrix * M1, const gsl_vector * V1)
{
    if(M1->size1 != V1->size && M1->size2 != V1->size){
        cout<<"matrix vector add diagonal dimension error"<<endl;
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i< V1->size; i++){
        gsl_matrix_set(M1,i,i,
                       gsl_matrix_get(M1,i,i) + gsl_vector_get(V1,i));
    }
    return 1;
}


/* ------------------------------------------------------------------ */
/* -----------------  namespace for MCMC functions ------------------ */
/* ------------------------------------------------------------------ */
namespace MCMCFuns{

void updateLambda(gsl_matrix * Lambda,
                  const gsl_matrix * EtaTEta, const gsl_matrix * EtaTY,
                  const gsl_matrix * Theta, const gsl_vector * Phi,
                  const gsl_vector * Z, const gsl_vector * Sigma2inv,
                  gsl_rng * r)
{
    int p = Lambda->size1;
    int k = Lambda->size2;
    double Djh;
    gsl_vector_view Vview;

    double zero = 0.0;
    double one = 1.0;
    gsl_matrix * SigmaLambdaj = gsl_matrix_alloc(k,k);
    gsl_vector * MuLambdaj = gsl_vector_alloc(k);

    for(int j = 0; j<p; j++){
        gsl_matrix_memcpy(SigmaLambdaj,EtaTEta);
        gsl_matrix_scale(SigmaLambdaj,gsl_vector_get(Sigma2inv,j));

        // Calculate SigmaLambdaj
        for(int h = 0; h<k; h++){
            if(gsl_vector_get(Z,h) == 1.0)
                Djh = 1.0/gsl_matrix_get(Theta,j,h);
            else
                Djh = 1.0/gsl_vector_get(Phi,h);
            gsl_matrix_set(SigmaLambdaj,h,h,
                           gsl_matrix_get(SigmaLambdaj,h,h) + Djh);
        }

        // Cholesky decomposition of SigmaLambdaj to calculate MulLambdaj
        gsl_linalg_cholesky_decomp(SigmaLambdaj);
        gsl_vector_const_view cVview = gsl_matrix_const_column(EtaTY,j);
        gsl_linalg_cholesky_solve(SigmaLambdaj,&cVview.vector,MuLambdaj);
        gsl_vector_scale(MuLambdaj,gsl_vector_get(Sigma2inv,j));

        // Cholesky decomposition of SigmaLambdaj to sample from MVN
        // get one sample from this MVN
        gsl_linalg_cholesky_invert(SigmaLambdaj);
        gsl_linalg_cholesky_decomp(SigmaLambdaj);
        Vview = gsl_matrix_row(Lambda,j);
        randomUniNormal(Vview.vector.data,k,&zero,1,&one,1,r);
        gsl_blas_dtrmv(CblasUpper,CblasTrans,CblasNonUnit,SigmaLambdaj,&Vview.vector);
        //cblas_dtrmv(CblasRowMajor,CblasUpper,CblasTrans,CblasNonUnit,
        //            k,SigmaLambdaj->data,k,Vview.vector.data,1);
        gsl_vector_add(&Vview.vector,MuLambdaj);
    }

    gsl_matrix_free(SigmaLambdaj);
    gsl_vector_free(MuLambdaj);
}


void updateEta(gsl_matrix * Eta, gsl_matrix * SigmaEta,
               const gsl_matrix * Lambda, const gsl_matrix * Y,
               const gsl_vector * Sigma2inv, gsl_rng * r)
{

    int n = Y->size2;
    int k = Eta->size1;
    int p = Lambda->size1;

    gsl_matrix * SigmaEtainv = gsl_matrix_alloc(k,k);
    gsl_matrix * LambdaTSigmaY = gsl_matrix_alloc(k,n);
    gsl_vector_view Vview;

    /* --------- using blas to calculate SigmaEta and LambdaTSigmaY --------- */
    // calculate Lambda * Sigma2inv
    gsl_matrix * LambdaTSigma = gsl_matrix_alloc(k,p);
    gsl_matrix_transpose_memcpy(LambdaTSigma,Lambda);
    for(int h = 0; h<k; h++){
        Vview = gsl_matrix_row(LambdaTSigma,h);
        gsl_vector_mul(&Vview.vector,Sigma2inv);
    }
    // calculate SigmaEta
    gsl_matrix_set_identity(SigmaEta);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,LambdaTSigma,Lambda,1.0,SigmaEta);
    // calculate LambdaTSigmaY
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,LambdaTSigma,Y,0.0,LambdaTSigmaY);


    double one = 1.0;
    double zero = 0.0;
    gsl_linalg_cholesky_decomp(SigmaEta); // chol of SigmaEta
    gsl_matrix_memcpy(SigmaEtainv,SigmaEta);
    gsl_linalg_cholesky_invert(SigmaEta);
    gsl_linalg_cholesky_decomp(SigmaEta);

    gsl_vector * MuEta = gsl_vector_alloc(k);
    gsl_vector * MVNunit = gsl_vector_alloc(k);
    for(int i = 0; i<n; i++){
        // calculate MuEta
        Vview = gsl_matrix_column(LambdaTSigmaY,i);
        gsl_linalg_cholesky_solve(SigmaEtainv,&Vview.vector,MuEta);

        // get one sample from this MVN
        Vview = gsl_matrix_column(Eta,i);
        randomUniNormal(MVNunit->data,k,&zero,1,&one,1,r);
        gsl_blas_dtrmv(CblasUpper,CblasTrans,CblasNonUnit,SigmaEta,MVNunit);
        //cblas_dtrmv(CblasRowMajor,CblasUpper,CblasTrans,CblasNonUnit,
        //            k,SigmaEta->data,k,MVNunit->data,1);
        gsl_vector_add(MVNunit,MuEta);
        gsl_vector_memcpy(&Vview.vector,MVNunit);

        //gsl_blas_dsymv(CblasLower,one,SigmaEta,&Vview.vector,one,MuEta);
        //gsl_vector_memcpy(&Vview.vector,MuEta);
    }

    gsl_matrix_free(SigmaEtainv);
    gsl_matrix_free(LambdaTSigmaY);
    gsl_vector_free(MuEta);
    gsl_vector_free(MVNunit);

    gsl_matrix_free(LambdaTSigma);
}


void updateThetaDeltaPhi(gsl_matrix * Theta, gsl_matrix * Delta,
                         gsl_vector * Phi,
                         const gsl_matrix * Lambda, const gsl_vector * Tau,
                         const gsl_vector * Z, double a, double b, double c,
                         gsl_rng * r, gig * mygig)
{
    int p = Theta->size1;
    int k = Theta->size2;

    double rgigR,rgigP,rgigA,rgigB;
    double gammaA, gammaB, gammaR;
    double phih,thetajh;
    double sumDeltah;
    gsl_vector_view Vview;
    for(int h = 0; h<k; h++){
        if(gsl_vector_get(Z,h) == 1.0){
            phih = gsl_vector_get(Phi,h);
            sumDeltah = 0.0;
            for(int j = 0; j<p; j++){
                rgigP = a - 0.5;
                rgigA = gsl_matrix_get(Delta,j,h) * 2.0;
                rgigB = gsl_matrix_get(Lambda,j,h) * gsl_matrix_get(Lambda,j,h);
                mygig->rgig(rgigR,&rgigP,&rgigA,&rgigB);
                gsl_matrix_set(Theta,j,h,rgigR);

                gammaA = a+b;
                gammaB = rgigR + phih;
                gammaR = gsl_ran_gamma(r,gammaA,1.0/gammaB);
                gsl_matrix_set(Delta,j,h,gammaR);

                sumDeltah += gammaR;
            }
            gammaA = p*b + c;
            gammaB = sumDeltah + gsl_vector_get(Tau,h);
            gammaR = gsl_ran_gamma(r,gammaA,1.0/gammaB);
            gsl_vector_set(Phi,h,gammaR);

        }
        else{
            phih = gsl_vector_get(Phi,h);
            rgigP = c - p*0.5;
            rgigA = 2.0*gsl_vector_get(Tau,h);
            rgigB = 0.0; // sum Lambda2h

            for(int j = 0; j<p; j++){
                rgigB += gsl_matrix_get(Lambda,j,h) * gsl_matrix_get(Lambda,j,h);
                // for Deltajh
                thetajh = gsl_matrix_get(Theta,j,h);
                gammaA = a+b;
                gammaB = thetajh + phih;
                gammaR = gsl_ran_gamma(r,gammaA,1.0/gammaB);
                gsl_matrix_set(Delta,j,h,gammaR);
            }
            /* ---- */
            mygig->rgig(rgigR,&rgigP,&rgigA,&rgigB);
            gsl_vector_set(Phi,h,rgigR);

            Vview = gsl_matrix_column(Theta,h);
            gsl_vector_set_all(&Vview.vector,rgigR);
        }
    }

}


void updateTauHGamma(gsl_vector * Tau, double * H, double * Gamma,
                     const gsl_vector * Phi, double c, double d, double e,
                     double f,double nu, gsl_rng * r)
{
    int k = Tau->size;
    double gammaA,gammaB,gammaR;
    gammaA = c+d;
    double sumTau = 0.0;
    for(int h = 0; h<k; h++){
        gammaB = gsl_vector_get(Phi,h) + *H;
        gammaR = gsl_ran_gamma(r,gammaA,1.0/gammaB);
        gsl_vector_set(Tau,h,gammaR);
        sumTau += gammaR;
    }

    gammaA = k*d + e;
    gammaB = *Gamma + sumTau;
    *H = gsl_ran_gamma(r,gammaA,1.0/gammaB);

    gammaA = e + f;
    gammaB = *H + nu;
    *Gamma = gsl_ran_gamma(r,gammaA,1.0/gammaB);
}


void updatePiZCollapse(double * Pi, gsl_vector * Z,
                       const gsl_matrix * Lambda, gsl_matrix * Theta,
                       gsl_matrix * Delta, const gsl_vector * Phi,
                       double a, double b, gsl_rng * r)
{
    int k = Z->size;
    int p = Lambda->size1;
    int Z1n, Z0n;
    Z1n = 0;
    for(int h=0; h<k; h++){
        Z1n += int(gsl_vector_get(Z,h));
    }
    Z0n = k - Z1n;
    *Pi = gsl_ran_beta(r,1.0+Z1n,1.0+Z0n);

    double sigma0,sigma1,dnorm,sumCol0,sumCol1;
    double logPmax, PZ1;
    double nx;
    int Zh;
    Z1n = 0;
    double gammaR,gammaA,gammaB;
    gammaA = a+b;
    gsl_vector_view Vview;
    for(int h = 0; h<k; h++){
        sigma0 = gsl_vector_get(Phi,h);
        sumCol0 = 0;
        sumCol1 = 0;
        for(int j = 0; j<p; j++){
            nx = gsl_matrix_get(Lambda,j,h);
            dnorm = -0.5*gsl_sf_log(2*M_PI*sigma0) - nx*nx/(2.0*sigma0);
            sumCol0 += dnorm;

            sigma1 = gsl_matrix_get(Theta,j,h);
            dnorm = -0.5*gsl_sf_log(2*M_PI*sigma1) - nx*nx/(2.0*sigma1);
            sumCol1 += dnorm;

            sumCol1 += (a-1)*gsl_sf_log(sigma1);
            sumCol1 -= (a+b)*gsl_sf_log(sigma1 + sigma0);
        }
        sumCol0 += gsl_sf_log(1 - *Pi);
        sumCol1 += gsl_sf_log(*Pi);
        sumCol1 += b*p*gsl_sf_log(sigma0);
        sumCol1 += p*(gsl_sf_lngamma(a+b) - gsl_sf_lngamma(a) - gsl_sf_lngamma(b));

        logPmax = max(sumCol0,sumCol1);
        sumCol0 -= logPmax;
        sumCol1 -= logPmax;
        if(sumCol0 < -500.0){
            sumCol0 = 0.0;
            sumCol1 = 1.0;
        }
        else if(sumCol1 < -500.0){
            sumCol0 = 1.0;
            sumCol1 = 0.0;
        }
        else{
            sumCol0 = gsl_sf_exp(sumCol0);
            sumCol1 = gsl_sf_exp(sumCol1);
        }
        PZ1 = sumCol1/(sumCol1 + sumCol0);
        Zh = gsl_ran_bernoulli(r,PZ1);

        gsl_vector_set(Z,h,double(Zh));
        Z1n += Zh;

        if(Zh == 1){
            for(int j = 0; j<p; j++){
                gammaB = gsl_matrix_get(Theta,j,h) + gsl_vector_get(Phi,h);
                gammaR = gsl_ran_gamma(r,gammaA,1.0/gammaB);
                gsl_matrix_set(Delta,j,h,gammaR);
            }
        }
        else{
            Vview = gsl_matrix_column(Theta,h);
            gsl_vector_set_all(&Vview.vector,gsl_vector_get(Phi,h));

            for(int j = 0; j<p; j++){
                gammaB = gsl_matrix_get(Theta,j,h) + gsl_vector_get(Phi,h);
                gammaR = gsl_ran_gamma(r,gammaA,1.0/gammaB);
                gsl_matrix_set(Delta,j,h,gammaR);
            }

        }
    }
    Z0n = k - Z1n;
}



void updatePiZnonCollapse(double * Pi, gsl_vector * Z,
                          const gsl_matrix * Lambda, gsl_matrix * Theta,
                          gsl_matrix * Delta, const gsl_vector * Phi,
                          double a, double b, gsl_rng * r)
{
    int k = Z->size;
    int p = Lambda->size1;
    int Z1n, Z0n;
    Z1n = 0;
    for(int h=0; h<k; h++){
        Z1n += int(gsl_vector_get(Z,h));
    }
    Z0n = k - Z1n;
    *Pi = gsl_ran_beta(r,1.0+Z1n,1.0+Z0n);

    double sigma0,sigma1,dnorm,dgamma,sumCol0,sumCol1;
    double ga,gb,gx;
    double nx;
    double logPmax, PZ1;
    int Zh;
    Z1n = 0;
    double gammaR,gammaA,gammaB;
    gammaA = a+b;
    gsl_vector_view Vview;
    for(int h = 0; h<k; h++){
        sigma0 = gsl_vector_get(Phi,h);
        sumCol0 = 0;
        sumCol1 = 0;
        for(int j = 0; j<p; j++){

            nx = gsl_matrix_get(Lambda,j,h);
            dnorm = -0.5*gsl_sf_log(2*M_PI*sigma0) - nx*nx/(2.0*sigma0);
            sumCol0 += dnorm;

            sigma1 = gsl_matrix_get(Theta,j,h);
            dnorm = -0.5*gsl_sf_log(2*M_PI*sigma1) - nx*nx/(2.0*sigma1);
            sumCol1 += dnorm;

            ga = a;
            gb = gsl_matrix_get(Delta,j,h);
            gx = gsl_matrix_get(Theta,j,h);
            dgamma = ga*gsl_sf_log(gb) - gsl_sf_lngamma(ga) +
                    (ga-1)*gsl_sf_log(gx) - gx*gb;
            sumCol1 += dgamma;

            ga = b;
            gb = gsl_vector_get(Phi,h);
            gx = gsl_matrix_get(Delta,j,h);
            dgamma = ga*gsl_sf_log(gb) - gsl_sf_lngamma(ga) +
                    (ga-1)*gsl_sf_log(gx) - gx*gb;
            sumCol1 += dgamma;
        }
        sumCol0 += gsl_sf_log(1 - *Pi);
        sumCol1 += gsl_sf_log(*Pi);

        logPmax = max(sumCol0,sumCol1);
        sumCol0 -= logPmax;
        sumCol1 -= logPmax;
        if(sumCol0 < -500.0){
            sumCol0 = 0.0;
            sumCol1 = 1.0;
        }
        else if(sumCol1 < -500.0){
            sumCol0 = 1.0;
            sumCol1 = 0.0;
        }
        else{
            sumCol0 = gsl_sf_exp(sumCol0);
            sumCol1 = gsl_sf_exp(sumCol1);
        }
        PZ1 = sumCol1/(sumCol1 + sumCol0);
        Zh = gsl_ran_bernoulli(r,PZ1);

        gsl_vector_set(Z,h,double(Zh));
        Z1n += Zh;

        if(Zh == 1){
            for(int j = 0; j<p; j++){
                gammaB = gsl_matrix_get(Theta,j,h) + gsl_vector_get(Phi,h);
                gammaR = gsl_ran_gamma(r,gammaA,1.0/gammaB);
                gsl_matrix_set(Delta,j,h,gammaR);
            }
        }
        else{
            Vview = gsl_matrix_column(Theta,h);
            gsl_vector_set_all(&Vview.vector,gsl_vector_get(Phi,h));
            for(int j = 0; j<p; j++){
                gammaB = gsl_matrix_get(Theta,j,h) + gsl_vector_get(Phi,h);
                gammaR = gsl_ran_gamma(r,gammaA,1.0/gammaB);
                gsl_matrix_set(Delta,j,h,gammaR);
            }

        }
    }
    Z0n = k - Z1n;
}



void updateSigma2inv(gsl_vector * Sigma2inv,
                     const gsl_matrix * Lambda,
                     const gsl_matrix * Eta, const gsl_matrix * Y,
                     const double a, const double b,
                     gsl_rng * r)
{
    int p = Sigma2inv->size;
    int k = Lambda->size2;
    int n = Y->size2;

    gsl_matrix * Ymean = gsl_matrix_alloc(p,n);
    mulM1M2(Ymean,Lambda,Eta);
    gsl_matrix_sub(Ymean,Y); // Ymean becomes error

    double gammaA,gammaB,gammaR;
    gammaA = n*.5+a;
    gsl_vector_view Vview;
    for(int j = 0; j<p; j++){
        gammaB = 0.0;
        Vview = gsl_matrix_row(Ymean,j);
        gsl_blas_ddot(&Vview.vector,&Vview.vector,&gammaB);

        gammaB = gammaB*.5 + b;
        gammaR = gsl_ran_gamma(r,gammaA,1.0/gammaB);
        gsl_vector_set(Sigma2inv,j,gammaR);
    }
    gsl_matrix_free(Ymean);
}


}


/* ----------------------------------------------------------------- */
/* ---------------- namespace for EM functions --------------------- */
/* ----------------------------------------------------------------- */

namespace EMFuns{

void updateLambda(gsl_matrix * Lambda, const gsl_vector_int * P, const gsl_matrix * Theta,
                  const gsl_matrix * Phi, const gsl_matrix * EZ,
                  const gsl_matrix * EEta, const gsl_matrix * Y,
                  const gsl_matrix * EEtaTEta, const gsl_vector * Sigma2inv)
{
    int p = Y->size1;
    int k = EEta->size1;
    int n = Y->size2;
    int v = P->size;

    double thetajh,phih;
    int zh;
    double Djh;
    gsl_vector_view Vview;

    gsl_matrix * YmLambdaEEta = gsl_matrix_alloc(p,n);
    gsl_matrix_memcpy(YmLambdaEEta,Y);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1.0,Lambda,EEta,1.0,YmLambdaEEta);


    for(int h = 0; h<k; h++){
        int pcul = 0;
        int pw;
        Vview = gsl_matrix_column(Lambda,h);
        gsl_vector_const_view cVview = gsl_matrix_const_row(EEta,h);
        gsl_blas_dger(1.0,&Vview.vector,&cVview.vector,YmLambdaEEta);

        for(int w = 0; w<v; w++){
            pw = gsl_vector_int_get(P,w);
            phih = gsl_matrix_get(Phi,w,h);
            zh = gsl_matrix_get(EZ,w,h);
            gsl_matrix_view Lambdaw = gsl_matrix_submatrix(Lambda,pcul,0,pw,k);
            Vview = gsl_matrix_column(&Lambdaw.matrix,h);
            if(zh == 0.0 && phih == 0.0)
                gsl_vector_set_zero(&Vview.vector);
            else{
                for(int j = pcul; j<pw+pcul; j++){
                    double top;
                    Vview = gsl_matrix_row(YmLambdaEEta,j);
                    gsl_blas_ddot(&Vview.vector,&cVview.vector,&top);

                    thetajh = gsl_matrix_get(Theta,j,h);
                    if(thetajh == 0.0){
                        gsl_matrix_set(Lambda,j,h,0.0);
                    }
                    else{
                        Djh = zh/thetajh + (1-zh)/phih;
                        Djh = Djh/gsl_vector_get(Sigma2inv,j) + gsl_matrix_get(EEtaTEta,h,h);
                        gsl_matrix_set(Lambda,j,h,top/Djh);
                    }
                }
            }
            pcul += pw;
        }
        Vview = gsl_matrix_column(Lambda,h);
        gsl_blas_dger(-1.0,&Vview.vector,&cVview.vector,YmLambdaEEta);

    }
    gsl_matrix_free(YmLambdaEEta);
}



// make sure that the empty factor check is before this step
// SigmaEta Lower triangle is stored EEtaTEta = EtaTEta + n*SigmaEta,

void updateEta(gsl_matrix * EEta, gsl_matrix * SigmaEta,
               const gsl_matrix * Lambda,
               const gsl_matrix * Y, const gsl_vector * Sigma2inv)
{
    int p = Y->size1;
    int n = Y->size2;
    int k = EEta->size1;

    // calculate Lambda * Sigma2inv
    gsl_vector_view Vview;
    gsl_matrix * LambdaTSigma = gsl_matrix_alloc(k,p);
    gsl_matrix_transpose_memcpy(LambdaTSigma,Lambda);
    for(int h = 0; h<k; h++){
        Vview = gsl_matrix_row(LambdaTSigma,h);
        gsl_vector_mul(&Vview.vector,Sigma2inv);
    }
    // calculate SigmaEta
    gsl_matrix_set_identity(SigmaEta);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,LambdaTSigma,Lambda,1.0,SigmaEta);

    gsl_matrix * SigmaEtainv = gsl_matrix_alloc(k,k);
    // cholesky decomposition
    gsl_linalg_cholesky_decomp(SigmaEta);
    gsl_matrix_memcpy(SigmaEtainv,SigmaEta);
    gsl_linalg_cholesky_invert(SigmaEta);

    // calculate LambdaT*Sigma2Inv*Y
    gsl_matrix * LambdaTSigmaY = gsl_matrix_alloc(k,n);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,LambdaTSigma,Y,0.0,LambdaTSigmaY);

    gsl_vector_view Vview2;
    for(int i = 0; i<n; i++){
        // calculate MuEta
        Vview = gsl_matrix_column(LambdaTSigmaY,i);
        //gsl_linalg_LU_solve(SigmaEtaNonZeroinv,Perm,&Vview.vector,MuEtaNonZero);
        Vview2 = gsl_matrix_column(EEta,i);
        gsl_linalg_cholesky_solve(SigmaEtainv,&Vview.vector,&Vview2.vector);
    }

    gsl_matrix_free(LambdaTSigma);
    gsl_matrix_free(SigmaEtainv);
    gsl_matrix_free(LambdaTSigmaY);
    //gsl_permutation_free(Perm);
}


// make sure that the empty factor check is before this step
// Lower part of SigmaEta is stored EEtaTEta = EtaTEta + n*SigmaEta;
void updateEtaPX(gsl_matrix * EEta, gsl_matrix * SigmaEta,
                 const gsl_matrix * Lambda,
                 const gsl_matrix * Y, const gsl_vector * Sigma2inv)
{
    int p = Y->size1;
    int n = Y->size2;
    int k = EEta->size1;

    // PXL-EM to calculate Lambda*A_L
    gsl_matrix * LambdaNew = gsl_matrix_alloc(p,k);
    gsl_matrix_memcpy(LambdaNew,Lambda);

    //gsl_matrix_set_identity(SigmaEta);
    //gsl_blas_dsyrk(CblasLower,CblasNoTrans,1.0/n,EEta,1.0,SigmaEta);

    gsl_matrix_scale(SigmaEta,1.0/n);
    gsl_linalg_cholesky_decomp(SigmaEta);
    gsl_blas_dtrmm(CblasRight,CblasUpper,CblasTrans,CblasNonUnit,1.0,SigmaEta,LambdaNew);


    // calculate Lambda * Sigma2inv
    gsl_vector_view Vview;
    gsl_matrix * LambdaTSigma = gsl_matrix_alloc(k,p);
    gsl_matrix_transpose_memcpy(LambdaTSigma,LambdaNew);
    for(int h = 0; h<k; h++){
        Vview = gsl_matrix_row(LambdaTSigma,h);
        gsl_vector_mul(&Vview.vector,Sigma2inv);
    }
    // calculate SigmaEta
    gsl_matrix_set_identity(SigmaEta);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,LambdaTSigma,LambdaNew,1.0,SigmaEta);

    gsl_matrix * SigmaEtainv = gsl_matrix_alloc(k,k);
    // cholesky decomposition
    gsl_linalg_cholesky_decomp(SigmaEta);
    gsl_matrix_memcpy(SigmaEtainv,SigmaEta);
    gsl_linalg_cholesky_invert(SigmaEta);

    // calculate LambdaT*Sigma2Inv*Y
    gsl_matrix * LambdaTSigmaY = gsl_matrix_alloc(k,n);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,LambdaTSigma,Y,0.0,LambdaTSigmaY);

    gsl_vector_view Vview2;
    for(int i = 0; i<n; i++){
        // calculate MuEta
        Vview = gsl_matrix_column(LambdaTSigmaY,i);
        //gsl_linalg_LU_solve(SigmaEtaNonZeroinv,Perm,&Vview.vector,MuEtaNonZero);
        Vview2 = gsl_matrix_column(EEta,i);
        gsl_linalg_cholesky_solve(SigmaEtainv,&Vview.vector,&Vview2.vector);
    }

    gsl_matrix_free(LambdaTSigma);
    gsl_matrix_free(SigmaEtainv);
    gsl_matrix_free(LambdaTSigmaY);
    gsl_matrix_free(LambdaNew);
    //gsl_permutation_free(Perm);
}


// no empty factor at this step too
void updateThetaDeltaPhi(gsl_matrix * Theta, gsl_matrix * Delta, gsl_vector * Phi,
                         const gsl_matrix * Lambda, const gsl_vector * Tau,
                         const gsl_vector * EZ, double a, double b, double c)
{
    int p = Theta->size1;
    int k = Theta->size2;

    double rgigP,rgigA,rgigB;
    double Eg,Eig,Elogg;
    double zh,phih;

    double gammaA,gammaB;
    double sumDeltah,sumLambda2h;
    double thetajh,deltajh;

    for(int h = 0; h<k; h++){
        phih = gsl_vector_get(Phi,h);
        if(phih == 0.0){
            continue;
        }
        zh = gsl_vector_get(EZ,h);
        rgigP = (a - 0.5);

        // for Delta
        sumDeltah = 0.0;
        sumLambda2h = 0.0;
        int numZero = 0;

        for(int j = 0; j<p; j++){
            // For Theta
            rgigA = gsl_matrix_get(Delta,j,h) * 2.0;
            rgigB = gsl_matrix_get(Lambda,j,h) * gsl_matrix_get(Lambda,j,h);
            //rgigP = rgigP * zh;
            //rgigB = rgigB * zh;
            //rgigA = rgigA * zh;
            //Egig(1,&rgigP,&rgigA,&rgigB,&Eg,&Eig,&Elogg);
            thetajh = ((rgigP - 1.0) + sqrt((rgigP - 1.0)*(rgigP - 1.0) + rgigA * rgigB))/rgigA;
            //thetajh = (rgigP + sqrt(rgigP * rgigP + rgigA * rgigB))/rgigA;

            // set thetajh to zero when too small
            // possible generating new empty factors
            if(thetajh < 1e-300){
                thetajh = 0.0;
                numZero ++;
            }
            gsl_matrix_set(Theta,j,h,thetajh);
            // For Delta
            gammaA = (a+b);
            gammaB = (phih + thetajh);
            deltajh = gammaA/gammaB;
            gsl_matrix_set(Delta,j,h,deltajh);
            //gsl_matrix_set(ElogDelta,j,h,gsl_sf_psi(gammaA*zh) - log(gammaB*zh));
            //gsl_matrix_set(ElogDelta,j,h,log(gammaA) - log(gammaB));

            sumDeltah += deltajh;
            sumLambda2h += rgigB;
        }
        if(numZero == p)
            phih = 0.0;
        else{
            rgigP = zh*p*b - (1.0-zh)*p*0.5 + c;
            rgigA = 2.0*(gsl_vector_get(Tau,h) + zh*sumDeltah);
            rgigB = (1.0-zh)*sumLambda2h;

            //Egig(1,&rgigP,&rgigA,&rgigB,&Eg,&Eig,&Elogg);
            phih = ((rgigP - 1.0) + sqrt((rgigP - 1.0)*(rgigP - 1.0) + rgigA * rgigB))/rgigA;
            if(phih < 1e-300)
                phih = 1e-300;
            //phih = (rgigP + sqrt(rgigP * rgigP + rgigA * rgigB))/rgigA;
        }
        gsl_vector_set(Phi,h,phih);
    }
}


void updateTauHGamma(gsl_vector * Tau, double * H, double * Gamma,
                     const gsl_vector * Phi, double c, double d, double e,
                     double f,double nu)
{
    int k = Tau->size;
    double gammaA,gammaB;
    gammaA = c+d;
    double sumTau = 0.0;
    double tauh;
    for(int h = 0; h<k; h++){
        gammaB = gsl_vector_get(Phi,h) + *H;
        tauh = gammaA/gammaB;
        gsl_vector_set(Tau,h,tauh);
        sumTau += tauh;
    }

    gammaA = k*d + e;
    gammaB = *Gamma + sumTau;
    *H = gammaA/gammaB;

    gammaA = e + f;
    gammaB = *H + nu;
    *Gamma = gammaA/gammaB;
}

void updatePiZCollapse(double * Pi, gsl_vector * EZ,
                       const gsl_matrix * Lambda, const gsl_matrix * Theta,
                       const gsl_vector * Phi, double a, double b)
{
    int k = EZ->size;
    int p = Lambda->size1;

    double sumCol0,sumCol1;
    double logPmax;
    double sumEZ = 0.0;
    double nx;
    double sigma0,sigma1;
    double dnorm;

    double prob;

    for(int h = 0; h<k; h++){
        sigma0 = gsl_vector_get(Phi,h);
        if(sigma0 == 0.0){
            gsl_vector_set(EZ,h,0.0);
            continue;
        }
        sumCol0 = 0;
        sumCol1 = 0;
        for(int j = 0; j<p; j++){
            // Theta = 0, Lambda = 0, then the density is inf, then the factor is sparse
            if(gsl_matrix_get(Theta,j,h) == 0.0){
                sumCol1 = INFINITY;
                break;
            }
            nx = gsl_matrix_get(Lambda,j,h);
            dnorm = -0.5*gsl_sf_log(2*M_PI*sigma0) - nx*nx/(2.0*sigma0);
            sumCol0 += dnorm;

            sigma1 = gsl_matrix_get(Theta,j,h);
            dnorm = -0.5*gsl_sf_log(2*M_PI*sigma1) - nx*nx/(2.0*sigma1);
            sumCol1 += dnorm;

            sumCol1 += (a-1)*gsl_sf_log(sigma1);
            sumCol1 -= (a+b)*gsl_sf_log(sigma1 + sigma0);
        }
        sumCol0 += gsl_sf_log(1 - *Pi);
        sumCol1 += gsl_sf_log(*Pi);
        sumCol1 += b*p*gsl_sf_log(sigma0);
        sumCol1 += p*(gsl_sf_lngamma(a+b) - gsl_sf_lngamma(a) - gsl_sf_lngamma(b));

        logPmax = max(sumCol0,sumCol1);
        sumCol0 -= logPmax;
        sumCol1 -= logPmax;
        if(sumCol0 < -500.0){
            sumCol0 = 0.0;
            sumCol1 = 1.0;
        }
        else if(sumCol1 < -500.0){
            sumCol0 = 1.0;
            sumCol1 = 0.0;
        }
        else{
            sumCol0 = gsl_sf_exp(sumCol0);
            sumCol1 = gsl_sf_exp(sumCol1);
        }
        prob = sumCol1/(sumCol1 + sumCol0);
        sumEZ += prob;
        gsl_vector_set(EZ,h,prob);
        //gsl_vector_set(EZ,h,0.99);
    }

    *Pi = (sumEZ+1)/(k+2);
}


void updatePiZnonCollapse(double *Pi, gsl_vector * EZ,
                          const gsl_matrix * Lambda, const gsl_matrix * Theta,
                          const gsl_matrix * Delta,
                          const gsl_vector * Phi, double a, double b)
{
    int k = EZ->size;
    int p = Lambda->size1;
    int Z1n, Z0n;
    Z1n = 0;

    double sigma0,sigma1,dnorm,dgamma,sumCol0,sumCol1;
    double ga,gb,gx;
    double nx;
    double logPmax, PZ1;
    int Zh;
    Z1n = 0;
    double gammaR,gammaA,gammaB;
    gammaA = a+b;
    double sumEZ = 0.0;
    for(int h = 0; h<k; h++){
        sigma0 = gsl_vector_get(Phi,h);
        if(sigma0 == 0.0){
            gsl_vector_set(EZ,h,0.0);
            continue;
        }
        sumCol0 = 0;
        sumCol1 = 0;
        for(int j = 0; j<p; j++){

            if(gsl_matrix_get(Theta,j,h) == 0.0){
                sumCol1 = INFINITY;
                break;
            }

            nx = gsl_matrix_get(Lambda,j,h);
            dnorm = -0.5*gsl_sf_log(2*M_PI*sigma0) - nx*nx/(2.0*sigma0);
            sumCol0 += dnorm;

            sigma1 = gsl_matrix_get(Theta,j,h);
            dnorm = -0.5*gsl_sf_log(2*M_PI*sigma1) - nx*nx/(2.0*sigma1);
            sumCol1 += dnorm;

            ga = a;
            gb = gsl_matrix_get(Delta,j,h);
            gx = gsl_matrix_get(Theta,j,h);
            dgamma = ga*gsl_sf_log(gb) - gsl_sf_lngamma(ga) +
                    (ga-1)*gsl_sf_log(gx) - gx*gb;
            sumCol1 += dgamma;

            ga = b;
            gb = gsl_vector_get(Phi,h);
            gx = gsl_matrix_get(Delta,j,h);
            dgamma = ga*gsl_sf_log(gb) - gsl_sf_lngamma(ga) +
                    (ga-1)*gsl_sf_log(gx) - gx*gb;
            sumCol1 += dgamma;
        }
        sumCol0 += gsl_sf_log(1 - *Pi);
        sumCol1 += gsl_sf_log(*Pi);

        logPmax = max(sumCol0,sumCol1);
        sumCol0 -= logPmax;
        sumCol1 -= logPmax;
        if(sumCol0 < -500.0){
            sumCol0 = 0.0;
            sumCol1 = 1.0;
        }
        else if(sumCol1 < -500.0){
            sumCol0 = 1.0;
            sumCol1 = 0.0;
        }
        else{
            sumCol0 = gsl_sf_exp(sumCol0);
            sumCol1 = gsl_sf_exp(sumCol1);
        }
        PZ1 = sumCol1/(sumCol1 + sumCol0);
        sumEZ += PZ1;
        gsl_vector_set(EZ,h,PZ1);
        Z1n += Zh;

    }
    Z0n = k - Z1n;

    *Pi = (sumEZ+1)/(k+2);

}

void updateSigma2inv(gsl_vector * Sigma2inv,
                     const gsl_matrix * Lambda,
                     const gsl_matrix * Eta, const gsl_matrix * Y,
                     const double a, const double b)
{
    int p = Sigma2inv->size;
    int k = Lambda->size2;
    int n = Y->size2;

    gsl_matrix * Ymean = gsl_matrix_alloc(p,n);
    mulM1M2(Ymean,Lambda,Eta);
    gsl_matrix_sub(Ymean,Y);

    double gammaA,gammaB;
    gammaA = n*.5+a;
    gsl_vector_view Vview;
    for(int j = 0; j<p; j++){
        gammaB = 0.0;
        Vview = gsl_matrix_row(Ymean,j);
        gsl_blas_ddot(&Vview.vector,&Vview.vector,&gammaB);
        gammaB = gammaB*.5 + b;
        gsl_vector_set(Sigma2inv,j,gammaA/gammaB);
    }
    gsl_matrix_free(Ymean);
}


}


/* ----------------------------------------------------------- */
/* ------------------- Handling factors ---------------------- */
/* ----------------------------------------------------------- */


int deleteNullFactors(gsl_matrix *& Lambda,
                      gsl_matrix *& Eta, gsl_matrix *& SigmaEta,
                      gsl_matrix *& EtaTEta, gsl_matrix *& EtaTY,
                      gsl_matrix *& Theta, gsl_matrix *& Delta,
                      gsl_matrix *& Phi, gsl_matrix *& Tau, gsl_matrix *& Z)
{
    int p = Lambda->size1;
    int k = Lambda->size2;
    int v = Phi->size1;
    int n = Eta->size2;

    int keepk = 0;
    int keepInd[k];
    for(int h = 0; h<k; h++){
        int sum = 0;
        for(int w=0; w<v; w++){
            sum += gsl_matrix_get(Phi,w,h) == 0.0;
        }

        if(sum != v){
            keepInd[keepk] = h;
            keepk ++;
        }
    }
    if(keepk == k){
        return 0;
    }
    else if(keepk == 0){
        cout<<"Loading matrix is empty, please re-run"<<endl;
        exit(EXIT_FAILURE);
    }

    gsl_vector_view Vview1;
    gsl_vector_view Vview2;
    gsl_matrix * EtaNew = gsl_matrix_alloc(keepk,n);
    gsl_matrix * SigmaEtaNew = gsl_matrix_calloc(keepk,keepk);
    gsl_matrix * EtaTEtaNew = gsl_matrix_alloc(keepk,keepk);
    gsl_matrix * EtaTYNew = gsl_matrix_alloc(keepk,p);

    gsl_matrix * LambdaNew = gsl_matrix_alloc(p,keepk);
    gsl_matrix * ThetaNew = gsl_matrix_alloc(p,keepk);
    gsl_matrix * DeltaNew = gsl_matrix_alloc(p,keepk);
    gsl_matrix * PhiNew = gsl_matrix_alloc(v,keepk);
    gsl_matrix * TauNew = gsl_matrix_alloc(v,keepk);
    gsl_matrix * ZNew = gsl_matrix_alloc(v,keepk);

    int preh;
    for(int h = 0; h<keepk; h++){
        preh = keepInd[h];
        // Lambda
        Vview1 = gsl_matrix_column(Lambda,preh);
        Vview2 = gsl_matrix_column(LambdaNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        // Theta
        Vview1 = gsl_matrix_column(Theta,preh);
        Vview2 = gsl_matrix_column(ThetaNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        // Delta
        Vview1 = gsl_matrix_column(Delta,preh);
        Vview2 = gsl_matrix_column(DeltaNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        // Phi
        Vview1 = gsl_matrix_column(Phi,preh);
        Vview2 = gsl_matrix_column(PhiNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        // Tau
        Vview1 = gsl_matrix_column(Tau,preh);
        Vview2 = gsl_matrix_column(TauNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        // Z
        Vview1 = gsl_matrix_column(Z,preh);
        Vview2 = gsl_matrix_column(ZNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        //Eta
        Vview1 = gsl_matrix_row(Eta,preh);
        Vview2 = gsl_matrix_row(EtaNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        //EtaTY
        Vview1 = gsl_matrix_row(EtaTY,preh);
        Vview2 = gsl_matrix_row(EtaTYNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        // SigmaEta, EtaTEta
        for(int h2 = 0; h2<keepk; h2++){
            int preh2 = keepInd[h2];
            gsl_matrix_set(SigmaEtaNew,h,h2,gsl_matrix_get(SigmaEta,preh,preh2));
            gsl_matrix_set(EtaTEtaNew,h,h2,gsl_matrix_get(EtaTEta,preh,preh2));
        }

    }
    gsl_matrix_free(Lambda); Lambda = LambdaNew;
    gsl_matrix_free(Theta); Theta = ThetaNew;
    gsl_matrix_free(Delta); Delta = DeltaNew;
    gsl_matrix_free(Phi); Phi = PhiNew;
    gsl_matrix_free(Tau); Tau = TauNew;
    gsl_matrix_free(Z); Z = ZNew;

    gsl_matrix_free(Eta); Eta = EtaNew;
    gsl_matrix_free(EtaTEta); EtaTEta = EtaTEtaNew;
    gsl_matrix_free(EtaTY); EtaTY = EtaTYNew;
    gsl_matrix_free(SigmaEta); SigmaEta = SigmaEtaNew;

    return 1;

}


// delete the loading that has only one element in it.
void pruneFactors(gsl_matrix *& Lambda,
                  gsl_matrix *& Eta, gsl_matrix *& SigmaEta,
                  gsl_matrix *& EtaTEta, gsl_matrix *& EtaTY,
                  gsl_matrix *& Theta, gsl_matrix *& Delta,
                  gsl_matrix *& Phi, gsl_matrix *& Tau,
                  gsl_matrix *& Z, gsl_vector *& Sigma2inv)
{
    int p = Lambda->size1;
    int k = Lambda->size2;
    int v = Phi->size1;
    int n = Eta->size2;

    int keepk = 0;
    int deletek = 0;
    int keepInd[k];
    int deleteInd[k];
    for(int h = 0; h<k; h++){
        int sum = 0;
        for(int j=0; j<p; j++){
            sum += gsl_matrix_get(Lambda,j,h) == 0.0;
        }

        if(sum != p-1){
            keepInd[keepk] = h;
            keepk ++;
        }
        else{
            deleteInd[deletek] = h;
            deletek ++;
        }
    }
    if(keepk == k){
        return;
    }
    else if(keepk == 0){
        cout<<"Loading matrix is empty, please re-run"<<endl;
        exit(EXIT_FAILURE);
    }

    gsl_vector_view Vview1;
    gsl_vector_view Vview2;
    gsl_matrix * EtaNew = gsl_matrix_alloc(keepk,n);
    gsl_matrix * SigmaEtaNew = gsl_matrix_calloc(keepk,keepk);
    gsl_matrix * EtaTEtaNew = gsl_matrix_alloc(keepk,keepk);
    gsl_matrix * EtaTYNew = gsl_matrix_alloc(keepk,p);

    gsl_matrix * LambdaNew = gsl_matrix_alloc(p,keepk);
    gsl_matrix * ThetaNew = gsl_matrix_alloc(p,keepk);
    gsl_matrix * DeltaNew = gsl_matrix_alloc(p,keepk);
    gsl_matrix * PhiNew = gsl_matrix_alloc(v,keepk);
    gsl_matrix * TauNew = gsl_matrix_alloc(v,keepk);
    gsl_matrix * ZNew = gsl_matrix_alloc(v,keepk);

    int preh;
    for(int h = 0; h<keepk; h++){
        preh = keepInd[h];
        // Lambda
        Vview1 = gsl_matrix_column(Lambda,preh);
        Vview2 = gsl_matrix_column(LambdaNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        // Theta
        Vview1 = gsl_matrix_column(Theta,preh);
        Vview2 = gsl_matrix_column(ThetaNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        // Delta
        Vview1 = gsl_matrix_column(Delta,preh);
        Vview2 = gsl_matrix_column(DeltaNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        // Phi
        Vview1 = gsl_matrix_column(Phi,preh);
        Vview2 = gsl_matrix_column(PhiNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        // Tau
        Vview1 = gsl_matrix_column(Tau,preh);
        Vview2 = gsl_matrix_column(TauNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        // Z
        Vview1 = gsl_matrix_column(Z,preh);
        Vview2 = gsl_matrix_column(ZNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        //Eta
        Vview1 = gsl_matrix_row(Eta,preh);
        Vview2 = gsl_matrix_row(EtaNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        //EtaTY
        Vview1 = gsl_matrix_row(EtaTY,preh);
        Vview2 = gsl_matrix_row(EtaTYNew,h);
        gsl_vector_memcpy(&Vview2.vector,&Vview1.vector);

        // SigmaEta, EtaTEta
        for(int h2 = 0; h2<keepk; h2++){
            int preh2 = keepInd[h2];
            gsl_matrix_set(SigmaEtaNew,h,h2,gsl_matrix_get(SigmaEta,preh,preh2));
            gsl_matrix_set(EtaTEtaNew,h,h2,gsl_matrix_get(EtaTEta,preh,preh2));
        }

    }


    for(int h = 0; h<deletek; h++){
        preh = deleteInd[h];
        double sum2 = 0;
        int j;
        for(j = 0; j<p; j++){
            if(gsl_matrix_get(Lambda,j,preh)!=0.0){
                sum2 = gsl_matrix_get(Lambda,j,preh)*gsl_matrix_get(Lambda,j,preh);
                break;
            }
        }
        double sigma2j = 1.0/gsl_vector_get(Sigma2inv,j);
        gsl_vector_set(Sigma2inv,j,1.0/(sigma2j+sum2));
    }

    gsl_matrix_free(Lambda); Lambda = LambdaNew;
    gsl_matrix_free(Theta); Theta = ThetaNew;
    gsl_matrix_free(Delta); Delta = DeltaNew;
    gsl_matrix_free(Phi); Phi = PhiNew;
    gsl_matrix_free(Tau); Tau = TauNew;
    gsl_matrix_free(Z); Z = ZNew;

    gsl_matrix_free(Eta); Eta = EtaNew;
    gsl_matrix_free(EtaTEta); EtaTEta = EtaTEtaNew;
    gsl_matrix_free(EtaTY); EtaTY = EtaTYNew;
    gsl_matrix_free(SigmaEta); SigmaEta = SigmaEtaNew;
}

int countNZeroLoadingElement(const gsl_matrix * Lambda)
{
    int num = 0;
    int p = Lambda->size1;
    int k = Lambda->size2;
    for(int j = 0; j<p; j++){
        for(int h = 0; h<k; h++){
            if(gsl_matrix_get(Lambda,j,h) != 0.0)
                num++;
        }
    }
    return num;
}

double calculateMatrixNorm(const gsl_matrix * Lambda)
{
    double sum = 0.0;
    int p = Lambda->size1;
    int k = Lambda->size2;
    for(int j = 0; j<p; j++){
        for(int h = 0; h<k; h++){
            sum += gsl_matrix_get(Lambda,j,h)*gsl_matrix_get(Lambda,j,h);
        }
    }
    return sum;
}


double calculateLogLikelihood(const gsl_matrix * Y,
                              const gsl_matrix * Lambda, const gsl_matrix * Eta,
                              const gsl_vector * Sigma2inv)
{
    int p = Y->size1;
    int n = Y->size2;

    double loglikelihood = 0.0;

    gsl_vector_view Vview;
    gsl_matrix * Ymean = gsl_matrix_alloc(p,n);
    mulM1M2(Ymean,Lambda,Eta);
    gsl_matrix_sub(Ymean,Y);
    double Sigmainvj;
    for(int j = 0; j<p; j++){
        Sigmainvj = sqrt(gsl_vector_get(Sigma2inv,j));
        Vview = gsl_matrix_row(Ymean,j);
        gsl_vector_scale(&Vview.vector,Sigmainvj);

        loglikelihood += Sigmainvj;
    }

    double dnorm;
    for(int i = 0; i<n; i++){
        Vview = gsl_matrix_column(Ymean,i);
        dnorm = gsl_blas_dnrm2(&Vview.vector);
        loglikelihood -= 0.5 * dnorm * dnorm;
    }

    gsl_matrix_free(Ymean);
    return loglikelihood;
}
