#ifndef FUNWRAPPER_H
#define FUNWRAPPER_H

#include "gig_par.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>


using namespace std;

/* ------------------------- ranodm number generator funs -------------------------- */
int randomGamma(double * rv, int lrv, const double * a, int la,
                const double * b, int lb, gsl_rng * r);
int randomUniNormal(double * rv, int lrv, const double * mu, int lm,
                    const double * sigma, int ls, gsl_rng *r);
int randomGIG(double * rv, int lrv, const double * p, int lp,
              const double * a, int la, const double * b, int lb,
              gig * mygig);
int randomBernoulli(double * rv, int lrv, const double * p, int lp,
                    gsl_rng * r);


/* -------------------------- matrix operation wrapper ----------------------------- */
int mulM1M2(gsl_matrix * Dest, const gsl_matrix * M1, const gsl_matrix * M2);
int mulM1V1byrow(gsl_matrix * M1, const gsl_vector * V1);
int divM1V1byrow(gsl_matrix * M1, const gsl_vector * V1);
int mulM1V1bycol(gsl_matrix * M1, const gsl_vector * V1);
int divM1V1bycol(gsl_matrix * M1, const gsl_vector * V1);
int mulM1diagV1M2(gsl_matrix * Dest, const gsl_matrix * M1,
                  const gsl_vector * V1, const gsl_matrix * M2);
int addM1V1diag(gsl_matrix * M1, const gsl_vector * V1);



/* ------------------------------------------------------------------ */
/* -----------------  namespace for MCMC functions ------------------ */
/* ------------------------------------------------------------------ */
namespace MCMCFuns{

void updateLambda(gsl_matrix * Lambda,
                  const gsl_matrix * EtaTEta, const gsl_matrix * EtaTY,
                  const gsl_matrix * Theta, const gsl_vector * Phi,
                  const gsl_vector * Z, const gsl_vector * Sigma2inv,
                  gsl_rng * r);

void updateEta(gsl_matrix * Eta, gsl_matrix *SigmaEta,
               const gsl_matrix * Lambda, const gsl_matrix * Y,
               const gsl_vector * Sigma2inv,
               gsl_rng * r);

void updateThetaDeltaPhi(gsl_matrix * Theta, gsl_matrix * Delta,
                         gsl_vector * Phi,
                         const gsl_matrix * Lambda, const gsl_vector * Tau,
                         const gsl_vector * Z, double a, double b, double c,
                         gsl_rng * r, gig * mygig);


void updateTauHGamma(gsl_vector * Tau, double *H, double *Gamma,
                     const gsl_vector * Phi, double c, double d, double e,
                     double f, double nu, gsl_rng * r);


void updatePiZCollapse(double * Pi, gsl_vector * Z,
                       const gsl_matrix * Lambda, gsl_matrix * Theta,
                       gsl_matrix * Delta, const gsl_vector * Phi,
                       double a, double b, gsl_rng * r);

void updatePiZnonCollapse(double * Pi, gsl_vector * Z,
                          const gsl_matrix * Lambda, gsl_matrix * Theta,
                          gsl_matrix * Delta, const gsl_vector * Phi,
                          double a, double b, gsl_rng * r);

void updateSigma2inv(gsl_vector * Sigma2inv,
                     const gsl_matrix * Lambda,
                     const gsl_matrix * Eta, const gsl_matrix * Y,
                     const double a, const double b,
                     gsl_rng * r);

}


/* ----------------------------------------------------------------- */
/* ---------------- namespace for EM functions --------------------- */
/* ----------------------------------------------------------------- */

namespace EMFuns{

void updateLambda(gsl_matrix * Lambda, const gsl_vector_int *P, const gsl_matrix * Theta,
                  const gsl_matrix *Phi, const gsl_matrix *EZ,
                  const gsl_matrix * EEta, const gsl_matrix * Y,
                  const gsl_matrix * EEtaTEta, const gsl_vector * Sigma2inv);

// make sure that the empty factor check is before this step
void updateEta(gsl_matrix * EEta, gsl_matrix * SigmaEta,
               const gsl_matrix * Lambda,
               const gsl_matrix * Y, const gsl_vector * Sigma2inv);

// make sure that the empty factor check is before this step
void updateEtaPX(gsl_matrix * EEta, gsl_matrix * SigmaEta,
               const gsl_matrix * Lambda,
               const gsl_matrix * Y, const gsl_vector * Sigma2inv);

// no empty factor at this step too
void updateThetaDeltaPhi(gsl_matrix * Theta, gsl_matrix * Delta, gsl_vector * Phi,
                        const gsl_matrix * Lambda, const gsl_vector * Tau,
                        const gsl_vector * EZ, double a, double b, double c);

void updateTauHGamma(gsl_vector * Tau, double *H, double *Gamma,
                     const gsl_vector * Phi, double c, double d, double e,
                     double f, double nu);

void updatePiZCollapse(double *Pi, gsl_vector * EZ,
                       const gsl_matrix * Lambda, const gsl_matrix * Theta,
                       const gsl_vector * Phi, double a, double b);

void updatePiZnonCollapse(double * Pi, gsl_vector * EZ,
                          const gsl_matrix * Lambda, const gsl_matrix * Theta,
                          const gsl_matrix * Delta,
                          const gsl_vector * Phi, double a, double b);

void updateSigma2inv(gsl_vector * Sigma2inv,
                     const gsl_matrix * Lambda,
                     const gsl_matrix * Eta, const gsl_matrix * Y,
                     const double a, const double b);

}


/* ----------------------------------------------------------- */
/* ------------------- Handling factors ---------------------- */
/* ----------------------------------------------------------- */

int deleteNullFactors(gsl_matrix *& Lambda,
                      gsl_matrix *& Eta, gsl_matrix *& SigmaEta,
                      gsl_matrix *& EtaTEta, gsl_matrix *& EtaTY,
                      gsl_matrix *& Theta, gsl_matrix *& Delta,
                      gsl_matrix *& Phi, gsl_matrix *& Tau, gsl_matrix *& Z);

void pruneFactors(gsl_matrix *& Lambda,
                 gsl_matrix *& Eta, gsl_matrix *& SigmaEta,
                  gsl_matrix *& EtaTEta, gsl_matrix *& EtaTY,
                 gsl_matrix *& Theta, gsl_matrix *& Delta,
                 gsl_matrix *& Phi, gsl_matrix *& Tau,
                 gsl_matrix *& Z, gsl_vector *& Sigma2inv);

int countNZeroLoadingElement(const gsl_matrix * Lambda);

double calculateMatrixNorm(const gsl_matrix * Lambda);

double calculateLogLikelihood(const gsl_matrix * Y,
                              const gsl_matrix * Lambda, const gsl_matrix * Eta,
                              const gsl_vector * Sigma2inv);

#endif // FUNWRAPPER_H
