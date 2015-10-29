#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <iostream>
#include <sstream>
#include <fstream>

#include "BASS.h"

using namespace std;

int main(int argc, char * argv[])
{
    const int maxV = 100;
    string fileY;
    int pArray[maxV];
    int n;
    int v;
    int k;
    int S, step;
    int MCMCS;
    int PXS;
    string outputDir;
    string sep;

    string s_y="-y", s_out="-out",s_itr="-itr";
    string s_k="-k", s_n="-n", s_v="-v";
    string s_pArray[maxV];
    for(int i = 1; i<=maxV; i++){
        stringstream ss;
        ss<<i;
        s_pArray[i-1] = "-p" + ss.str();
    }
    string s_mcmcs="-mcmc";
    string s_step="-step";
    string s_sep="-sep";
    string s_pxs="-px";
    for(int i = 0; i<argc; i++){
        if(s_y.compare(argv[i])==0) {fileY = argv[i+1]; cout<<"fileY: "<<fileY<<endl;}
        if(s_sep.compare(argv[i])==0) {sep = argv[i+1]; cout<<"sep: "<<sep<<endl;}
        if(s_k.compare(argv[i])==0) {k = atoi(argv[i+1]); cout<<"k: "<<k<<endl;}
        if(s_v.compare(argv[i])==0) {v = atoi(argv[i+1]); cout<<"v: "<<v<<endl;}
        for(int w=0; w<maxV; w++){
            if(s_pArray[w].compare(argv[i])==0)
            {pArray[w] = atoi(argv[i+1]); cout<<s_pArray[w] + ": "<<pArray[w]<<endl;}
        }
        if(s_n.compare(argv[i])==0) {n = atoi(argv[i+1]); cout<<"n: "<<n<<endl;}
        if(s_itr.compare(argv[i])==0) {S = atoi(argv[i+1]); cout<<"S: "<<S<<endl;}
        if(s_out.compare(argv[i])==0) {outputDir = argv[i+1]; cout<<"outputDir: "<<outputDir<<endl;}
        if(s_step.compare(argv[i])==0) {step = atoi(argv[i+1]); cout<<"step: "<<step<<endl;}
        if(s_mcmcs.compare(argv[i])==0) {MCMCS = atoi(argv[i+1]); cout<<"MCMCS: "<<MCMCS<<endl;}
        if(s_pxs.compare(argv[i])==0) {PXS = atoi(argv[i+1]); cout<<"PXS: "<<PXS<<endl;}
    }
    if(PXS+MCMCS>=S){ cout<<"please increase iteration number"<<endl; exit(EXIT_FAILURE);}

    int cs = MCMCS/2;

    /*
    fileY = "RData/Ystd2.sim";
    sep = "space";
    pArray[0] = 70;
    pArray[1] = 50;
    v = 2;
    n = 30;
    k = 7;
    S = 500;
    MCMCS = 200;
    cS = MCMCS-10;
    step = 10;
    outputDir = "cppRunTest/rep1";

    string outputDir2 = outputDir + "/TraceMCMC";
    string outputDir3 = outputDir + "/TraceEM";
    */

    // make output directory
    // ------------------------------------------

    int tag = mkdir(outputDir.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(tag != 0){
        cout<<"Error in creating output directory." <<endl;
        cout<<"Possible reasons:"<<endl;
        cout<<"   it has been created;" <<endl;
        cout<<"   its parent directory does not exist; "<<endl;
        exit(EXIT_FAILURE);
    }

    /*
    tag = mkdir(outputDir2.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(tag != 0){
        cout<<"Error in creating output directory." <<endl;
        cout<<"Possible reasons:"<<endl;
        cout<<"   it has been created;" <<endl;
        cout<<"   its parent directory does not exist; "<<endl;
        exit(EXIT_FAILURE);
    }

    tag = mkdir(outputDir3.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(tag != 0){
        cout<<"Error in creating output directory." <<endl;
        cout<<"Possible reasons:"<<endl;
        cout<<"   it has been created;" <<endl;
        cout<<"   its parent directory does not exist; "<<endl;
        exit(EXIT_FAILURE);
    }
    */

    // construct BASS class
    // ------------------------------------------
    BASS myBASS(pArray,v,n,k);
    myBASS.LoadData(fileY,sep);
    myBASS.InitializePara();

    // ---------------------------------------------
    // -------------- Write info to file -----------
    ofstream initoutput;
    string initoutputFile = outputDir + "/InitFile";
    initoutput.open(initoutputFile.c_str());
    initoutput <<"fileY: "<<fileY<<endl;
    initoutput <<"sep: "<<sep<<endl;
    initoutput <<"k: "<<k<<endl;
    initoutput <<"v: "<<v<<endl;
    for(int w=0; w<v; w++){
        stringstream ss;
        ss<<w+1;
        initoutput <<"p" + ss.str() + ": "<<pArray[w]<<endl;
    }
    initoutput <<"n: "<<n<<endl;
    initoutput <<"S: "<<S<<endl;
    initoutput <<"outputDir: "<<outputDir<<endl;
    initoutput <<"step: "<<step<<endl;
    initoutput <<"MCMCS: "<<MCMCS<<endl;
    initoutput <<"PXS: "<<PXS<<endl;
    initoutput <<"Seed: "<<myBASS.seed<<endl;
    initoutput.close();

    // define additional parameters
    // -------------------------------------------
    gsl_vector * LoadingN0num = gsl_vector_alloc(S - MCMCS - PXS);
    gsl_vector * LoadingFNorm = gsl_vector_alloc(S - MCMCS - PXS);
    gsl_vector * LogLikelihood = gsl_vector_alloc(S - MCMCS - PXS);

    // start iteration
    // ------------------------------------------
    // MCMC run
    myBASS.a = myBASS.b = 0.5;
    myBASS.c = myBASS.d = 0.5;
    myBASS.e = myBASS.f = 0.5;
    int collapse = 1;
    for(int s = 0; s<MCMCS; s++){
        if (s>MCMCS/2) collapse = 0;
        cout << s <<"th"<<" iteration of MCMC... "<<endl;
        myBASS.LambdaUpdateMCMC();
        myBASS.EtaUpdateMCMC();
        myBASS.HyperParamsUpdateMCMC(collapse);
        myBASS.Sigma2invUpdateMCMC();
    }
    //myBASS.Write2File(outputDir);
    //exit(0);
    // EM run
    int ss;
    myBASS.a = myBASS.b = 20;
    myBASS.c = myBASS.d = 10;
    myBASS.e = myBASS.f = 5;
    for(int s = MCMCS; s<MCMCS+PXS; s++){

        myBASS.a -= 1/PXS; myBASS.b -= 1/PXS;
        myBASS.c -= 1/PXS; myBASS.d -= 1/PXS;
        myBASS.e -= 1/PXS; myBASS.f -= 1/PXS;
        cout<<s<<"th"<<" iteration..."<<endl;
        //myBASS.DeleteEmptyFactor();
        myBASS.LambdaUpdateEM();
        myBASS.HyperParamsUpdateEM(0);
        myBASS.Sigima2invUpdateEM();

        myBASS.EtaUpdateEM(1);
    }
    if(PXS != 0){ myBASS.Rotate();}
    //myBASS.Write2File(outputDir);
    //exit(0);
    myBASS.a = myBASS.b = 0.5;
    myBASS.c = myBASS.d = 0.5;
    myBASS.e = myBASS.f = 0.5;

    for(int s = MCMCS+PXS; s<S; s++){
        cout << s <<"th"<<" iteration, ";
        cout << "k: "<< myBASS.k <<endl;
        myBASS.DeleteEmptyFactor();
        myBASS.LambdaUpdateEM();
        myBASS.HyperParamsUpdateEM(1);
        myBASS.Sigima2invUpdateEM();

        myBASS.EtaUpdateEM(0);

        //if(s%step == 0)
        //myBASS.Write2File(outputDir3);

        ss = s - MCMCS - PXS;
        gsl_vector_set(LoadingN0num, ss, myBASS.GetN0num());
        cout<<"non-zero loadings: " << gsl_vector_get(LoadingN0num,ss);
        gsl_vector_set(LoadingFNorm, ss, myBASS.GetFNorm());
        gsl_vector_set(LogLikelihood,ss,myBASS.GetLogLikelihood());
        cout<<"  logLikelihood: " << gsl_vector_get(LogLikelihood,ss)<<endl;
        if(ss >= step){

            int N0diff = gsl_vector_get(LoadingN0num,ss) -
                    gsl_vector_get(LoadingN0num,ss-step);
            double FNormdiff = abs(gsl_vector_get(LoadingFNorm,ss) -
                                   gsl_vector_get(LoadingFNorm,ss-step));
            double LogLikediff = abs(gsl_vector_get(LogLikelihood,ss) -
                                     gsl_vector_get(LogLikelihood,ss-step));
            

            if(ss % step == 0){ // output every step iterations
                stringstream strstr;
                strstr<<ss;
                string outputDirss = outputDir + "/iter" + strstr.str();
                int tag = mkdir(outputDirss.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                myBASS.DeleteEmptyFactor();
                myBASS.PruneFactors();
                //myBASS.WriteRotation(outputDirss);
                myBASS.Write2File(outputDirss);
            }
            if(N0diff == 0 && LogLikediff < 10e-5){  // converge
                myBASS.DeleteEmptyFactor();
                myBASS.PruneFactors();
                //myBASS.WriteRotation(outputDir);
                myBASS.Write2File(outputDir);
                break;
            }
        }
    }

    gsl_vector_free(LoadingN0num);
    gsl_vector_free(LoadingFNorm);
    gsl_vector_free(LogLikelihood);


}

