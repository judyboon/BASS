#!/bin/bash                                                                                                                                       
# an example of running BASS on sim5 with n = 100

k=15;  # initial k
v=10;  # number of observations

n=100  # number of samples
S=20000  # maximum iteration
sep="space"  # separation in fileY
step=50  # step size to record results
MCMCS=0  # MCMC iteration
PXS=20  # PX-EM iteration

fileY="../simData/Yn100sim5.sim";

# following to generate command to run BASS
string="./main -y "${fileY}" -sep "${sep}
string+=" -k "$k" -v "$v

for (( i=1; i<=$v; i++ ))
do
    string+=" -p"$i" 50"
done

string+=" -n "$n" -itr "$S" -step "$step" -mcmc "$MCMCS" -px "$PXS

echo $string
fileOut="../output/MCMC${MCMCS}PX${PXS}ik${k}n${n}"
string=$string" -out "$fileOut
eval $string

