/*
  Bayesian group factor Analysis with Structured Sparsity (BASS)
  Copyright (C) 2015 Shiwen Zhao
  
  This source code is downloaded from 
  http://jonaswallin.github.io/articles/2013/07/simulation-of-gig-distribution/
  
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

 
#ifndef RNGSTREAM_H
#define RNGSTREAM_H
 
#include <string>

class RngStream
{
public:

RngStream (const char *name = "");


static bool SetPackageSeed (const unsigned long seed[6]);


void ResetStartStream ();


void ResetStartSubstream ();


void ResetNextSubstream ();


void SetAntithetic (bool a);


void IncreasedPrecis (bool incp);


bool SetSeed (const unsigned long seed[6]);


void AdvanceState (long e, long c);


void GetState (unsigned long seed[6]) const;


void WriteState () const;


void WriteStateFull () const;


double RandU01 ();


int RandInt (int i, int j);

unsigned long int RandUInt (unsigned long int i, unsigned long int j);

private:

double Cg[6], Bg[6], Ig[6];


bool anti, incPrec;


std::string name;


static double nextSeed[6];


double U01 ();


double U01d ();


};
 
#endif
 

