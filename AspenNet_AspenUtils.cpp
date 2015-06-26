/*
 * The purpose of this file is to provide cross-compilable
 * aspen-based utils to interface with the C code in the 
 * CODES_AspenNet. In this way the CODES model can use runtime estimates.
 * Author: Mark Blanco
 * Date: 26 June 2015
 * Credit to Jeff Vetter and Jeremy Meredith (Oakridge National Labs) for Aspen code
 * CODES can be found at http://www.mcs.anl.gov/research/projects/codes/
 * (Credit to Argonne National Laboratory and Rensselaer Polytechnic Institute)
 */

#include <iostream>
#include <deque>
#include <cstdio>
#include <map>

#include "model/ASTAppModel.h"
#include "model/ASTMachModel.h"
#include "parser/Parser.h"
#include "walkers/RuntimeCounter.h"
#include "walkers/RuntimeExpression.h"
// Going to need to figure out how to get most of these includes to work...