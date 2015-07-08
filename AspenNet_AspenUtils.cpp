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

// Note: this is based heavily off of the "newruntime" tool in Aspen's source tree

#include <iostream>
#include <deque>
#include <cstdio>
#include <map>

#include "model/ASTAppModel.h"
#include "model/ASTMachModel.h"
#include "parser/Parser.h"
#include "walkers/RuntimeCounter.h"
#include "walkers/RuntimeExpression.h"

#include <sys/time.h>

extern "C"{
double runtimeCalc(char *a, char *m, char * socket);
int getSockets(char *m, char*** buf);


double runtimeCalc(char *a, char *m, char * socket)
{
    if (a && m && socket){
        
        ASTAppModel *app = LoadAppModel(a);
        ASTMachModel *mach = LoadMachineModel(m);
        
        struct timeval tv;
        gettimeofday(&tv, NULL);
        srand(tv.tv_usec);
        
        RuntimeCounter *t = new RuntimeCounter(app, mach, socket);
        t->SetCacheExecutionBlockExpressions(false);
        app->kernelMap["main"]->Traverse(t);
        fprintf(stderr, "INFO: Returning Aspen result of %f.\n", t->GetResult());       
        return t->GetResult();
        
    }
    else {
        
        return -1;
        // Will need to check for this error condition when \
        receiving the value in the C code
        
    }
}

int getSockets(char *m, char*** buf){
    if (m){

        ASTMachModel *mach = LoadMachineModel(m);
        fprintf(stderr,"INFO: Planning to return %d sockets.\n",  mach->socketlist.size());
        *buf = (char**) malloc(sizeof(char*) * mach->socketlist.size());
        
        // Copy the strings over as char*'s
        for (int i = 0; i < mach->socketlist.size(); i++){
            (*buf)[i] = (char*) malloc(sizeof(char) * 30);
            for (int j = 0; j < 30; j++){
                if (j < mach->socketlist[i].size()){
                    (*buf)[i][j] = mach->socketlist[i][j];
                }
                else {
                    (*buf)[i][j] = '\0';
                }
            }
            
        }
        return mach->socketlist.size();
    }
    return -1;
}
}
