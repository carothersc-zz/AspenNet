#include <string.h>
#include <assert.h>
#include <ross.h>

#include "codes/lp-io.h"
#include "codes/codes.h"
#include "codes/codes_mapping.h"
#include "codes/configuration.h"
#include "codes/model-net.h"
#include "codes/lp-type-lookup.h"
#include "codes/local-storage-model.h"

double runtimeCalc(char *a, char *m, char * socket);
void getSockets(char *m, char*** buf, int * size);



int main(){
    printf("%f", runtimeCalc("models/fft/1D_FFT.aspen", "models/machine/TestRig.aspen", "amd_HD5770"));
    return 0;
}
