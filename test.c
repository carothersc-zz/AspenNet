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

// Prototypes for extern C functions:
double runtimeCalc(char *a, char *m, char * socket);
int getSockets(char *m, char*** buf);


int main(){
    char **buf = NULL;
    int size = -1;

    printf("%e\n", runtimeCalc("./models/fft/1D_FFT.aspen", "./models/machine/TestRig.aspen", "amd_HD5770"));

    size = getSockets("./models/machine/TestRig.aspen", &buf);

    if (size == -1){
        fprintf(stderr, "ERROR: No int was stored!\n");
    }
    else {
        printf("%d sockets should be returned.\n", size);
        if (!buf){
            fprintf(stderr, "ERROR: No sockets were returned!\n");
        }
        else {
            for (int i = 0; i < size; i++){
                printf("%s\n", buf[i]);
            }
        }   
    }
    return 0;
}
