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
#include "aspenc.h"

int main(){
    int i;
    char **buf = NULL;
    int size = -1;
    char* app = "./models/matmul/matmul.aspen";
    char* machine = "./models/machine/simple.aspen";
    
    AppModel_p app_model = Aspen_LoadAppModel(app);
    printf("Loaded model: %s\n", AppModel_GetName(app_model));
    
    Kernel_p kernel = AppModel_GetMainKernel(app_model);
    printf("Main kernel is: %s\n", Kernel_GetName(kernel));

    MachModel_p mach_model = Aspen_LoadMachModel(machine);
    MachComponent_p mach = MachModel_GetMachine(mach_model);
    printf("Machine Model: %s\n", MachComponent_GetName(mach));

    Expression_p runtime_expr = Kernel_GetTimeExpression(kernel, app_model, mach_model, "SimpleCPU");
    ParamMap_p appParams = AppModel_GetParamMap(app_model);
    ParamMap_p machParams = MachModel_GetParamMap(mach_model);

    Expression_p rt_exp1 = Expression_Expanded(runtime_expr, ParamMap_Create("n", 277));
    // Above: "Runtime expanded by n = 277"
    Expression_p rt_exp2 = Expression_Expanded(Expression_Expanded(rt_exp1, appParams), machParams);

    printf("Expression as value: %lf\n", Expression_Evaluate(rt_exp2));

    return 0;
}
