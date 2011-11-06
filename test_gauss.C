#include "IB.h"

#include "debug.h"
DEBUG_HEADER

using std::string;

int main(int argc, char *argv[]) {

/*  ------------------------ DEBUG INITIALIZATION --------------------------  */
    debug_init();

    // set default as current domain //
    debug_set_domain("default");
    char *log_domain_environment = getenv(LOG_ENVIRONMENT_DOMAIN);
    if (log_domain_environment)
        debug_set_active_domain(log_domain_environment);
    else
        debug_set_active_domain(LOG_DEFAULT_DOMAIN_NAME);

    log_domain_environment = getenv(LOG_ENVIRONMENT_OUTPUT);
    if (log_domain_environment)
        debug_set_active_output(log_domain_environment);
    else
        debug_set_active_output("stdout");
/*  ------------------------------------------------------------------------  */

    IB_set_default_geometry(NULL);

    IBInt3d dim(100,100,100);
    IBVoxCollection voxels(dim);
    voxels.init_lambda(0);

    voxels.get(IBInt3d(50,50,50))->density = 1.;

    IBInt3d fdim(5,5,5);
    IBVoxFilter3d_Spherical sphere_filter(fdim, IB_FILTER_SHAPE_GAUSS);

    printf("voxels mean before filter = %f\n",voxels.mean());
    sphere_filter.run(voxels);
    printf("voxels mean after filter = %f\n",voxels.mean());
    voxels.export_to_vtk("filter_gauss_test.vtk");

    return 0;
}

