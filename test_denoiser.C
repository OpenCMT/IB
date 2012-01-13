#include "IB.h"

#include "debug.h"
DEBUG_HEADER

using std::string;

int main(int argc, char *argv[]) {

    /* WARNING: SV CODE DIFFERENCE IN GEOMETRY CONVENTION */
    /* X->phiview Z->thetaview Y->along distance from ch  */
    struct
    {
        int run_num;
        uint vxnum; // phi view
        uint vynum; // chamber distance
        uint vznum; // theta view
        uint kevent_num;
        float lambda_0;
        
    } parameters = {
        200,          // default run num is not set
        102,         // default phi view vox number
        53,         // default y vox numbers
        80,         // default theta view vox number
        1000,       // one milion events
        0.1E-6      // default lambda starting value
    };


/*  ------------------------ DEBUG INITIALIZATION --------------------------  */
    debug_init();
    IBMuonCollection_debug_init();
    IBVoxCollection_debug_init();
    IBAnalyzerEM_debug_init();

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


/*  -------------    READING SIMPLE COMMANDLINE PARAMETERS .................  */
    if (argc == 7) {
    	parameters.run_num = atoi(argv[1]);
    	parameters.vxnum = atoi(argv[2]);
    	parameters.vynum = atoi(argv[3]);
    	parameters.vznum = atoi(argv[4]);
    	parameters.kevent_num = atoi(argv[5]);
        parameters.lambda_0 = atof(argv[6]);
    }
    else {
        perror("using DEFAULT parameters");
    }
/*  -----------------------------------------------------------------------   */



    /* READING DATA */
    debug_message("reading file");
    TFile* hbfile = new TFile("/data_03/Radmufit_r829_0_18000000ev_PR.root");
    //forse errore nel costruttore (valgrind error)
    if (hbfile->IsZombie()) {
        debug_error("Error opening file");
        exit(-1);
    }
    debug_message("retriving RADMU TTree...");
    TTree* tree = (TTree*) hbfile->Get("RADMU");
    assert(tree);


    /* COLLECT MUONS DATA */
    IB_set_default_geometry(NULL);
    debug_message("reading muons from TTree");

    IBMuonCollection muons(
                           //comment here unwanted cuts
                           MCOLL_FLAG_SINGLECHCUT |
                           MCOLL_FLAG_NUMHITS |
                           //MCOLL_FLAG_POSCUT |
                           MCOLL_FLAG_SLOPECUT |
                           0);
    muons.addMuonFromTTree(tree,1000,parameters.kevent_num);
    
    printf("\rsono rimasti %d muoni ..\n",muons.size());

    /* SETUP VOXELS CONTAINER */
    IBInt3d dim(parameters.vxnum,
                parameters.vynum,
                parameters.vznum);

    IBVoxCollection voxels(dim);
    voxels.init_lambda(parameters.lambda_0);

    
    // set current domain to VoxCollection if you want to debug raytracer    //
    // raytracer dumps bounding_vtk and raytrace_vtk if current domain is    //
    // set to voxels.                                                        //
    //    debug_set_domain("IBVoxCollection");

    IBAnalyzerEM em(&muons, &voxels);

    printf("Performing cuts\n");
    em.run(1, 1);
    em.SijCut(2);
    em.run(10,1);
    em.chi2Cut(0,40);

    printf("Begin analyze\n");
    /* BEGIN ANALYZE */
    voxels.init_lambda(0);
    IBAnalyzerPoca pan(&muons,&voxels);
    pan.run(1,1);
    voxels *= 5;
    voxels.export_to_vtk("poca_analyzer_x5.vtk");

    IBInt3d fdim(4, 4, 4);
    IBVoxFilter3d_Spherical sphere_filter(fdim, IB_FILTER_SHAPE_GAUSS);
    IBVoxFilter3d_abtrim abtrim((IBVoxFilter3d &)sphere_filter);
    abtrim.atrim = 0;
    abtrim.btrim = 2;

    IBVoxCollection filtered = voxels;
    IBVoxCollection tmp = voxels;

    char out_file[50];
    int step = parameters.run_num;
    for (int index = 0; index < 200 / parameters.run_num; index++) {
        em.run(step, 1);
        sprintf(out_file, "voxel_it%d.vti\0", (index + 1) * step);
        voxels.export_to_vtk_xml(out_file);
        filtered = voxels;
        tmp = voxels;
        abtrim.run(filtered);
        sprintf(out_file, "filtered_it%d.vti\0", (index + 1) * step);
        filtered.export_to_vtk_xml(out_file);
        filtered *= 30;
        tmp -= filtered;
        sprintf(out_file, "vox-filt_it%d.vti\0", (index + 1) * step);
        tmp.export_to_vtk_xml(out_file);
        IBVoxThreshold threshold(tmp);
        threshold.run(filtered.mean());
        sprintf(out_file, "noise_it%d.vti\0", (index + 1) * step);
        tmp.export_to_vtk_xml(out_file);
        voxels -= tmp;
        sprintf(out_file, "denoised_it%d.vti\0", (index + 1) * step);
        voxels.export_to_vtk_xml(out_file);
    }

        


    
    delete (hbfile);
    return 0;
}

