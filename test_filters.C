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
        200,       //  events
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

    /* SETUP VOXELS CONTAINER */
    IBInt3d dim(parameters.vxnum,
                parameters.vynum,
                parameters.vznum);

    IBVoxCollection voxels(dim);
    voxels.init_lambda(parameters.lambda_0);


   

   IBInt3d fdim(5, 1, 5);
    float array[] = {
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
        1, 1, 1, 1, 1,
    };

    IBVoxFilter3d_Median median(fdim);
    median.set_xzy_numeric(voxels, array);


    fdim.set(5,5,5);
    IBVoxFilter3d_Spherical sphere_filter(fdim, IB_FILTER_SHAPE_GAUSS);
    IBVoxFilter3d_abtrim abtrim((IBVoxFilter3d &)sphere_filter);
    abtrim.atrim = 0;
    abtrim.btrim = 3;
//    abtrim.set_xzy_numeric(voxels, array);

    fdim.set(5,5,5);
    IBVoxFilter3d_Spherical sphere_filter5(fdim, IB_FILTER_SHAPE_GAUSS);
    // standard main ...............


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
    muons.addMuonFromTTree(tree,parameters.kevent_num);
    
    printf("\rsono rimasti %d muoni ..\n",muons.size());


    
    // set current domain to VoxCollection if you want to debug raytracer    //
    // raytracer dumps bounding_vtk and raytrace_vtk if current domain is    //
    // set to voxels.                                                        //
    //    debug_set_domain("IBVoxCollection");

    IBAnalyzerEM em(&muons, &voxels);



    printf("Performing cuts\n");

    em.run(1, 1);
    em.SijCut(2);
    printf("[%d] Muons survived\n",em.size());

    em.run(10,1);

    em.chi2Cut(0,30);
    printf("[%d] Muons survived\n",em.size());



    printf("Begin analyze\n");
    /* BEGIN ANALYZE */
    voxels.init_lambda(parameters.lambda_0);


////////////////////////////////////////////////////////////////////////////////
////////////////// POCA           //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    IBAnalyzerPoca pan(&muons,&voxels);
    pan.run(1,1);
    voxels *= 25;
    voxels.export_to_vtk("poca_analyzer_x5.vtk");


////////////////////////////////////////////////////////////////////////////////
////////////////// EM SIMPLE      //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//    voxels.init_lambda(parameters.lambda_0);
//    IBAnalyzerEM_simple ems(em);
    //ems.run(1,1);
    //voxels.export_to_vtk("ems1.vtk");




////////////////////////////////////////////////////////////////////////////////
////////////////// EM STANDARD    //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


//    voxels.init_lambda(parameters.lambda_0);

    IBVoxCollection filtered = voxels;
    char file_name[50];
 //   for (int i = 0; i < parameters.run_num; i++) {
        em.run(200, 1);

        sprintf(file_name,"poca_em_%dev_it%d.vtk",em.size(),400);
        voxels.export_to_vtk(file_name);

        filtered = voxels;
        abtrim.run(filtered);
        sprintf(file_name,"poca_em_%dev_abtrim_it%d.vtk",em.size(),400);
        filtered.export_to_vtk(file_name);

        sphere_filter5.run(filtered);
        sprintf(file_name,"poca_em_%dev_filtered_it%d.vtk",em.size(),400);
        filtered.export_to_vtk(file_name);
//        IBVoxThreshold th(filtered);
//        printf("mean -> %f -> filtered = %f\n", voxels.mean()*1E6, filtered.mean()*1E6);
//        printf("count -> th.count(%f) = %d\n",10.,th.count(10E-6));
    



    delete (hbfile);
    return 0;
}

