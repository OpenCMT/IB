#include "IB.h"
#include "debug.h"

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
        0,          // default run num is not set
        44,         // default phi view vox number
        23,         // default y vox numbers
        34,         // default theta view vox number
        100,       // K events
        0.1E-6      // default lambda starting value
    };


/*  ------------------------ DEBUG INITIALIZATION --------------------------  */
    debug_init();
    IBMuonCollection_debug_init();
    IBVoxCollection_debug_init();
    IBAnalyzerEM_debug_init();

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


    if (argc == 2) {
        parameters.run_num = atoi(argv[1]);
    }

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

    muons.addMuonFromTTree(tree,parameters.kevent_num);



    /* SETUP VOXELS CONTAINER */
    IBVoxCollection voxels(parameters.vxnum,
                           parameters.vynum,
                           parameters.vznum);
    voxels.init_lambda(parameters.lambda_0);

    debug_set_domain("IBMuonCollection");
    debug_current("count_poca","Poca count",&muons,&voxels);

    int start = parameters.run_num;
    int end = start + 1;
    debug_set_domain("IBMuonCollection");
    vector <IBMuon *> v(muons.muons().begin()+start,muons.muons().begin()+end);
    debug_current("muons_vtk","muon_plot.vtk",&muons,&v);
    debug_current("poca_vtk", "poca_plot.vtk", &v);
    debug_current("ch_position_vtk", "muon_ch_position.vtk", 
		  &muons, &muons.muons());

    debug_set_domain("IBVoxCollection");
    debug_current("container_vtk","container_box.vtk",&voxels);
    debug_current("muon_entry_vtk","muon_entry.vtk",&voxels,&muons,&v);
    debug_current("muon_entry_vtk","muons_entry.vtk",
		  &voxels,&muons,&muons.muons());
 //   for (int i = start; i < end ; i++)
    IBRay *ray = voxels.raytrace(muons.mu(start));
    printf("raytrace length = %f\n",ray->length);
//    for(int i=0;i<ray->L->size();i++)
//        printf("%d -> Lij = %f\n", i, ray->L->at(i));


    debug_set_domain("default");
    debug_message("Creating Analyzer: muon -> EMmuon");
    IBAnalyzerEM em(&muons,&voxels);

    voxels.init_lambda(parameters.lambda_0);
    em.print_ray_parameters(0);


 
    delete (hbfile);
    return 0;
}

