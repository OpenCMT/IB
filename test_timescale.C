
/* TEST PER LA MISURA DELLO SCALING DI TEMPO NELLA DIMENSIONE DEI VOXELS
 *
 * la macro TPROBES definisce il numero di tentativi da mediare
 * per ottenere la stima del tempo di calcolo per ogni passo del
 * processo. questo test va invocato con il comando:
 * ./test_timescale Nkev dimv1 dimv2 dimv3 dimv4 ...
 *
 *  OUTPUT: produce un output in timing_results.csv
 * le colonne rappresentano i tempi e gli errori per:
        size del voxel,
        root_loader,root_loader_err,
        EM_constructor,EM_constructor_err,
        Sijcut,Sijcut_err,
        Chi2cut,Chi2cut_err,
        poca_iter,poca_iter_err,
        iterations,iterations_err,
        filesave,filesave_err,
 *
 * NOTA BENE !
 * in questo test non sono veramente eseguiti i cuts in modo da non modificare
 * lo scaling delle iterazione al variare della dimensione dei voxels
 *
 */





#include "IB.h"
#include <time.h>
#include <omp.h>
#include <math.h>
#include "debug.h"
DEBUG_HEADER

using std::string;

int main(int argc, char *argv[]) {

    struct    {
        uint kevent_num;
        float Sij;
        float Chi2;
        vector<float> vox_size;
    } parameters = {
        500, // ~ 250k events
        2.,
        30.
    };


/*  ------------------------ TIMING FUNCTIONS AND MACROS -------------------  */

    struct  {
        unsigned int muons_aq, muons_left_Sij, muons_left_Chi2;
        double root_loader,root_loader_err;
        double EM_constructor,EM_constructor_err;
        double Sijcut,Sijcut_err;
        double Chi2cut,Chi2cut_err;
        double poca_iter,poca_iter_err;
        double iterations,iterations_err;
        double filesave,filesave_err;
    } timing_results = {
        0, 0, 0,
        0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0
    };

    ////// NUMBER OF PROBES
#define TPROBES 10
    double tstart, tend, timing_probes[TPROBES];

    ////// TSTART MACRO
#define TSTART \
for(int _timing_id = 0 ; _timing_id < TPROBES ; _timing_id++ ) {\
tstart = omp_get_wtime();

    ////// TEND MACRO
#define TEND \
tend = omp_get_wtime(); \
timing_probes[_timing_id] = (double)(tend-tstart); \
}

    ////// TSTORE MACRO
#define TSTORE(item) \
for(int _timing_id = 0 ; _timing_id < TPROBES ; _timing_id++ ) {\
timing_results.item += timing_probes[_timing_id];\
} \
timing_results.item /= TPROBES; \
for(int _timing_id = 0 ; _timing_id < TPROBES ; _timing_id++ ) {\
timing_results.item##_err += \
(timing_probes[_timing_id] - timing_results.item) * \
(timing_probes[_timing_id] - timing_results.item);  \
} \
timing_results.item##_err /= (TPROBES - 1);\
timing_results.item##_err = sqrt(timing_results.item##_err);

    
    ////// OUTPUT FILE /////////////////////////////////////////
    FILE *tf = fopen("timing_results.csv", "w");




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
    if (argc >= 2) {
    	parameters.kevent_num = atoi(argv[1]);
        for (int i = 2; i < argc; i++)
            parameters.vox_size.push_back(atof(argv[i]));
    }
    else {
        printf("Please usage ./test_timescale k_ev_number sizes ... \n");
        exit(1);
    }
/*  -----------------------------------------------------------------------   */


    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  ACQUIRE DATA  /////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

/* ------ DUMMY READING DATA ----------------------------------------------   */
    TSTART
    IBMuonCollection muons(
                           //comment here unwanted cuts
                           MCOLL_FLAG_SINGLECHCUT |
                           MCOLL_FLAG_NUMHITS |
                           //MCOLL_FLAG_POSCUT |
                           MCOLL_FLAG_SLOPECUT |
                           0);

    ////////// READING OF ROOT FILE ////////////////////////////////////////////
    TFile * hbfile;
    printf("Reading ROOT file\n");
    hbfile = new TFile("/data_03/Radmufit_r829_0_18000000ev_PR.root");
    //forse errore nel costruttore (valgrind error)
    if (hbfile->IsZombie()) {
        debug_error("Error opening file");
        exit(-1);
    }
    TTree* tree = (TTree*) hbfile->Get("RADMU");
    assert(tree);

    /////////// COLLECT MUONS DATA /////////////////////////////////////////////
    IB_set_default_geometry(NULL);
    muons.addMuonFromTTree(tree, 0, parameters.kevent_num);
    delete (hbfile);
    TEND TSTORE(root_loader);

    
/* -------- REAL DATA LOADER   --------------------------------------------   */
    IBMuonCollection muons(
                           //comment here unwanted cuts
                           MCOLL_FLAG_SINGLECHCUT |
                           MCOLL_FLAG_NUMHITS |
                           //MCOLL_FLAG_POSCUT |
                           MCOLL_FLAG_SLOPECUT |
                           0);
    TFile * hbfile;
    printf("Reading ROOT file\n");
    hbfile = new TFile("/data_03/Radmufit_r829_0_18000000ev_PR.root");
    //forse errore nel costruttore (valgrind error)
    if (hbfile->IsZombie()) {
        debug_error("Error opening file");
        exit(-1);
    }
    TTree* tree = (TTree*) hbfile->Get("RADMU");
    assert(tree);

    /////////// COLLECT MUONS DATA /////////////////////////////////////////////
    IB_set_default_geometry(NULL);
    muons.addMuonFromTTree(tree, 0, parameters.kevent_num);
    delete (hbfile);
/* ------------------------------------------------------------------------   */






    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  LOOP FOR VOXELS SCALING  //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    while (!parameters.vox_size.empty()) {

        float size = parameters.vox_size.back();
        parameters.vox_size.pop_back();

        uint vxnum, vynum, vznum;
        vxnum = IB_default_geometry->main_roi_size.x / size;
        vynum = IB_default_geometry->main_roi_size.y / size;
        vznum = IB_default_geometry->main_roi_size.z / size;
        /* SETUP VOXELS CONTAINER */
        IBInt3d dim(vxnum, vynum, vznum);
        printf("////////////////////////////////////////////////////\n"
               "size = %f -> Building VoxCollection [%d %d %d ]\n", size,
               vxnum, vynum, vznum);
        IBVoxCollection voxels(dim);
        voxels.init_lambda(1.E-6);

        //////// EM CONSTRUCTOR ////////////////////////////
        /// this execute TPROBE TIMES a dummy constructor //
        TSTART
        IBAnalyzerEM em(&muons, &voxels);
        TEND TSTORE(EM_constructor);

        // this is real constructor //
        IBAnalyzerEM em(&muons, &voxels);


        ////// CUTS ///////////////////////////////////////
        TSTART
        printf("Sijcut\n");
        em.run(1, 1);
        //em.SijCut(parameters.Sij);
        TEND TSTORE(Sijcut);

        TSTART
        printf("Chi2cut\n");
        em.run(10,1);
        //em.chi2Cut(0,parameters.Chi2);
        TEND TSTORE(Chi2cut)


        ////// POCA ////////////////////////////////////////
        voxels.init_lambda(0);
        TSTART
        IBAnalyzerPoca pan(&muons, &voxels);
        pan.run(1, 1);
        TEND TSTORE(poca_iter);
        voxels *= 5;

        ////// EM //////////////////////////////////////////
        TSTART
        em.run(200, 1);
        TEND TSTORE(iterations);


        ////// SAVE ////////////////////////////////////////
        TSTART
        voxels.export_to_vtk("poca_em200.vtk");
        TEND TSTORE(filesave);


        ////// DUMP TO CSV ////////////////////////////////
        fprintf(tf,
                "%.3f,\t"
                "%.3f,\t%.3f,\t%.3f,\t%.3f,\t%.3f,\t%.3f,\t%.3f\t"
                "%.3f,\t%.3f,\t%.3f,\t%.3f,\t%.3f,\t%.3f,\t%.3f\n",
                size,
                timing_results.root_loader,
                timing_results.root_loader_err,
                timing_results.EM_constructor,
                timing_results.EM_constructor_err,
                timing_results.Sijcut,
                timing_results.Sijcut_err,
                timing_results.Chi2cut,
                timing_results.Chi2cut_err,
                timing_results.poca_iter,
                timing_results.poca_iter_err,
                timing_results.iterations,
                timing_results.iterations_err,
                timing_results.filesave,
                timing_results.filesave_err
        );
        fflush(tf);
    }
    fclose(tf);

    return 0;
}

