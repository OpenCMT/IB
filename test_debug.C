#include "IB.h"
#include "stdio.h"

#include "debug.h"
DEBUG_HEADER



DEBUG_INIT(demo) {
    DEBUG_LOCAL_REF = debug_set_domain("demo");
}


int main(int argc, char *argv[]) {

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

    P("entering message")

//
//    debug_handle("triggered_message","Message var:%d",123); // create trigger
//    debug_message("\nThis is a triggered message:");
//    debug_handle("triggered_message",NULL); // shot
    debug_message("hello: this is debug test");
    
//    demo_debug_init();
//    debug_set_active_domain("demo");
//    printf("domain demo active flag = %d\n",
//           debug_test_active(DEBUG_LOCAL_REF, 5));


    IBVoxCollection a(10,10,10);

    return 0;
}



