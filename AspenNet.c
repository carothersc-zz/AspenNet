/*
 * Aspen-Codes integration model
 * Author: Mark Blanco
 * Date: 24 June 2015
 * Note: this model is directly based off of the codes-base example (example.c)
 * They may look very similar for a while
 * The CODES project can be found at http://www.mcs.anl.gov/research/projects/codes/
 * The ROSS project, which codes is based on, can be found at 
 * https://github.com/carothersc/ROSS/wiki
 * Finally, ASPEN can be found at https://ft.ornl.gov/node/1566... sort of.
 */

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

#include "AspenNet.h"

/* Global variable for the aspen config file (ROSS/network config handle is
 * a global preallocated/declared by configuration.c */
ConfigHandle aspen_config;
/* NOTE: the aspen config file should contain the paths to the aspen kernel model(s)
 * and hardware model, as well as the socket on which the kernel model should be
 * evaluated. */

/* Global for path and name of Aspen model */
char Aspen_Mach_Path[100];
/* Global array for paths and names of Aspen application/kernels */
char **Aspen_App_Path = NULL;
/* Global array for names of sockets to be used */
char **Aspen_Socket = NULL;
// TODO: Do something better than having a hard-coded length...
/* Global variables to specify how destination servers should be found */
char network_traffic_type[100];
static int traffic_pattern_number = 0;
/* number of requests sent by each server (read from config) */
static int num_reqs = 0;
/* size of simulated data payload, bytes (read from config) */
static int payload_sz = 0;
/* number of computation-simulation rounds to perform (read from config) */
static int num_rounds = 0;
/* flag for debug output. (read from config) */
static int debug_output = 0;

/* model-net ID, can be either simple-net, dragonfly or torus (more will be added) */
static int net_id = 0;
/* Global for total number of Aspen LPs in the simulation */
static int num_servers = 0;
static int offset = 2;
/* Global for total number of LPs in the simulation */
static int ttl_lps = 0;

/* Expected key for ROSS/network config file */
static char *network_conf_key = "network_conf_file";
/* Expected LP group name in configure files for this program */
static char *group_name = "ASPEN_SERVERS";
/* Expected parameter group name for rounds of communication */
static char *param_group_nm = "server_pings";
/* Expected misc parameters group name */
static char *misc_param_gp_nm = "PARAMS";
/* Expected name of conf group for Aspen file path parameters */
static char *aspen_group_nm = "ASPEN_PARAMS";
/* Number of network requests to be sent by each LP per round */
static char *num_reqs_key = "num_reqs";
/* The size of each message's payload (in KB) */
static char *payload_sz_key = "payload_sz";
/* Template names for aspen file path keys and socket choices: */
static char aspen_app_key[] = "aspen_app_path000";
static char *aspen_mach_key = "aspen_mach_path";
static char aspen_socket_key[] = "socket_choice000";
/* The number of network-computation rounds to be performed in the
 * simulation (set from config file) */
static char *num_rounds_key = "num_rounds";
/* Static char for name of traffic type parameter: */
static char *traffic_type_key = "network_traffic_pattern";
static char *aspen_svr_name = "server";

/* Global value to keep track of total runtime */
tw_stime totalRuntime = 0;
/* Global value to specifically keep track of ASPENCOMP rollbacks */
unsigned int computationRollbacks = 0;
/* Global value to keep track of the # of network-computation rounds performed */
unsigned int roundsExecuted = 0;

/* Prototypes for two extern C functions:
 *  NOTE: in order for the functions (and the whole simulation)
 *  to function, this file must be linked with AspenNet_AspenUtils.cpp,
 *  which contains the function definitions. AspenNet_AspenUtils.cpp
 *  #includes several Aspen source files from which its functionality
 *  is derived. */
extern double runtimeCalc(char *a, char *m, char * socket);
extern int getSockets(char *m, char*** buf);

/* Main function:
 * handles getting configuration options from config file and setup of ROSS/CODES
 * simulation environment (including LPs) */
int main(int argc, char *argv[])
{
    int nprocs;
    int rank;
    int num_nets, *net_ids;
    
    g_tw_lookahead = 0.5;
    
    g_tw_ts_end = s_to_ns(60*60*24*365); /* one year, in nsecs */
    
    /* ROSS initialization function calls */
    tw_opt_add(app_opt); /* add user-defined args */
    /* initialize ROSS and parse args. NOTE: tw_init calls MPI_Init */
    tw_init(&argc, &argv);
    if ( !aspen_conf_file_name[0] )
    {
        fprintf(stderr, "Expected \"--conf\" option with AspenNet configuration file.\n");
        MPI_Finalize();
        return 1;
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    /* loading the aspen config file into the codes-mapping utility, giving us the
     * parsed aspen_config object in return. */
    if ( configuration_load(aspen_conf_file_name, MPI_COMM_WORLD, &aspen_config) )
    {
        fprintf(stderr, "Error loading AspenNet config file %s.\n", aspen_conf_file_name);
        MPI_Finalize();
        return 1;
    }
    /* Retrieve the name of the network config file, stored in the aspen config file */
    configuration_get_value(&aspen_config, aspen_group_nm, network_conf_key, NULL,\
                            network_conf_file_name, 256);
    /* "config" is a global var defined by codes-mapping (which is now only used
     * for network-specific settings) */
    if ( configuration_load(network_conf_file_name, MPI_COMM_WORLD, &config) )
    {
        fprintf(stderr, "Error loading network config file %s.\n", network_conf_file_name);
        MPI_Finalize();
        return 1;
    }
    /* retreive debug_output flag from conf file */
    configuration_get_value_int(&aspen_config, aspen_group_nm, "debug_output", NULL, \
                                &debug_output);
    if (debug_output && !g_tw_mynode)
    {
        printf("\n***** DEBUG LEVEL %d: WILL PRINT OUTPUT *****\n", debug_output);
    }
    /* register model-net LPs with ROSS */
    model_net_register();
    
    /* register the server LP type with ROSS */
    svr_add_lp_type();
    
    /* Setup takes the global config object, the registered LPs, and
     * generates/places the LPs as specified in the configuration file.
     * This should only be called after ALL LP types have been registered in
     * codes */
    codes_mapping_setup();
    
    /* Setup the model-net parameters specified in the global config object,
     * returned are the identifier(s) for the network type. In this example, we
     * only expect one*/
    net_ids = model_net_configure(&num_nets);
    assert(num_nets==1);
    net_id = *net_ids;
    free(net_ids);
    
    /* calculate the number of servers in this simulation, and total LPs
     * ignoring annotations */
    num_servers = codes_mapping_get_lp_count(group_name, 0, aspen_svr_name, NULL, 1);
    MPI_Allreduce(&g_tw_nlp, &ttl_lps, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    
    /* for this example, we read from a separate configuration group for
     * server message parameters. Since they are constant for all LPs,
     * go ahead and read them prior to running */
    configuration_get_value_int(&aspen_config, param_group_nm, num_reqs_key, NULL, &num_reqs);
    configuration_get_value_int(&aspen_config, param_group_nm, payload_sz_key, NULL, &payload_sz);
    /* Now, if this is the 0th MPI rank, read the Aspen configuration paths and filenames,
     * as well as the total number of simulation rounds to run. */
    configuration_get_value(&aspen_config, aspen_group_nm, traffic_type_key, NULL, \
                            &network_traffic_type, 100);
    if (g_tw_mynode == 0)
    {
        int i, j, temp;
        /* Get the number of rounds: */
        configuration_get_value_int(&aspen_config, aspen_group_nm, num_rounds_key, \
                                    NULL, &num_rounds);
        
        Aspen_App_Path = calloc(num_rounds+1, sizeof(char*));
        Aspen_Socket = calloc(num_rounds+1, sizeof(char*));
        
        if (num_rounds > 1)
        {
            char **buffer = malloc(sizeof(char*));
            for (i = 0; i < num_rounds; i++)
            {
                // Load the app path for each round:
                temp = int_to_array(i, buffer);
                for (j = 1; j <= temp; j++)
                {
                    aspen_app_key[17 - j] = (*buffer)[temp - j];
                }
                free(*buffer);
                Aspen_App_Path[i] = malloc(100 * sizeof(char));
                configuration_get_value(&aspen_config, aspen_group_nm, aspen_app_key, \
                                        NULL, Aspen_App_Path[i], 100);
                
                // Load the socket choice for each round:
                temp = int_to_array(i, buffer);
                for (j = 1; j <= temp; j++)
                {
                    aspen_socket_key[16 - j] = (*buffer)[temp - j];
                }
                free(*buffer);
                Aspen_Socket[i] = malloc(100 * sizeof(char*));
                configuration_get_value(&aspen_config, aspen_group_nm, aspen_socket_key, \
                                        NULL, Aspen_Socket[i], 100);
            }
            free(buffer);
        }
        else
        {
            Aspen_App_Path[0] = calloc(100, sizeof(char));
            configuration_get_value(&aspen_config, aspen_group_nm, aspen_app_key, NULL, \
                                    Aspen_App_Path[0], 100);
            Aspen_Socket[0] = malloc(100 * sizeof(char));
            configuration_get_value(&aspen_config, aspen_group_nm, aspen_socket_key, NULL, \
                                    Aspen_Socket[i], 100);
        }
        /* Finally, load the aspen machine path, which should be the same for all rounds. */
        configuration_get_value(&aspen_config, aspen_group_nm, aspen_mach_key, NULL, \
                                &Aspen_Mach_Path, 100);
        // TODO: remove hard-coded length of 100!
        if (debug_output)
        {
            if (network_traffic_type && network_traffic_type[0] != '\0')
            {
                printf("INFO: Loaded network traffic type %s.\n", network_traffic_type);
            }
            else
            {
                fprintf(stderr, "WARNING: No network traffic pattern loaded from conf file!\n"\
                        "Defaulting to nearest-neighbor.\n");
                memcpy(network_traffic_type, "nearest-neighbor", 16);
            }
            printf("INFO: Will execute %d network-computation rounds.\n", num_rounds);
            for (i = 0; i < num_rounds; i++)
            {
                printf("\tAspen app model path loaded: %s\n"\
                       "\tAspen socket choice loaded: %s\n",\
                       Aspen_App_Path[i], Aspen_Socket[i]);
            }
            printf("INFO: Aspen machine model path loaded: %s\n", Aspen_Mach_Path);
        }
    }
    
    /* Set traffic pattern number (global): */
    if (network_traffic_type && network_traffic_type[0] != '\0')
    {
        if (!strcmp(network_traffic_type, "nearest-neighbor"))
        {
            traffic_pattern_number = 0;
        }
        else if (!strcmp(network_traffic_type, "random"))
        {
            traffic_pattern_number = 1;
        }
        else
        {
            traffic_pattern_number = 0;
            fprintf(stderr, "WARNING: Invalid traffic pattern entered! Defaulting to nearest-neighbor.\n");
        }
    }
    else
    {
        traffic_pattern_number = 0;
        fprintf(stderr, "WARNING: Invalid traffic pattern entered! Defaulting to nearest-neighbor.\n");
    }
    
    /* begin simulation */
    tw_run();
    
    if (g_tw_mynode == 0)
    {
        printf("\nFINAL REPORT: The final runtime for the application "\
               "kernel(s) is %f seconds.\n", totalRuntime);
        if (debug_output)
        {
            printf("INFO: Aspen computation was rolled back %u times.\n", computationRollbacks);
            printf("INFO: Aspen computation was performed %u times.\n", roundsExecuted);
        }
        /* Sanity check for optimized scheduler runs: */
        if (roundsExecuted != num_rounds + computationRollbacks)
        {
            fprintf(stderr, "ERROR: Aspen computation was performed an incorrect number of times!\n");
        }
    }
    tw_end();
    /* Memory cleanup on the 0th MPI rank: */
    if (g_tw_mynode == 0)
    {
        int i;
        for (i = 0; i < num_rounds; i++)
        {
            free(Aspen_App_Path[i]);
            free(Aspen_Socket[i]);
        }
        free(Aspen_App_Path);
        free(Aspen_Socket);
    }
    return 0;
}

const tw_lptype* svr_get_lp_type()
{
    return(&svr_lp);
}

static void svr_add_lp_type()
{
    /* lp_type_register should be called exactly once per process per
     * LP type */
    lp_type_register(aspen_svr_name, svr_get_lp_type());
}

static void aspen_svr_init(
                           aspen_svr_state * ns,
                           tw_lp * lp)
{
    tw_event *e;
    aspen_svr_msg *m;
    tw_stime kickoff_time;
    
    memset(ns, 0, sizeof(*ns));
    ns->end_ts = 0; // Set this to 0 in order to use it as a flag later
    /* each server sends a dummy event to itself that will kick off the real
     * simulation
     */
    
    /* skew each kickoff event slightly to help avoid event ties later on */
    kickoff_time = g_tw_lookahead + tw_rand_unif(lp->rng);
    
    /* first create the event (time arg is an offset, not absolute time) */
    e = codes_event_new(lp->gid, kickoff_time, lp);
    /* after event is created, grab the allocated message and set msg-specific
     * data */
    m = tw_event_data(e);
    m->aspen_svr_event_type = KICKOFF;
    m->src = lp->id;
    /* event is ready to be processed, send it off */
    tw_event_send(e);
    
    return;
}

/* event processing entry point
 * - simply forward the message to the appropriate handler */
static void aspen_svr_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    switch (m->aspen_svr_event_type)
    {
        case REQ:
            handle_req_event(ns, b, m, lp);
            break;
        case ACK:
            handle_ack_event(ns, b, m, lp);
            break;
        case KICKOFF:
            handle_kickoff_event(ns, b, m, lp);
            break;
        case LOCAL:
            handle_local_event(ns, b, m, lp);
            break;
        case DATA:
            handle_data_event(ns, b, m, lp);
            break;
        case ASPENCOMP:
            handle_computation_event(ns, b, m, lp);
            break;
        case RESTART:
            handle_restart_event(ns, b, m, lp);
            break;
        default:
            fprintf(stderr, "\n Invalid message type %d ", m->aspen_svr_event_type);
            assert(0);
            break;
    }
}

/* reverse event processing entry point
 * - simply forward the message to the appropriate handler */
static void aspen_svr_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    switch (m->aspen_svr_event_type)
    {
        case REQ:
            handle_req_rev_event(ns, b, m, lp);
            break;
        case ACK:
            handle_ack_rev_event(ns, b, m, lp);
            break;
        case KICKOFF:
            handle_kickoff_rev_event(ns, b, m, lp);
            break;
        case LOCAL:
            handle_local_rev_event(ns, b, m, lp);
            break;
        case DATA:
            handle_data_rev_event(ns, b, m, lp);
            break;
        case ASPENCOMP:
            handle_computation_rev_event(ns, b, m, lp);
            break;
        case RESTART:
            handle_restart_rev_event(ns, b, m, lp);
            break;
        default:
            fprintf(stderr, "\n Invalid reverse message type %d ", m->aspen_svr_event_type);
            assert(0);
            break;
    }
    
    return;
}

/* once the simulation is over, do some output */
static void aspen_svr_finalize(
    aspen_svr_state * ns,
    tw_lp * lp)
{
    if (debug_output >= 2)
    {
        printf("server %llu recvd %d bytes in %lf seconds, "\
               "%lf MiB/s sent_count %d recvd_count %d local_count %d \n", \
               (unsigned long long)(lp->gid/2),\
               payload_sz*ns->msg_recvd_count,\
               ns_to_s(ns->end_ts-ns->start_ts),\
               ((double)(payload_sz*num_reqs)/(double)(1024*1024)/ns_to_s(ns->end_ts-ns->start_ts)),\
               ns->msg_sent_count,\
               ns->msg_recvd_count,\
               ns->local_recvd_count);
    }
    return;
}

/* convert ns to seconds */
static tw_stime ns_to_s(tw_stime ns)
{
    return(ns / (1000.0 * 1000.0 * 1000.0));
}

/* convert seconds to ns */
static tw_stime s_to_ns(tw_stime ns)
{
    return(ns * (1000.0 * 1000.0 * 1000.0));
}

/* see declaration for more general info */
tw_lpid get_next_server(tw_lp *lp)
{
    tw_lpid rtn_id;
    /* first, get callers LP and group info from codes-mapping. Caching this
     * info in the LP struct isn't a bad idea for preventing a huge number of
     * lookups */
    char grp_name[MAX_NAME_LENGTH], lp_type_name[MAX_NAME_LENGTH],
    annotation[MAX_NAME_LENGTH];
    int  lp_type_id, grp_id, grp_rep_id, offset_num, num_reps;
    int  dest_rel_id;
    /* Grab necessary LP information about the sender: */
    codes_mapping_get_lp_info(lp->gid, grp_name, &grp_id, lp_type_name,
                              &lp_type_id, annotation, &grp_rep_id, &offset_num);
    /* Obtain the original server's relative id, increment by one, and convert back to lpid.
     * Also use modulo to avoid going out-of-bounds on the last relative id */
    switch (traffic_pattern_number){
        case NEXTNEIGHBOR:
            dest_rel_id = (codes_mapping_get_lp_relative_id(lp->gid, 0, 0) + 1) % num_servers;
            break;
        case RANDOM:
            dest_rel_id = tw_rand_integer(lp->rng, 0, num_servers - 2);
            if (dest_rel_id >= (codes_mapping_get_lp_relative_id(lp->gid, 0, 0)))
            {
                dest_rel_id ++;
            }
            assert(dest_rel_id < num_servers);
            break;
        default:
            fprintf(stderr, "ERROR: Invalid traffic pattern slipped through the shields, captain!\n");
            assert(0);
    }
    rtn_id = codes_mapping_get_lpid_from_relative(dest_rel_id, grp_name, lp_type_name, NULL, 0);
    /* Return the nearest neighbor lpid */
    return rtn_id;
}

/* handle initial event */
static void handle_kickoff_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    int dest_id;
    /* normally, when using ROSS, events are allocated as a result of the event
     * creation process. However, since we are now asking model-net to
     * communicate with an entity on our behalf, we need to generate both the
     * message to the recipient and an optional callback message
     * - thankfully, memory need not persist past the model_net_event call - it
     *   copies the messages */
    aspen_svr_msg m_local;
    aspen_svr_msg m_remote;
    
    m_local.aspen_svr_event_type = LOCAL;
    m_local.src = lp->gid;
    m_remote.aspen_svr_event_type = REQ;
    m_remote.src = lp->gid;
    
    /* record when transfers started on this server */
    ns->start_ts = tw_now(lp);
    
    dest_id = get_next_server(lp);
    
    /* model-net needs to know about
     * (1) higher-level destination LP which is usually a neighboring server in this case
     * (2) struct and size of remote message and
     * (3) struct and size of local message (a local message can be null) */
    model_net_event(net_id, "test", dest_id, payload_sz, 0, sizeof(aspen_svr_msg),\
                    (const void*)&m_remote, sizeof(aspen_svr_msg), (const void*)&m_local, lp);
    ns->msg_sent_count++;
}

static void handle_restart_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    int dest_id;
    /* normally, when using ROSS, events are allocated as a result of the event
     * creation process. However, since we are now asking model-net to
     * communicate with an entity on our behalf, we need to generate both the
     * message to the recipient and an optional callback message
     * - thankfully, memory need not persist past the model_net_event call - it
     *   copies the messages */
    aspen_svr_msg m_local;
    aspen_svr_msg m_remote;
    
    m_local.aspen_svr_event_type = LOCAL;
    m_local.src = lp->gid;
    m_remote.aspen_svr_event_type = REQ;
    m_remote.src = lp->gid;
    
    /* record when transfers restarted on this server */
    m->start_ts = ns->start_ts;
    ns->start_ts = tw_now(lp);
    /* adding this reset in for multi-round simulations */
    ns->msg_sent_count = 1;
    ns->msg_recvd_count = 0;
    
    dest_id = get_next_server(lp);
    
    model_net_event(net_id, "test", dest_id, payload_sz, 0, sizeof(aspen_svr_msg),\
                    (const void*)&m_remote, sizeof(aspen_svr_msg), (const void*)&m_local, lp);
}



/* at the moment, no need for local callbacks from model-net, so we maintain a
 * count for debugging purposes */
static void handle_local_event(
		aspen_svr_state * ns,
		tw_bf * b,
		aspen_svr_msg * m,
		tw_lp * lp)
{
    ns->local_recvd_count++;
}

/* handle recving ack
 * for this simulation, we repeatedly ping the destination server until num_reqs
 * of size payload_sz have been satisfied - we begin the next req when we
 * receive an ACK from the destination server */
static void handle_ack_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    /* the ACK actually doesn't come from the NIC on the other server -
     * model-net "hides" the NIC LP from us so we only see the original
     * destination server */
    
    if(ns->msg_sent_count < num_reqs)
    {
        m->incremented_flag = 1;
        /* again, allocate our own msgs so model-net can transmit on our behalf */
        aspen_svr_msg m_local;
        aspen_svr_msg m_remote;
        
        m_local.aspen_svr_event_type = LOCAL;
        m_local.src = lp->gid;
        m_remote.aspen_svr_event_type = REQ;
        m_remote.src = lp->gid;
        
        /* send another request */
        model_net_event(net_id, "test", m->src, payload_sz, 0, sizeof(aspen_svr_msg),\
                        (const void*)&m_remote, sizeof(aspen_svr_msg), (const void*)&m_local, lp);
        ns->msg_sent_count++;
    }
    else
    {
        /* threshold count reached, stop sending messages */
        m->incremented_flag = 0;
        m->end_ts = ns->end_ts;
        ns->end_ts = tw_now(lp);
        /* Send a message to LP 0 conatining your start and end times: */
        tw_event *e;
        aspen_svr_msg *msg;
        tw_stime data_time;
        
        // skew each data event slightly to help avoid event ties later on
        data_time = g_tw_lookahead + tw_rand_unif(lp->rng);
        
        // first create the event (time arg is an offset, not absolute time)
        e = codes_event_new(0, data_time, lp);
        // after event is created, grab the allocated message and set msg-specific\
        * data
        msg = tw_event_data(e);
        msg->src = lp->gid;
        msg->aspen_svr_event_type = DATA;
        msg->start_ts = ns->start_ts;
        msg->end_ts = ns->end_ts;
        // event is ready to be processed, send it off
        tw_event_send(e);
    }
    return;
}

/* handle receiving request */
static void handle_req_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    aspen_svr_msg m_local;
    aspen_svr_msg m_remote;
    
    m_local.aspen_svr_event_type = LOCAL;
    m_local.src = lp->gid;
    m_remote.aspen_svr_event_type = ACK;
    m_remote.src = lp->gid;
    
    ns->msg_recvd_count++;
    
    /* send ack back */
    /* simulated payload of 4 MiB */
    /* also trigger a local event for completion of payload msg */
    /* remote host will get an ack event */
    // TODO: see if the offset (0.000001) can be 0 or not.
    model_net_event(net_id, "test", m->src, payload_sz, 0, sizeof(aspen_svr_msg),\
                    (const void*)&m_remote, sizeof(aspen_svr_msg), (const void*)&m_local, lp);
    return;
}

static void handle_data_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    /* Make sure that the lp receiving this event is the 0 LP */
    assert(!lp->gid && !g_tw_mynode);
    if (debug_output >= 2)
    {
        printf("INFO: LP %lu received data event. (%u)\n", lp->gid, ns->data_recvd + 1);
    }
    if (m->start_ts < ns->start_global)
    {
        b->c0 = 1;
        swap_start(ns, m);
    }
    else
    {
        b->c0 = 0;
    }
    if (m->end_ts > ns->end_global)
    {
        b->c1 = 1;
        swap_end(ns, m);
    }
    else
    {
        b->c1 = 0;
    }
    ns->data_recvd ++;
    // When the last one has been received, send a self message for aspen computation
    if (ns->data_recvd == num_servers)
    {
        if (debug_output >= 2)
        {
            printf("\tLP %lu has received the last timestamp pair. Preparing for Aspen Comp.\n",\
                   lp->gid);
        }
        tw_event *e;
        aspen_svr_msg *msg;
        tw_stime compute_time;
        compute_time = g_tw_lookahead;
        
        // first create the event (time arg is an offset, not absolute time)
        e = codes_event_new(0, compute_time, lp);
        // after event is created, grab the allocated message and set msg-specific\
        * data
        msg = tw_event_data(e);
        msg->aspen_svr_event_type = ASPENCOMP;
        msg->src = lp->gid;
        // event is ready to be processed, send it off
        tw_event_send(e);
    }
    return;
}


static void handle_computation_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    tw_stime delta = 0;
    // Non-master LPs should never receive or send this event, so exit if they do.
    assert(!g_tw_mynode && !lp->gid && m->src == lp->gid);
    printf("INFO: Master LP %lu is now performing Aspen Computation %d\n",\
           lp->gid, roundsExecuted-computationRollbacks);
    
    /* Proceed with the computation: */
    m->incremented_flag = 0;
    delta += ns_to_s(ns->end_global - ns->start_global);
    
    printf("\tINFO: The network time elapsed is: %f s\n"\
           "\tThe start and end values are: %f ns and %f ns\n",\
           delta, ns->start_global, ns->end_global);
    
    delta += runtimeCalc(Aspen_App_Path[roundsExecuted - computationRollbacks],\
                         Aspen_Mach_Path, Aspen_Socket[roundsExecuted - computationRollbacks]);
    
    // Backup totalRuntime and then update it with the change from the past round
    m->end_ts = totalRuntime;
    totalRuntime += delta;
    
    printf("\tINFO: The final calculated runtime (up to this round) is %f seconds.\n",\
           totalRuntime);
    
    roundsExecuted ++;
    ns->data_recvd = 0;
    if (roundsExecuted < num_rounds + computationRollbacks)
    {
        if (debug_output >= 2){
            printf("INFO: sending restart messages to LPs now!\n");
        }
        m->incremented_flag = 1;
        /* There are more rounds to simulate, so send kickoffs to all LPs
         * Note that current_lpid is a relative lp id */
        tw_lpid global_id;
        /* Grab necessary LP information about the sender: */
        char grp_name[MAX_NAME_LENGTH], lp_type_name[MAX_NAME_LENGTH],
        annotation[MAX_NAME_LENGTH];
        int  lp_type_id, grp_id, grp_rep_id, offset_num;
        codes_mapping_get_lp_info(lp->gid, grp_name, &grp_id, lp_type_name,
                                  &lp_type_id, annotation, &grp_rep_id, &offset_num);
        /* Now send all the restarts: */
        tw_lpid relative_lpid = 0;
        for (; relative_lpid < num_servers; relative_lpid++)
        {
            tw_event *e;
            aspen_svr_msg *msg;
            tw_stime kickoff_time;
            /* skew each kickoff event slightly to help avoid event ties later on */
            kickoff_time = g_tw_lookahead + tw_rand_unif(lp->rng);
            /* Convert the relative id to global: */
            global_id = codes_mapping_get_lpid_from_relative(relative_lpid, grp_name, lp_type_name, NULL, 0);
            /* create the event (time arg is an offset, not absolute time) */
            e = codes_event_new(global_id, kickoff_time, lp);
            /* after event is created, grab the allocated message and set msg-specific
             * data */
            msg = tw_event_data(e);
            msg->aspen_svr_event_type = RESTART;
            msg->src = lp->gid;
            /* event is ready to be processed, send it off */
            tw_event_send(e);
            if (debug_output >= 2)
            {
                printf("INFO: Sent restart to lp %lu.\n", global_id);
            }
        }
    }
    return;
}


/* for us, reverse events are very easy, the only LP state that needs to be
 * rolled back are the counts.
 * for more complex simulations, this will not be the case (e.g., state
 * containing queues) */

/* reverse handler for kickoff */
static void handle_kickoff_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    ns->msg_sent_count--;
    model_net_event_rc(net_id, lp, payload_sz);
    switch (traffic_pattern_number){
        case NEXTNEIGHBOR:
            break;
        case RANDOM:
            tw_rand_reverse_unif(lp->rng);
            break;
        default:
            break;
    }
    return;
}

/* reverse handler for restart */
static void handle_restart_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    model_net_event_rc(net_id, lp, payload_sz);
    switch (traffic_pattern_number){
        case NEXTNEIGHBOR:
            break;
        case RANDOM:
            tw_rand_reverse_unif(lp->rng);
            break;
        default:
            break;
    }
    ns->start_ts = m->start_ts;
    ns->msg_sent_count = num_reqs;
    ns->msg_recvd_count = num_reqs;
}

static void handle_local_rev_event(
	       aspen_svr_state * ns,
	       tw_bf * b,
	       aspen_svr_msg * m,
	       tw_lp * lp)
{
    ns->local_recvd_count--;
}

/* reverse handler for ack*/
static void handle_ack_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    if(m->incremented_flag)
    {
        ns->msg_sent_count--;
        model_net_event_rc(net_id, lp, payload_sz);
    }
    else
    {
        ns->end_ts = m->end_ts;
        tw_rand_reverse_unif(lp->rng);
    }
    return;
}

/* reverse handler for req event */
static void handle_req_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    /* model-net has its own reverse computation support */ 
    model_net_event_rc(net_id, lp, payload_sz);
    ns->msg_recvd_count--;
    return;
}


/* reverse handler for data passing: */
static void handle_data_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    /* There's really not much to do here...just roll back to the \
     * previous start and end times and decrement the data_recvd counter. */
    assert(!lp->gid && !g_tw_mynode);
    if (debug_output >=2) printf("ROLLBACK: reversing data event. (%u)\n", ns->data_recvd - 1);
    if (b->c0)
    {
        swap_start(ns, m);
    }
    if (b->c1)
    {
        swap_end(ns, m);
    }
    ns->data_recvd --;
    if (ns->data_recvd == 0)
    {
        if (debug_output >= 2) printf("\tAll data events have been reversed.\n");
    }
    return;
}

/* reverse handler for aspen computation */
static void handle_computation_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp)
{
    assert(!lp->gid && !g_tw_mynode);
    
    printf("ROLLBACK: Performing reverse aspen computation.\n"\
           "\tCurrent value is: %f\n", totalRuntime);
    
    computationRollbacks ++;
    totalRuntime = m->end_ts;
    ns->data_recvd = num_servers; 
    if (totalRuntime < 0)
    {
        fprintf(stderr, "\tWARNING: after rollback totalRuntime was less than zero."\
                "Setting to zero.\n");
        totalRuntime = 0;
    }
    printf("\tAfter rollback, runtime value is: %f\n", totalRuntime);
    if (m->incremented_flag)
    {
        int i = 0;
        for ( ; i < num_servers; i++)
        {
            tw_rand_reverse_unif(lp->rng);
        }
    }
    return;
}


/*
 * Local variables:
 *  c-indent-level: 4
 *  c-basic-offset: 4
 * End:
 *
 * vim: ft=c ts=8 sts=4 sw=4 expandtab
 */
