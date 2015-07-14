/*
 * Aspen-Codes integration model
 * Author: Mark Blanco
 * Date: 24 June 2015
 * Note: this model is directly based off of the codes-base example (example.c)
 * They may look very similar for a while
 * The CODES project can be found at http://www.mcs.anl.gov/research/projects/codes/
 * The ROSS project, which codes is based on, can be found at 
 * https://github.com/carothersc/ROSS/wiki
 * Finally, ASPEN can be found at https://ft.ornl.gov/node/1566...sort of.
 */

#ifndef ASPEN_HEADER
#define ASPEN_HEADER

#include <string.h>
#include <assert.h>
#include <ross.h>

// NOTE: the config file should contain the paths to the aspen kernel model(s)
// and hardware model, as well as the socket on which the kernel model should be
// evaluated.

static int socket;     // Global int for the socket to be used (perhaps this should be configure in the conf file?)
char Aspen_Mach_Path[100];      // Global for path and name of Aspen model
char **Aspen_App_Path = NULL;        // Global array for paths and names of Aspen application/kernels
char **Aspen_Socket;         // Global array for names of sockets to be used
// TODO: Do something better than having a hard-coded length...
static int num_reqs = 0;/* number of requests sent by each server (read from config) */
static int payload_sz = 0; /* size of simulated data payload, bytes (read from config) */
static int num_rounds = 0; /* number of computation-simulation rounds to perform (read from config) */

/* model-net ID, can be either simple-net, dragonfly or torus (more may be
 * added) */
static int net_id = 0;
static int num_servers = 0;
static int offset = 2;

static unsigned int ttl_lps = 0;

/* expected LP group name in configure files for this program */
static char *group_name = "ASPEN_SERVERS";
/* expected parameter group name for rounds of communication */
static char *param_group_nm = "server_pings";
static char *misc_param_gp_nm = "PARAMS";
static char *num_reqs_key = "num_reqs";
static char *payload_sz_key = "payload_sz";
static char *aspen_group_nm = "ASPEN_PARAMS";
static char aspen_app_key[] = "aspen_app_path000";
static char *aspen_mach_key = "aspen_mach_path";
static char aspen_socket_key[] = "socket_choice000";
static char *num_rounds_key = "num_rounds";

typedef struct svr_msg aspen_svr_msg;
typedef struct svr_state aspen_svr_state;

/* types of events that will constitute server activities */
enum svr_event
{
    KICKOFF,    /* initial event */
    RESTART,    /* multi-round kickoff */
    REQ,        /* request event */
    ACK,        /* ack event */
    LOCAL,      /* local event */
    ASPENCOMP,  /* event during which Aspen will be called */
    DATA        /* Send start and end timestamps to the 0 LP */
};

/* this struct serves as the ***persistent*** state of the LP representing the 
 * server in question. This struct is setup when the LP initialization function
 * ptr is called */
struct svr_state
{
    int msg_sent_count;   /* requests sent */
    int msg_recvd_count;  /* requests recvd */
    int local_recvd_count; /* number of local messages received */
    tw_stime start_ts;    /* time that we started sending requests */
    tw_stime end_ts;      /* time that last request finished */
    tw_stime start_global;  /* global earliest start time */
    tw_stime end_global;    /* global latest end time */
    /* Note that the globals above are obviously not global,
     * but need to be stored here for purposes of reverse event handling */ 
    unsigned int data_recvd; /* counter for data sent to LP 0 */
};

/* this struct serves as the ***temporary*** event data, which can be thought
 * of as a message between two LPs. */
struct svr_msg
{
    enum svr_event aspen_svr_event_type;
    tw_lpid src;          /* source of this request or ack */

    int incremented_flag; /* helper for reverse computation */
    tw_stime start_ts;    /* storage of start time for data */
    tw_stime end_ts;      /* storage of end time for data */
};

/* ROSS expects four functions per LP:
 * - an LP initialization function, called for each LP
 * - an event processing function
 * - a *reverse* event processing function (rollback), and
 * - a finalization/cleanup function when the simulation ends
 */
static void aspen_svr_init(
    aspen_svr_state * ns,
    tw_lp * lp);
static void aspen_svr_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void aspen_svr_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void aspen_svr_finalize(
    aspen_svr_state * ns,
    tw_lp * lp);

/* set up the function pointers for ROSS, as well as the size of the LP state
 * structure (NOTE: ROSS is in charge of event and state (de-)allocation) */
tw_lptype svr_lp = {
    (init_f) aspen_svr_init,
    (pre_run_f) NULL,
    (event_f) aspen_svr_event,
    (revent_f) aspen_svr_rev_event,
    (final_f)  aspen_svr_finalize,
    (map_f) codes_mapping,
    sizeof(aspen_svr_state),
};

extern const tw_lptype* svr_get_lp_type();
static void svr_add_lp_type();
static tw_stime ns_to_s(tw_stime ns);
static tw_stime s_to_ns(tw_stime ns);

/* as we only have a single event processing entry point and multiple event
 * types, for clarity we define "handlers" for each (reverse) event type */
static void handle_kickoff_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void handle_restart_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void handle_ack_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void handle_req_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void handle_local_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
   tw_lp * lp);
static void handle_data_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void handle_computation_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void handle_local_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
   tw_lp * lp);
static void handle_kickoff_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void handle_restart_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void handle_ack_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void handle_req_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void handle_data_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);
static void handle_computation_rev_event(
    aspen_svr_state * ns,
    tw_bf * b,
    aspen_svr_msg * m,
    tw_lp * lp);


/* for this simulation, each server contacts its neighboring server in an id.
 * this function shows how to use the codes_mapping API to calculate IDs when
 * having to contend with multiple LP types and counts. Note that in this simple
 * example codes_mapping is overkill. */
static tw_lpid get_next_server(tw_lpid sender_id);

/* arguments to be handled by ROSS - strings passed in are expected to be
 * pre-allocated */
static char conf_file_name[256] = {0};
/* this struct contains default parameters used by ROSS, as well as
 * user-specific arguments to be handled by the ROSS config sys. Pass it in
 * prior to calling tw_init */

/* two value-swapper functions for processing the start and end 
 * timestamps received in data events */
inline void swap_start(aspen_svr_state * ns, aspen_svr_msg * m)
{
    tw_stime temp = ns->start_global;
    ns->start_global = m->start_ts;
    m->start_ts = temp;
}

inline void swap_end(aspen_svr_state * ns, aspen_svr_msg * m)
{
    tw_stime temp = ns->end_global;
    ns->end_global = m->end_ts;
    m->end_ts = temp;
}

const tw_optdef app_opt [] =
{
	TWOPT_GROUP("Model net test case" ),
        TWOPT_CHAR("conf", conf_file_name, "name of codes configuration file"),
	TWOPT_END()
};

/* Helper function to return a stringified version of an int */
int int_to_array(int num, char** array){
    int temp = num;
    int count = 0;
    while (temp > 0){
        temp /= 10;
        count ++;
    }
    if (count == 0) count = 1;
    *array = calloc(count, sizeof(char));
    temp = count;
    while (count > 0){
        count --;
        (*array)[count] = (num % 10) + '0';
        num /= 10;
    }
    return temp;
}
#endif
