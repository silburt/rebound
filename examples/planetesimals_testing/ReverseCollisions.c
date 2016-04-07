/**
 * A.S. This is my planetesimal disk with special integrator for close encounters.
 *      Particle id's: 0 = star, 1 = massive body, 2 = planetesimal, 3 = CLOSE ENCOUNTER
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double c_angle(struct reb_simulation* r, double r_ps, int i, int incoming);
double c_dcom(struct reb_simulation* r);
void c_orb(struct reb_simulation* r, int planet_index, double* a, double* e, double* w, double* sinf, double* cosf, double* ey);
void apsidal_precession(struct reb_simulation* r, int index, double dwdt);

double E0, t_output, t_log_output, xyz_t = 0;
int xyz_counter = 0, numdt = 20;
char* mercury_dir; char* swifter_dir;

//temp
int output_xyz = 0; //switch to 0 for no outputs
time_t t_ini;
int N_prev;
char* argv4;
char output_name[100] = {0};
int *in_mini;
int N_CE = 0;
int L_CE = 0;
struct reb_particle comr0;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation_from_binary("ReverseCollisions.bin");
    
    strcat(output_name,"output/ReverseCollisions");
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HYBARID;
    r->ri_hybarid.CE_radius = 20.;         //X*radius
    r->testparticle_type = 1;
    r->heartbeat	= heartbeat;
    r->ri_hybarid.switch_ratio = 6;        //Hill radii
    r->dt = -0.001;
    double tmax = 0;

    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->collisions_track_dE = 1;     //switch to track the energy from collisions/ejections
    r->usleep = 10000;
    
    t_log_output = 1.00048;
    t_output = r->t;

    r->N_active = 2;
    reb_move_to_com(r);
    comr0 = reb_get_com(r);
    
    //energy
    E0 = reb_tools_energy(r);
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,argv[3]); strcat(syss,"*");
    system(syss);
    
    t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    N_prev = r->N;
    
    //Integrate!
    reb_integrate(r, tmax);
    
}

void heartbeat(struct reb_simulation* r){
    if(r->t < t_output){//log output
        t_output = r->t/t_log_output;
        
        double E = reb_tools_energy(r) + r->collisions_dE;
        //double dE = fabs((E-E0)/E0);
        double dE = (E-E0)/E0;
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
        
        time_t t_curr = time(NULL);
        struct tm *tmp2 = gmtime(&t_curr);
        double time = t_curr - t_ini;
        
        int N_mini = 0;
        //if(r->integrator == REB_INTEGRATOR_HYBARID) N_mini = r->ri_hybarid.mini->N - r->ri_hybarid.mini->N_active;
        
        FILE *append;
        
        char outputt[100] = {0}; strcat(outputt,output_name);
        append = fopen(strcat(outputt,".txt"), "a");
        fprintf(append, "%.16f,%.16f,%d,%d,%.1f,%d,%e,%e,%e,0,0,%e\n",r->t,dE,r->N,N_mini,time,N_CE,E-E0,E,E0,c_dcom(r));
        fclose(append);
    }
    
    //record collisions in mini
    if(r->N < N_prev){
        char removed[100] = {0}; strcat(removed,output_name);
        FILE* append = fopen(strcat(removed,"_removed.txt"),"a");
        fprintf(append,"Collision,%.5f\n",r->t);
        fclose(append);
        
        N_prev = r->N;
    }
}

double c_dcom(struct reb_simulation* r){
    struct reb_particle com = reb_get_com(r);
    double dx = comr0.x - com.x;
    double dy = comr0.y - com.y;
    double dz = comr0.z - com.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}
