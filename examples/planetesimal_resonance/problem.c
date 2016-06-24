/**
 * Planetesimal Disk Migration
 *
 * This example integrates a star, 2 planet, N planetesimal disk system, with the
 * outer planet at the inner edge of the planetesimal disk. If the system is
 * integrated for at least 10^5 years outward migration by the outer planet in
 * the planetesimal disk will be observed. By default, the semi-major axis of both
 * planets along with the fractional energy error are printed to energy.txt.
 *
 * The ideal integrator choice for this problem is HERMES due to the large number
 * of close encounters.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include <time.h>
double tout = 0;
void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    char binary[100] = {0}; strcat(binary, argv[1]); strcat(binary,".bin");
    struct reb_simulation* r = reb_create_simulation_from_binary(binary);
    struct reb_particle p = r->particles[0];
    reb_reset_temporary_pointers(r);
    printf("%f %f %f %f %f %f \n",p.m,p.y,p.z,p.vx,p.vy,p.vz);
    
	// Simulation Setup
    r->integrator	= REB_INTEGRATOR_LEAPFROG;
    r->additional_forces = NULL;
    r->N_active = 2;
    r->gravity_ignore_10 = 0;
    r->dt = 0.01;
    double tmax = 1e8;
    r->usleep = 1e5;
    r->heartbeat = heartbeat;
  
  
    reb_move_to_com(r);
    reb_remove(r,2,0);
    
    // Integrate!
    reb_integrate(r, tmax);
}

void heartbeat(struct reb_simulation* r){
struct reb_particle p = r->particles[1];
printf("%f %f %f %f %f %f \n",p.m,p.y,p.z,p.ax,p.ay,p.az);
}
