/**
 * Migration and other orbit modifications
 *
 * This example shows how to add migration, eccentricity damping
 * and pericenter precession to a REBOUND simulation.  If you have
 * GLUT installed for visualization, press 'w' to see the orbits
 * as wires.  You can zoom out by holding shift, holding down the mouse
 * and dragging.  Press 'c' to better see migration/e-damping.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* sim);

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Setup constants
    sim->dt             = 0.012;        // initial timestep.
    sim->heartbeat = heartbeat;

    struct reb_particle p = {0}; 
    p.m     = 1.;   
    reb_add(sim, p); 

    double m = 0.;
    double a1 = 1.;
    double a2 = 2.;
    double e = 0.4;
    double omega = 0.;
    double f = 0.;

    struct reb_particle p1 = reb_tools_orbit2d_to_particle(sim->G, p, m, a1, e, omega, f);
    struct reb_particle p2 = reb_tools_orbit2d_to_particle(sim->G, p, m, a2, e, omega, f);
    reb_add(sim,p1);
    reb_add(sim,p2);
    reb_move_to_com(sim);

    struct rebx_extras* rebx = rebx_init(sim);

    // There are two options for how to modify orbits.  You would only choose one (comment the other out).  
    // You can't set precession separately with modify_orbits_forces (eccentricity and inclination damping induce pericenter and nodal precession).

    struct rebx_effect* params = rebx_add(rebx, "modify_orbits_direct");    					// directly update particles' orbital elements each timestep
    //struct rebx_effect* params = rebx_add(rebx, "modify_orbits_forces");    					// add forces that orbit-average to give exponential a and e damping

    // Set the timescales for each particle.  
    double tmax = 5.e4;
    
	double* tau_a = rebx_add_param(&sim->particles[1], "tau_a", REBX_TYPE_DOUBLE); 			// add semimajor axis damping on inner planet (e-folding timescale) 
	double* tau_omega = rebx_add_param(&sim->particles[1], "tau_omega", REBX_TYPE_DOUBLE);	// add linear precession (set precession period). Won't do anything for modify_orbits_forces
	double* tau_e = rebx_add_param(&sim->particles[2], "tau_e", REBX_TYPE_DOUBLE);			// add eccentricity damping on particles[2] (e-folding timescale)

    *tau_a = -tmax;
    *tau_omega = -tmax/10.;
    *tau_e = -tmax/10.;

    printf("Semimajor axis damping timescale for inner planet is %f.\n", fabs(*tau_a));
    printf("Precession timescale for inner planet is %f.\n", fabs(*tau_omega));
    printf("Eccentricity damping timescale for outer planet is %f.\n", fabs(*tau_e));

    /* One can also adjust a coupling parameter between eccentricity and semimajor axis damping.  We use the parameter p
     * as defined by Deck & Batygin (2015).  The default p=0 corresponds to no coupling, while p=1 corresponds to e-damping
     * at constant angular momentum.  This is only implemented for modify_orbits_direct.
	 * modify_orbits_forces damps eccentricity at constant angular momentum.
     *
     * Additionally, the damping by default is done in Jacobi coordinates.  If you'd prefer to use barycentric 
     * coordinates, or coordinates referenced to a particular particle, set a coordinates parameter in the effect
     * parameters returned by rebx_add to REBX_COORDINATES_BARYCENTRIC or REBX_COORDINATES_PARTICLE.  If the latter, add a 'primary' flag
     * to the reference particle (not neccesary for barycentric):
     */

	int* coordinates = rebx_add_param(params, "coordinates", REBX_TYPE_INT);
	*coordinates = REBX_COORDINATES_PARTICLE;

	int* primary = rebx_add_param(&sim->particles[0], "primary", REBX_TYPE_INT);
	*primary = 1;

    reb_integrate(sim, tmax);
    rebx_free(rebx);    // Free all the memory allocated by rebx
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    // output a e and pomega (Omega + omega) of inner body
    if(reb_output_check(sim, 5.e2)){
        const struct reb_particle sun = sim->particles[0];
        const struct reb_orbit orbit = reb_tools_particle_to_orbit(sim->G, sim->particles[1], sun); // calculate orbit of particles[1]
        printf("%f\t%f\t%f\t%f\n",sim->t,orbit.a, orbit.e, orbit.pomega);
    }
}
