/**
 * Radiation forces on a debris disk
 * 
 * This example shows how to integrate a debris disk around a Sun-like star, with
 * dust particles under the action of radiation forces using WHFast. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* r);

double tmax = 3e12;                     // in sec, ~ 10^5 yrs

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Setup constants
    double AU = 1.5e11;                 // in meters
    sim->integrator     = REB_INTEGRATOR_WHFAST;
    sim->G              = 6.674e-11;    // Use SI units
    sim->dt             = 1e8;          // At ~100 AU, orbital periods are ~1000 yrs, so here we use ~1% of that, in sec
    sim->N_active       = 1;            // The dust particles do not interact with one another gravitationally
    sim->heartbeat      = heartbeat;
    
    // sun
    struct reb_particle sun = {0};
    sun.m  = 1.99e30;                   // mass of Sun in kg
    reb_add(sim, sun);
    
    struct rebx_extras* rebx = rebx_init(sim);
    struct rebx_effect* rad_params = rebx_add(rebx, "radiation_forces");
    double* c = rebx_add_param(rad_params, "c", REBX_TYPE_DOUBLE);
    *c = 3.e8;                      // speed of light in SI units

    /* The effect assumes the radiation source is particles[0].  You can set it to another one by adding a radiation_source flag to it:
     * int* source = rebx_add_param(&sim->particles[3], "radiation_source", REBX_TYPE_INT);
     * *source = 1;*/

    /* Initialize dust particles
     * We idealize a perfectly coplanar debris disk with particles that have semimajor axes between 100 and 120 AU.
     * We initialize particles with 0.01 eccentricity, random pericenters and azimuths, and uniformly distributed
     * semimajor axes in the above range.  We take all particles to have a beta parameter (ratio of radiation 
     * force to gravitational force from the star) of 0.1.
     *
     ***** Only particles that have their beta parameter set will feel radiation forces. *****/

    double amin = 100.*AU;
    double awidth = 20.*AU;
    double e = 0.1;
    double Ndust = 1000;                // Number of dust particles

    int seed = 3;                       // random number generator seed
    srand(seed);
    double a, pomega, f;
    struct reb_particle p;
    double* beta;
    for(int i=1; i<=Ndust; i++){
                                        // first we set up the orbit from sets of uniformly drawn orbital elements.  
        a = amin + awidth*(double)rand() / (double)RAND_MAX;
        pomega = 2*M_PI*(double)rand() / (double)RAND_MAX;
        f = 2*M_PI*(double)rand() / (double)RAND_MAX;
        
        double m=0;                     // We treat the dust grains as massless.
        p = reb_tools_orbit2d_to_particle(sim->G, sim->particles[0], m, a, e, pomega, f); 
        reb_add(sim, p); 
                                        // Only particles with beta set will feel radiation forces 
        beta = rebx_add_param(&sim->particles[i], "beta", REBX_TYPE_DOUBLE);    
        *beta = 0.1;
    }

    reb_move_to_com(sim);

    reb_integrate(sim, tmax);

    /* Note that the debris disk will seem to undergo oscillations at first.  
     * This is actually the correct behavior.  We set up the
     * particles with orbital elements (that assume only gravity is acting).  When we turn on radiation,
     * all of a sudden particles are moving too fast for the reduced gravity they feel (due to radiation
     * pressure), so they're all at pericenter.  All the particles therefore move outward in phase.  It
     * takes a few cycles before the differential motion randomizes the phases. */

    rebx_free(rebx);                                /* free memory allocated by REBOUNDx */
    reb_free_simulation(sim);
}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 1.e8)){
        reb_output_timing(sim, tmax);
    }
}
