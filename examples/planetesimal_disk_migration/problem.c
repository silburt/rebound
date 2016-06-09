/**
 * Planetesimal Disk Migration
 *
 * This example integrates a star, 2 planet, N planetesimal disk system, with the
 * outer planet at the inner edge of the planetesimal disk. If the system is
 * integrated for at least 10^4 years, outward migration of the outer planet in
 * the planetesimal disk will be observed. By default, the semi-major axis of both
 * planets along with the fractional energy error are printed to energy.txt.
 *
 * Since many close encounters occur the ideal integrator choice for this problem 
 * is HERMES.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double calc_a(struct reb_simulation* r, int index);
double E0;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
	// Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->heartbeat	= heartbeat;
    r->ri_hermes.hill_switch_factor = 3;
    r->ri_hermes.radius_switch_factor = 20.;
    r->testparticle_type = 1;
    
    // Collisions
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    
    // Boundaries
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 12;
    reb_configure_box(r,boxsize,2,2,1);
    
    srand(12);
    
	// Star
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
    
    // Planet 1 - inner massive planet to scatter planetesimals out
    double a1=2, m1=5e-4, e1=0, inc1=reb_random_normal(0.00001);
    struct reb_particle p1 = {0};
    p1 = reb_tools_orbit_to_particle(r->G, star, m1, a1, e1, inc1, 0, 0, 0);
    p1.r = 0.000467;
    reb_add(r, p1);
    
    // Planet 2 - outer smaller planet to migrate in the disk
    double m_earth = 0.000003003;
    double a2=5, m2=2.3*m_earth, e2=0, inc2=reb_random_normal(0.00001);
    struct reb_particle p2 = {0};
    p2 = reb_tools_orbit_to_particle(r->G, star, m2, a2, e2, inc2, 0, 0, 0);
    p2.r = 0.0000788215;
    reb_add(r, p2);
    
    r->N_active = r->N;
    r->dt = pow(a1,1.5)/30;
    
    // Planetesimal disk parameters
    double total_disk_mass = m2*10.;
    int N_planetesimals = 1000;
    double planetesimal_mass = total_disk_mass/N_planetesimals;
    double amin = a2, amax = a2 + 2;                //planet at edge of disk
    double powerlaw = 0;
    
    // Generate Planetesimal Disk
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a    = reb_random_powerlaw(amin,amax,powerlaw);
        double e    = reb_random_rayleigh(0.005);
        double inc  = reb_random_rayleigh(0.005);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, r->testparticle_type?planetesimal_mass:0., a, e, inc, Omega, apsis, phi);
		pt.r 		= 0.00000934532;
		reb_add(r, pt);
    }

    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    system("rm -f energy.txt");
    
    // Integrate!
    reb_integrate(r, INFINITY);
}

double tout = 0.1;
void heartbeat(struct reb_simulation* r){
    if (tout <r->t){
        tout *=1.01;
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        FILE* f = fopen("energy.txt","a+");
        int N_mini = 0;
        if (r->ri_hermes.mini_active){
            N_mini = r->ri_hermes.mini->N;
        }
        fprintf(f,"%e,%e,%f,%f,%d,%d\n",r->t,relE,calc_a(r,1),calc_a(r,2),r->N,N_mini);
        fclose(f);
    }
    
    if (reb_output_check(r, 100.*r->dt)){
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("%e",relE);
    }
}

double calc_a(struct reb_simulation* r, int index){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = reb_get_com(r);
    struct reb_particle p = particles[index];
    const double mu = r->G*(com.m + p.m);
    const double dvx = p.vx-com.vx;
    const double dvy = p.vy-com.vy;
    const double dvz = p.vz-com.vz;
    const double dx = p.x-com.x;
    const double dy = p.y-com.y;
    const double dz = p.z-com.z;
    
    const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
    const double d = sqrt(dx*dx + dy*dy + dz*dz);    //distance
    const double dinv = 1./d;
    const double a = -mu/(v2 - 2.*mu*dinv);
    
    return a;
}