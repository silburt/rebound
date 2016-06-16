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

void heartbeat(struct reb_simulation* r);
void migration_forces(struct reb_simulation* r);
void calc_resonant_angles(struct reb_simulation* r, FILE* f);

double E0;
char output_name[100] = {0};
time_t t_ini;
double* tau_a; 	/**< Migration timescale in years for all particles */
double* tau_e; 	/**< Eccentricity damping timescale in years for all particles */
double* omega;
double* lambda;
double* a;
double* e;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    double tmax = 1e5;
    srand(atoi(argv[1]));
    strcat(output_name,argv[2]);
    
	// Simulation Setup
	r->integrator	= REB_INTEGRATOR_WHFAST;
    r->heartbeat	= heartbeat;
    r->additional_forces = migration_forces;
    r->force_is_velocity_dependent = 1;
    r->dt = 0.01;
    
    //Migration constants
    double t_mig = 2e4;
    double K = 100.0;       //Lee & Peale (2002) K.
    double e_ini = 0.01;
    
	// Star
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
    
    // Planet 1
    {
        double a=1, m=5e-4, inc=reb_random_normal(0.00001);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e_ini, inc, 0, 0, 0);
        p.r = 0.000467;
        reb_add(r, p);
    }
    
    //Planet 2
    {
        double a=1.75, m=5e-4, inc=reb_random_normal(0.00001);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e_ini, inc, 0, 0, 0);
        p.r = 0.000467;
        reb_add(r, p);
    }
    
    r->N_active = r->N;
    
    //migration stuff
    tau_a = calloc(sizeof(double),r->N);
    tau_e = calloc(sizeof(double),r->N);
    tau_a[2] = 2.*M_PI*t_mig;
    tau_e[2] = 2.*M_PI*t_mig/K;
    a = calloc(sizeof(double),r->N);
    e = calloc(sizeof(double),r->N);
    omega = calloc(sizeof(double),r->N);
    lambda = calloc(sizeof(double),r->N);
    
    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    
    //naming
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,output_name); strcat(syss,"*");
    system(syss);
    strcat(output_name,".txt");
    
    // Integrate!
    reb_integrate(r, tmax);
    
    printf("\nSimulation complete. \n\n");
}

double tout = 0;
void heartbeat(struct reb_simulation* r){
    if (tout <r->t){
        tout += 100;
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        int N_mini = 0;
        if (r->ri_hermes.mini_active){
            N_mini = r->ri_hermes.mini->N;
        }
        
        FILE* f = fopen(output_name, "a");
        fprintf(f,"%e,%e,%d,%d,",r->t,relE,r->N,N_mini);
        calc_resonant_angles(r,f);
        fclose(f);
    }
    
    if (reb_output_check(r, 100.*r->dt)){
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("%e",relE);
    }
}

void calc_resonant_angles(struct reb_simulation* r, FILE* f){
    struct reb_particle com = reb_get_com(r);
    for(int i=1;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        const double mu = r->G*(r->particles[0].m + p.m);
        const double dvx = p.vx-com.vx;
        const double dvy = p.vy-com.vy;
        const double dvz = p.vz-com.vz;
        const double dx = p.x-com.x;
        const double dy = p.y-com.y;
        const double dz = p.z-com.z;
        
        const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
        const double d = sqrt ( dx*dx + dy*dy + dz*dz );
        const double vr = (dx*dvx + dy*dvy + dz*dvz)/d;
        const double ex = 1./mu*( (v*v-mu/d)*dx - d*vr*dvx );
        const double ey = 1./mu*( (v*v-mu/d)*dy - d*vr*dvy );
        const double ez = 1./mu*( (v*v-mu/d)*dz - d*vr*dvz );
        e[i] = sqrt( ex*ex + ey*ey + ez*ez );   // eccentricity
        a[i] = -mu/(v*v - 2.*mu/d);
        const double rdote = dx*ex + dy*ey + dz*ez;
        const double cosf = rdote/(e[i]*d);
        
        omega[i] = atan2(ey,ex);
        if(ey < 0.) omega[i] += 2*M_PI;
        double cosE = (a[i] - d)/(a[i]*e[i]);
        double E;
        if(cosf > 1. || cosf < -1.){
            E = M_PI - M_PI*cosE;
        } else {
            E = acos(cosE);
        }
        if(vr < 0.) E = 2.*M_PI - E;
        double MA = E - e[i]*sin(E);
        lambda[i] = MA + omega[i];
    }
    double phi = 0, phi2 = 0, phi3 = 0;
    phi = 2.*lambda[2] - lambda[1] - omega[1];
    phi2 = 2.*lambda[2] - lambda[1] - omega[2];
    phi3 = omega[1] - omega[2];
    while(phi >= 2*M_PI) phi -= 2*M_PI;
    while(phi < 0.) phi += 2*M_PI;
    while(phi2 >= 2*M_PI) phi2 -= 2*M_PI;
    while(phi2 < 0.) phi2 += 2*M_PI;
    while(phi3 >= 2*M_PI) phi3 -= 2*M_PI;
    while(phi3 < 0.) phi3 += 2*M_PI;
    
    fprintf(f,"%f,%f,%f,%f,%f,%f,%f\n",a[1],e[1],a[2],e[2],phi,phi2,phi3);
    
}

void migration_forces(struct reb_simulation* r){
    const double G = r->G;
    const int N = r->N;
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = reb_get_com(r);
    for(int i=1;i<N;i++){
        if (tau_e[i]!=0||tau_a[i]!=0){
            struct reb_particle* p = &(particles[i]);
            const double dvx = p->vx-com.vx;
            const double dvy = p->vy-com.vy;
            const double dvz = p->vz-com.vz;
            
            if (tau_a[i]!=0){ 	// Migration
                p->ax -=  dvx/(2.*tau_a[i]);
                p->ay -=  dvy/(2.*tau_a[i]);
                p->az -=  dvz/(2.*tau_a[i]);
            }
            
            if (tau_e[i]!=0){ 	// Eccentricity damping
                const double mu = G*(com.m + p->m);
                const double dx = p->x-com.x;
                const double dy = p->y-com.y;
                const double dz = p->z-com.z;
                
                const double hx = dy*dvz - dz*dvy;
                const double hy = dz*dvx - dx*dvz;
                const double hz = dx*dvy - dy*dvx;
                const double h = sqrt ( hx*hx + hy*hy + hz*hz );
                const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
                const double r = sqrt ( dx*dx + dy*dy + dz*dz );
                const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
                const double ex = 1./mu*( (v*v-mu/r)*dx - r*vr*dvx );
                const double ey = 1./mu*( (v*v-mu/r)*dy - r*vr*dvy );
                const double ez = 1./mu*( (v*v-mu/r)*dz - r*vr*dvz );
                const double e = sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity
                const double a = -mu/( v*v - 2.*mu/r );			// semi major axis
                const double prefac1 = 1./(1.-e*e) /tau_e[i]/1.5;
                const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))  /tau_e[i]/1.5;
                p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
                p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
                p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;
            }
        }
        com = reb_get_com_of_pair(com,particles[i]);
    }
}