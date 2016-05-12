#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double calc_a(struct reb_simulation* r, int i);

double E0;
char output_name[100] = {0};
double log_constant, tlog_output, lin_constant, tlin_output;
time_t t_ini;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    int seed = 14;
    srand(seed);
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HYBARID;
    r->ri_hybarid.switch_ratio = 2;  //Hill radii
    r->ri_hybarid.CE_radius = 15.;          //X*radius
    r->testparticle_type = 1;
	r->heartbeat	= heartbeat;
    r->dt = 0.03 * 6.2832;
    double tmax = 1000 * 6.283;
    
    //for this example if the boundaries are enforced then the itegration fails
    //r->usleep = 1000;
    //r->boundary	= REB_BOUNDARY_OPEN;
    //const double boxsize = 100;
    //reb_configure_box(r,boxsize,1,1,1);
    
    //r->collision = REB_COLLISION_DIRECT;
    //r->collision_resolve = reb_collision_resolve_merge;
    //r->collisions_track_dE = 1;
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
    
    double mm = 50;
    {
        double a=5.2, m=0.0009543*mm, e=0, inc = reb_random_normal(0.00001);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0.85);
        p1.r = 0.00046732617;
        p1.id = r->N;
        reb_add(r, p1);
    }
    {
        double a=9.5, m=0.0002857*mm, e=0.0, inc=reb_random_normal(0.00001);
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p2.r = 0.000389256877;
        p2.id = r->N;
        reb_add(r, p2);
    }
    {
        double a=19.2, m=0.00004365*mm, e=0.0, inc=reb_random_normal(0.00001);
        struct reb_particle p3 = {0};
        p3 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p3.r = 0.000169534499;
        p3.id = r->N;
        reb_add(r, p3);
    }
    {
        double a=30.1, m=5e-5*mm, e=0.0, inc=reb_random_normal(0.00001);
        struct reb_particle p4 = {0};
        p4 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p4.r = 1.6e-4;
        p4.id = r->N;
        reb_add(r, p4);
    }
    r->N_active = r->N;
    reb_move_to_com(r);
    
    int n_output = 25000;
    log_constant = pow(tmax + 1, 1./(n_output - 1));
    tlog_output = r->dt;
    lin_constant = tmax/n_output;
    tlin_output = r->dt;
    
    E0 = reb_tools_energy(r);
    
    //naming stuff
    char seedstr[15];
    sprintf(seedstr, "%d", seed);
    strcat(output_name,"output/outerSS"); strcat(output_name,"_sd"); strcat(output_name,seedstr);
    char timeout[200] = {0};
    strcat(timeout,output_name);
    strcat(output_name,".txt");
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,output_name);
    system(syss);
    time_t t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    
    //Integrate!
    reb_integrate(r, tmax);
    
    time_t t_fini = time(NULL);
    struct tm *tmp2 = gmtime(&t_fini);
    double time = t_fini - t_ini;
    strcat(timeout,"_elapsedtime.txt");
    FILE* outt = fopen(timeout,"w");
    fprintf(outt,"\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    fclose(outt);
    printf("\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    
}

void heartbeat(struct reb_simulation* r){
    if(r->t > tlog_output || r->t > tlin_output){//log output or linear output!!
        if(r->t > tlog_output)tlog_output = r->t*log_constant; else tlin_output = r->t+lin_constant;
        
        double E = reb_tools_energy(r);
        double dE = fabs((E-E0)/E0);
        
        time_t t_curr = time(NULL);
        struct tm *tmp2 = gmtime(&t_curr);
        double time = t_curr - t_ini;
        FILE *append;
        append = fopen(output_name, "a");
        fprintf(append, "%.16f,%.16f,%d,%.2f",r->t,dE,r->N,time);
        struct reb_particle* particles = r->particles;
        struct reb_particle p0 = particles[0];
        for(int i=0;i<r->N;i++){
            //double a = calc_a(r,i);
            struct reb_particle p = particles[i];
            double dx = p.x - p0.x;
            double dy = p.y - p0.y;
            double dz = p.z - p0.z;
            double a = sqrt(dx*dx + dy*dy + dz*dz);
            fprintf(append,",%.8f",a);
        }
        fprintf(append,"\n");
        fclose(append);
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
    }
    
}

double calc_a(struct reb_simulation* r, int i){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = reb_get_com(r);
    struct reb_particle p = particles[i]; //output planet only.
    const double mu = r->G*(com.m + p.m);
    const double dvx = p.vx-com.vx;
    const double dvy = p.vy-com.vy;
    const double dvz = p.vz-com.vz;
    const double dx = p.x-com.x;
    const double dy = p.y-com.y;
    const double dz = p.z-com.z;
    
    const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
    const double d = sqrt( dx*dx + dy*dy + dz*dz );    //distance
    const double dinv = 1./d;
    //const double muinv = 1./mu;
    //const double vr = (dx*dvx + dy*dvy + dz*dvz)*dinv;
    //const double term1 = v2-mu*dinv;
    //const double term2 = d*vr;
    //const double ex = muinv*( term1*dx - term2*dvx );
    //const double ey = muinv*( term1*dy - term2*dvy );
    //const double ez = muinv*( term1*dz - term2*dvz );
    //const double e = sqrt(ex*ex + ey*ey + ez*ez);   // eccentricity
    const double a = -mu/( v2 - 2.*mu*dinv );
    
    return a;
}
