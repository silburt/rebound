#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double calc_a(struct reb_simulation* r);
double draw_ainv_powerlaw(double min, double max);

double E0;
char filename[100] = {0}, simname[300] = {0};
double log_constant, tlog_output, lin_constant, tlin_output, tbin_constant, tbin;
time_t t_ini;

int main(int argc, char* argv[]){
    strcat(filename,argv[1]);
    strcat(simname,filename); strcat(simname,".bin");
    struct reb_simulation* r;
    
    double tmax = 7e4 * 6.283;
    
    if(access(simname,F_OK) != -1){//if file exists
        r = reb_create_simulation_from_binary(simname);
        r->integrator	= REB_INTEGRATOR_HERMES;
        r->heartbeat	= heartbeat;
        r->testparticle_type = 1;
        printf("Found binary. Loaded snapshot at t=%.16f.\n",r->t);
    } else {
        printf("No binary file found. Creating new simulation.\n");
        r = reb_create_simulation();
        double m_earth = 0.000003003;
        int seed = atoi(argv[2]);
        srand(seed);
        
        //Simulation Setup
        r->integrator	= REB_INTEGRATOR_HERMES;
        r->heartbeat	= heartbeat;
        r->testparticle_type = 1;
        //r->dt = 2 * 6.283;  //planet's period = 125 years
        r->dt = 8;
        
        // Initial conditions
        struct reb_particle star = {0};
        star.m 		= 1;
        star.r		= 0.005;        // Radius of particle is in AU!
        reb_add(r, star);
        
        //planet 1
        double a1=25, m1=2.3*m_earth, e1=0, inc1=reb_random_normal(0.00001);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m1, a1, e1, inc1, 0, 0, 0);
        p1.r = 0.0000788215;       //radius of particle using 2g/cm^3 (AU)
        //p1.r = 5e-4;
        reb_add(r, p1);
        
        r->N_active = r->N;
        reb_move_to_com(r);
        
        //planetesimals
        double planetesimal_mass = m1/600.;     //each planetesimal = 1/600th of planet mass
        int N_planetesimals = 230.*m_earth/planetesimal_mass;
        //int N_planetesimals = 5000;
        double amin = a1 - 10.5, amax = a1 + 10.5;   //10.5AU on each side of the planet
        while(r->N<N_planetesimals + r->N_active){
            struct reb_particle pt = {0};
            double a	= reb_random_powerlaw(amin,amax,1);
            //double a = draw_ainv_powerlaw(amin,amax);
            double e = reb_random_rayleigh(0.01);   //rayleigh dist
            double inc = reb_random_rayleigh(0.005);
            double Omega = reb_random_uniform(0,2.*M_PI);
            double apsis = reb_random_uniform(0,2.*M_PI);
            double phi 	= reb_random_uniform(0,2.*M_PI);
            pt = reb_tools_orbit_to_particle(r->G, star, r->testparticle_type?planetesimal_mass:0., a, e, inc, Omega, apsis, phi);
            pt.r 		= 0.00000934532;
            reb_add(r, pt);
        }
    }
    
    //output rates
    int n_output = 10000;
    log_constant = pow(tmax + 1, 1./(n_output - 1));
    tlog_output = r->dt;
    lin_constant = tmax/n_output;
    tlin_output = r->dt;
    tbin_constant = tmax/1000;
    tbin = r->dt;
    
    //collisions
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
    // Boundaries
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 50;
    reb_configure_box(r,boxsize,3,3,1);
    
    //naming
    char timeout[300] = {0};
    strcat(timeout, filename); strcat(timeout,"_elapsedtime.txt");
    strcat(filename,".csv");
    
    E0 = reb_tools_energy(r);
    
    //Integrate!
    t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    heartbeat(r);
    reb_integrate(r, tmax);
    heartbeat(r);
    
    time_t t_fini = time(NULL);
    struct tm *tmp2 = gmtime(&t_fini);
    double time = t_fini - t_ini;

    FILE* outt = fopen(timeout,"a");
    fprintf(outt,"\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    fclose(outt);
    printf("\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    
}

void heartbeat(struct reb_simulation* r){
    if(r->t > tlog_output || r->t > tlin_output){//log output or linear output!!
        if(r->t > tlog_output)tlog_output = r->t*log_constant; else tlin_output = r->t+lin_constant;
        
        double E = reb_tools_energy(r);
        double dE = fabs((E-E0)/E0);
        double a1 = calc_a(r);
        
        time_t t_curr = time(NULL);
        struct tm *tmp2 = gmtime(&t_curr);
        double time = t_curr - t_ini;
        
        FILE *append;
        append = fopen(filename, "a");
        fprintf(append, "%f,%e,%e,%d,%.1f,%f\n",r->t,dE,a1,r->N,time,r->ri_hermes.hill_switch_factor);
        fclose(append);
        
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
    }
    
    if(r->t > tbin){
        tbin += tbin_constant;
        reb_output_binary(r, simname);
    }

}

double calc_a(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = reb_get_com(r);
    struct reb_particle p = particles[1]; //output planet only.
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

//returns value randomly drawn from P(x) = 1/x distribution
double draw_ainv_powerlaw(double min, double max){
    double y = reb_random_uniform(0., 1.);
    return exp(y*log(max/min) + log(min));
}
