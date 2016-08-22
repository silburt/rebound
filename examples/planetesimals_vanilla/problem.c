//This is essentially the Kirsch example, using it as a vanilla example, need to vary various parameters and make sure the scaling relations of Armitage work.  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double calc_a(struct reb_simulation* r, int index);
double draw_ainv_powerlaw(double min, double max);

double E0;
int N_prev;
char output_name[100] = {0};
char removed[200] = {0};
char* argv4;
double log_constant, tlog_output, lin_constant, tlin_output;
time_t t_ini;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    double m_earth = 0.000003003;
    //double m_neptune = 0.00005149;
    
    double powerlaw = 0;
    int seed = 10;
    srand(seed);
    strcat(output_name,"output/test");
    argv4 = "output/test";
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->ri_hermes.hill_switch_factor = 6;         //Hill radii
    r->ri_hermes.adaptive_hill_switch_factor = 0;
    r->ri_hermes.solar_switch_factor = 20.;     //X*radius
    r->testparticle_type = 1;
	r->heartbeat	= heartbeat;
    double tmax = 1e5 * 6.283;
    
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 25;
    reb_configure_box(r,boxsize,2,2,1);
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
    
    //does order of adding bodies matter? Yes, apparently it does says Hanno. Jacobi coords.
    
    //inner massive planet to scatter planetesimals out
    double a2=3, m2=5e-4, e2=0, inc2=reb_random_normal(0.00001);
    struct reb_particle p1 = {0};
    p1 = reb_tools_orbit_to_particle(r->G, star, m2, a2, e2, inc2, 0, 0, 0);
    p1.r = 0.000467;       //radius of Jupiter (AU)
    p1.hash = r->N;
    reb_add(r, p1);
    
    //planet 2
    double a1=5, m1=2.3*m_earth, e1=0, inc1=reb_random_normal(0.00001);
    struct reb_particle p2 = {0};
    p2 = reb_tools_orbit_to_particle(r->G, star, m1, a1, e1, inc1, 0, 0, 0);
    p2.r = 0.0000788215;       //radius of particle using 2g/cm^3 (AU)
    //p2.r = 5e-4;
    p2.hash = r->N;
    reb_add(r, p2);

    r->N_active = r->N;
    r->dt = pow(a2,1.5)/50;
    
    //planetesimals-what's a reasonable ini? Perfectly cold disk seems unlikely...
    //double planetesimal_mass = m1/600;     //each planetesimal = 1/600th of planet mass
    //int N_planetesimals = 230.*m_earth/planetesimal_mass;
    double total_disk_mass = m1;
    int N_planetesimals = 100;
    double planetesimal_mass = total_disk_mass / N_planetesimals;
    //double amin = a1 - 2, amax = a1 + 2;          //planet in center of disk
    //double amin = a1 - 0.5, amax = a1 + 2.5;      //planet in asymmetric disk
    double amin = a1-0.25, amax = a1 + 0.25;                //planet at edge of disk
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(amin,amax,powerlaw);
        //double a = draw_ainv_powerlaw(amin,amax);
        double e = reb_random_rayleigh(0.2);   //rayleigh dist
        double inc = reb_random_rayleigh(0.005);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, r->testparticle_type?planetesimal_mass:0., a, e, inc, Omega, apsis, phi);
		pt.r 		= 0.00000934532;
        pt.hash = r->N;
		reb_add(r, pt);
    }
    
    int n_output = 25000;
    log_constant = pow(tmax + 1, 1./(n_output - 1));
    tlog_output = r->dt;
    lin_constant = tmax/n_output;
    tlin_output = r->dt;

    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    N_prev = r->N;
    
    //naming stuff
    char timeout[200] = {0};
    strcat(timeout,output_name);
    strcat(removed,output_name); strcat(removed,"_removed.txt");
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,output_name); strcat(syss,"*");
    system(syss);
    strcat(output_name,".txt");
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
    //Center if ejected particle
    if(r->N < N_prev){
        N_prev = r->N;
        double E = reb_tools_energy(r);
        reb_move_to_com(r);
        r->energy_offset += E - reb_tools_energy(r);
    }
    
    if(r->t > tlog_output || r->t > tlin_output){//log output or linear output!!
        if(r->t > tlog_output)tlog_output = r->t*log_constant; else tlin_output = r->t+lin_constant;
        
        double E = reb_tools_energy(r);
        double dE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
        double a1 = calc_a(r,1);
        double a2 = 0;
        if(r->N_active > 2) a2 = calc_a(r,2);
        
        int N_mini = 0;
        if(r->integrator == REB_INTEGRATOR_HERMES) N_mini = r->ri_hermes.mini->N - r->ri_hermes.mini->N_active;
        
        time_t t_curr = time(NULL);
        struct tm *tmp2 = gmtime(&t_curr);
        double time = t_curr - t_ini;
        FILE *append;
        append = fopen(output_name, "a");
        fprintf(append, "%e,%e,%e,%f,%d\n",r->t,dE,time,r->ri_hermes.hill_switch_factor,r->N);
        fclose(append);
        
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
    }
}

double calc_a(struct reb_simulation* r, int index){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = reb_get_com(r);
    struct reb_particle p = particles[index]; //output planet only.
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
    const double a = -mu/( v2 - 2.*mu*dinv );
    
    return a;
}

//returns value randomly drawn from P(x) = 1/x distribution
double draw_ainv_powerlaw(double min, double max){
    double y = reb_random_uniform(0., 1.);
    return exp(y*log(max/min) + log(min));
}


