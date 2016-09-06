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

//temp
char* argv4;
char output_name[100] = {0};
time_t t_ini;
int N_prev;
int N_ini;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    strcat(output_name,"output/test.txt");
    int seed = 15;
    srand(seed);
    
    //Simulation Setup
    r->integrator	= REB_INTEGRATOR_HERMES;
    //r->ri_hermes.hill_switch_factor = 1;         //Hill radii
    //r->ri_hermes.solar_switch_factor = 20.;     //X*radius
    r->testparticle_type = 1;
    r->heartbeat	= heartbeat;
    double tmax = 1e6 * 6.283;
    
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 12;
    reb_configure_box(r,boxsize,2,2,1);
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
    
    //massive bodies
    double J2AU = 0.00046731951;    //conversion factor from jupiter radius to AU
    double d2r = 0.0174533;
    /*
     {
         double m=4.84568e-5, rad=1.*J2AU;
         double a=9.703237608614e-2, e=2.626647049513e-2, inc=4.767341227875e-1*d2r;
         double g=1.249385857281e2*d2r, n=2.397754395330e2*d2r, M=1.719987629785e2*d2r;
         struct reb_particle p1 = {0};
         p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, n, g, M);
         p1.r = rad;
         p1.hash = r->N;
         reb_add(r, p1);
     }
    {
        double m=4.84877e-5, rad=1.*J2AU;
        double a=1.986664416121e-1, e=1.284618351955e-2, inc=5.235482090217e-1*d2r;
        double g=9.207489467772e1*d2r, n=1.247973464163e2*d2r, M=2.632085892267e2*d2r;
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, n, g, M);
        p1.r = rad;
        p1.hash = r->N;
        reb_add(r, p1);
     
    }*/
    {
        double m=1.74984e-6, rad=1.*J2AU;
        double a=1.9, e=0.1, inc=5.026524644629e-1*d2r;
        double g=1.874461205587e2*d2r, n=2.946138921044e2*d2r, M=reb_random_uniform(0,2.*M_PI);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, n, g, M);
        p1.r = rad;
        p1.hash = r->N;
        reb_add(r, p1);
        
        //set dt
        r->dt = pow(a,1.5)/30;
    }
    {
        double m=3.04157e-6, rad=1.*J2AU;
        double a=2.515, e=0.2, inc=8.856892519209e-1*d2r;
        double g=3.023854601867e2*d2r, n=1.421158173436e2*d2r, M=reb_random_uniform(0,2.*M_PI);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, n, g, M);
        p1.r = rad;     //planet radius, AU
        p1.hash = r->N;
        reb_add(r, p1);
        
    }
    {
        double m=9.71267e-5, rad=1.*J2AU;
        double a=3.311503223441, e=0.1, inc=2.017989263302e-1*d2r;
        double g=1.740259990382e2*d2r, n=1.389578681864e2*d2r, M=reb_random_uniform(0,2.*M_PI);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, n, g, M);
        p1.r = rad;
        p1.hash = r->N;
        reb_add(r, p1);
    }
    {
        double m=5.08788e-6, rad=1.*J2AU;
        double a=4.383122343529, e=0.25, inc=1.310554328953e-1*d2r;
        double g=1.773477245097e2*d2r, n=1.405393882599e2*d2r, M=reb_random_uniform(0,2.*M_PI);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, n, g, M);
        p1.r = rad;
        p1.hash = r->N;
        reb_add(r, p1);
    }

    r->N_active = r->N;
    reb_move_to_com(r);
    N_prev = r->N;
    N_ini = r->N;
    
    char sys[200] = {0}; strcat(sys,"rm -f "); strcat(sys,"output/test*.txt");
    system(sys);
    
    //calculate initial energy
    E0 = reb_tools_energy(r);
    
    t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    
    //Integrate!
    reb_integrate(r, tmax);
    
}

double tout = 0.1;
void heartbeat(struct reb_simulation* r){
    //Center if ejected particle
    if(r->N < N_prev){
        N_prev = r->N;
        double E = reb_tools_energy(r);
        reb_move_to_com(r);
        r->energy_offset += E - reb_tools_energy(r);
    }
    
    if (tout <r->t){
        tout += 100;
        //tout += 1;
        double E = reb_tools_energy(r);
        double dE = fabs((E-E0)/E0);
        FILE* f = fopen(output_name,"a+");
        
        time_t t_curr = time(NULL);
        struct tm *tmp2 = gmtime(&t_curr);
        double time = t_curr - t_ini;
        
        int N_mini = 0;
        if (r->ri_hermes.mini_active){
            N_mini = r->ri_hermes.mini->N;
        }
        fprintf(f,"%e,%e,%d,%d,%f,%f",r->t,dE,N_mini,r->N,time,r->ri_hermes.hill_switch_factor);
        
        double a = 0;
        for(int i=1;i<N_ini;i++){
            if(i<r->N) a = calc_a(r,i);
            fprintf(f,",%e",a);
            a = 0;
        }
        fprintf(f,"\n");
        fclose(f);
    }
    
    if (reb_output_check(r, 100.*r->dt)){
        double E = reb_tools_energy(r);
        double dE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        int N_mini = 0;
        if (r->ri_hermes.mini_active){
            N_mini = r->ri_hermes.mini->N;
        }
        printf("dE=%e,N_mini=%d",dE,N_mini);
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