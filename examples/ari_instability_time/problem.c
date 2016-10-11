/**
 * Instability Time
 *
 * I'm curious to find out whether the time of first Hill sphere crossing and the time of first 
 * collision is the same. People use Hill sphere since it's easy to simulate, just use WH, but 
 * is this a good assumption? Use Dan's setup as the default.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include <time.h>

void heartbeat(struct reb_simulation* r);
void calc_resonant_angles(struct reb_simulation* r, FILE* f);
double calc_a(struct reb_simulation* r, int index);

double E0, tout;
char output_name[200] = {0};
time_t t_ini;
int N_prev;     //keep track of ejected particles

//eia snapshot
char CE[250] = {0};

//binary save
char binary_out[100] = {0};

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    //args
    strcat(output_name,argv[1]); strcat(output_name,"_sd"); strcat(output_name,argv[2]);
    int seed = atoi(argv[2]);
    srand(seed);
    
	// Standard Simulation Setup
    r->integrator	= REB_INTEGRATOR_HERMES;
    r->heartbeat	= heartbeat;
    r->testparticle_type = 1;
    
    //Important parameters
    r->ri_hermes.hill_switch_factor = 3;
    r->ri_hermes.solar_switch_factor = 20.;
    r->dt = 0.001;
    double tmax = 5e7;
    
    // Collisions
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
    // Boundaries
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 1;
    reb_configure_box(r,boxsize,1,1,1);
    
	// Star
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
    star.hash = 0;
	reb_add(r, star);
    
    //params
    double m_planet = 5.0*3e-6;
    double r_planet = 2*0.00004258689;
    double imax = 1.*M_PI/180.;
    double emax = 2e-2;
    double gamma = pow((2./3.*m_planet),1./3.);
    double betamin = 5;
    double betamax = 10;
    double a_inner = 0.05;
    double a_middle = 0;
    
    // Planet 1
    {
        double a=a_inner, m=m_planet, inc=reb_random_uniform(0,imax), e=reb_random_uniform(0,emax);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, reb_random_uniform(0,2*M_PI), reb_random_uniform(0,2*M_PI), reb_random_uniform(0,2*M_PI));
        p.r = r_planet;
        p.hash = r->N;
        reb_add(r, p);
    }
    
    //Planet 2
    {
        double a=a_inner+gamma*reb_random_uniform(betamin,betamax), m=m_planet, inc=reb_random_uniform(0,imax), e=reb_random_uniform(0,emax);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, reb_random_uniform(0,2*M_PI), reb_random_uniform(0,2*M_PI), reb_random_uniform(0,2*M_PI));
        p.r = r_planet;
        p.hash = r->N;
        reb_add(r, p);
        a_middle = a;
    }
    
    //Planet 3
    {
        double a=a_middle+gamma*reb_random_uniform(betamin,betamax), m=m_planet, inc=reb_random_uniform(0,imax), e=reb_random_uniform(0,emax);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, reb_random_uniform(0,2*M_PI), reb_random_uniform(0,2*M_PI), reb_random_uniform(0,2*M_PI));
        p.r = r_planet;
        p.hash = r->N;
        reb_add(r, p);
    }
    
    r->N_active = r->N;
    
    N_prev = r->N;
    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    tout = r->dt;
    
    //naming
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,output_name); strcat(syss,"*");
    char binary[200] = {0}; strcat(binary,output_name); strcat(binary, "t=0.bin");
    strcat(CE,output_name); strcat(CE,"_CE.csv");
    system(syss);
    strcat(output_name,".csv");
    
    // Integrate!
    reb_output_binary(r, binary);
    reb_integrate(r, tmax);

}

void heartbeat(struct reb_simulation* r){
    //output values to file
    if (tout <r->t){
        tout *= 1.005;
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        int N_mini = 0;
        if (r->ri_hermes.mini_active){
            N_mini = r->ri_hermes.mini->N;
        }
        
        FILE* f = fopen(output_name, "a");
        fprintf(f,"%e,%e,%d,%d,%e,",r->t,relE,r->N,N_mini,r->ri_hermes.current_hill_switch_factor);
        calc_resonant_angles(r,f);
        fclose(f);
    }
    
    //Check if close encounters
    double gamma2 = pow(r->particles[1].m/3.,2./3.);
    struct reb_particle p0 = r->particles[0];
    for(int i=1;i<r->N;i++){
        struct reb_particle pi = r->particles[i];
        const double dxi = p0.x - pi.x;
        const double dyi = p0.y - pi.y;
        const double dzi = p0.z - pi.z;
        const double r0i2 = dxi*dxi + dyi*dyi + dzi*dzi;
        double rhi2 = r0i2*gamma2;
        for(int j=i+1;j<r->N;j++){
            struct reb_particle pj = r->particles[j];
            const double dxj = p0.x - pj.x;
            const double dyj = p0.y - pj.y;
            const double dzj = p0.z - pj.z;
            const double r0j2 = dxj*dxj + dyj*dyj + dzj*dzj;
            double rhj2 = r0j2*gamma2;
            
            const double dx = pi.x - pj.x;
            const double dy = pi.y - pj.y;
            const double dz = pi.z - pj.z;
            const double rij2 = dx*dx + dy*dy + dz*dz;
            if(rij2<(rhj2+rhi2)){
                FILE* f = fopen(CE, "a");
                fprintf(f,"%e,%d,%d,%e,%e\n",r->t,i,j,rhi2,rhj2);
                fclose(f);
            }
        }
    }
    
    //Exit if ejected or collided particle
    if(r->N < N_prev){
        r->status = REB_EXIT_ESCAPE;
    }
    
    //output on screen
    if (reb_output_check(r, 100.*r->dt)){
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("%e",relE);
    }
}

void calc_resonant_angles(struct reb_simulation* r, FILE* f){
    double m[4] = {0};
    double e[4] = {0}; //Hardcoding, probably should change in the future.
    double a[4] = {0};
    double omega[4] = {0};
    double lambda[4] = {0};
    struct reb_particle com = reb_get_com(r);
    for(int i=1;i<r->N_active;i++){
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
        m[i] = p.m;
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
    
    for(int i=1;i<r->N;i++){
        fprintf(f,"%e,%f,%f,",m[i],a[i],e[i]);
    }
    fprintf(f,"%f,%f,%f\n",phi,phi2,phi3);
    
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

