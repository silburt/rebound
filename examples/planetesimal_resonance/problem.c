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
double calc_a(struct reb_simulation* r, int index);
void calc_resonant_angles(struct reb_simulation* r, FILE* f);

double E0;
char output_name[100] = {0};
time_t t_ini;
double tout = 0;

//binary save
char binary_output_name[100] = {0};
double binary_output_time;

int main(int argc, char* argv[]){
    char binary[100] = {0}; strcat(binary, argv[1]); strcat(binary,".bin");
    struct reb_simulation* r = reb_create_simulation_from_binary(binary);
    int N_planetesimals = 5000;
    srand(atoi(argv[2]));
    strcat(output_name,argv[1]); strcat(output_name,"planetesimals_sd"); strcat(output_name,argv[2]);
    
	// Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->heartbeat	= heartbeat;
    r->ri_hermes.hill_switch_factor = 3;
    r->ri_hermes.radius_switch_factor = 20.;
    r->testparticle_type = 1;
    //r->gravity_ignore_10 = 0; //Use if created binary with WHFAST but using non-WHFAST now.
    r->dt = 0.01;
    double tmax = 1e6;
    tout = r->t;
    
    // Collisions
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
    // Boundaries
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 5;
    reb_configure_box(r,boxsize,2,2,1);
    
    // Planetesimal disk parameters (Planets already added)
    double total_disk_mass = r->particles[1].m;
    double planetesimal_mass = total_disk_mass/N_planetesimals;
    printf("%e,%e\n",total_disk_mass,planetesimal_mass);
    double amin = calc_a(r, 1) - 0.5, amax = calc_a(r, 2) + 0.5;
    double powerlaw = 0;
    
    // Generate Planetesimal Disk
    struct reb_particle star = r->particles[0];
    while(r->N<(N_planetesimals + r->N_active)){
		struct reb_particle pt = {0};
		double a    = reb_random_powerlaw(amin,amax,powerlaw);
        double e    = reb_random_rayleigh(0.005);
        double inc  = reb_random_rayleigh(0.005);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, r->testparticle_type?planetesimal_mass:0., a, e, inc, Omega, apsis, phi);
		pt.r 		= 0.00000934532;
        pt.hash = r->N;
		reb_add(r, pt);
    }

    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    
    //binary (temp)
    strcat(binary_output_name, output_name); strcat(binary_output_name, "_t=");
    binary_output_time = r->t;
    
    //naming
    char timeout[200] = {0};
    strcat(timeout,output_name);
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,output_name); strcat(syss,"*");
    system(syss);
    strcat(output_name,".txt");
    time_t t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    
    // Integrate!
    reb_integrate(r, tmax);
    
    //final outputs
    time_t t_fini = time(NULL);
    struct tm *tmp2 = gmtime(&t_fini);
    double time = t_fini - t_ini;
    strcat(timeout,"_elapsedtime.txt");
    FILE* outt = fopen(timeout,"w");
    fprintf(outt,"\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    fprintf(outt, "Simulation Details:\n");
    fprintf(outt, "Setup Parmaeters: HSF=%.2f,RSF=%.1f,dt=%e,tmax=%f\n",r->ri_hermes.hill_switch_factor,r->ri_hermes.radius_switch_factor,r->dt,tmax);
    fprintf(outt, "Planetesimals: N_pl=%d,M_pl_tot=%e,amin_pl=%f,amax_pl=%f,powerlaw=%.2f\n",N_planetesimals,total_disk_mass,amin,amax,powerlaw);
    fclose(outt);
    printf("\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
}

void heartbeat(struct reb_simulation* r){
    if (tout <r->t){
        //tout += 0.01;
        tout += 10;
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
    
    //output binary, temp
    if(binary_output_time < r->t){
        char out_time[10] = {0}; sprintf(out_time,"%.0f",r->t);
        char out[200] = {0}; strcat(out, binary_output_name); strcat(out, out_time); strcat(out, ".bin");
        reb_output_binary(r, out);
        binary_output_time += 1;
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

void calc_resonant_angles(struct reb_simulation* r, FILE* f){
    double e[3] = {0};
    double a[3] = {0};
    double omega[3] = {0};
    double lambda[3] = {0};
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