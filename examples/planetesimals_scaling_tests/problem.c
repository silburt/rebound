/**
 * A.S. This is my planetesimal disk with special integrator for close encounters.
 *      Particle id's: 0 = star, 1 = massive body, 2 = planetesimal, 3 = CLOSE ENCOUNTER
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define PI 3.14159265358979323846

void heartbeat(struct reb_simulation* r);
double E0, t_output, t_log_output;
char* mercury_dir; char* swifter_dir;

//temp
int output_xyz = 0, numdt=1, xyz_counter=0; double xyz_t = 36.76;//switch to 0 for no outputs
time_t t_ini;
int N_prev;
char* argv4;
char output_name[100] = {0};
int warning_message = 0;
double *in_mini;
int N_CE = 0;

//swifter/mercury compare
void output_to_mercury_swifter(struct reb_simulation* r, double HSR, double tmax, int n_output);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    //args
    int test_dt = 0;
    double HSF_or_dt = atof(argv[1]);  //HSF or dt depending on test_dt value
    strcat(output_name,argv[2]); strcat(output_name,".txt"); argv4=argv[2];
    double theta = atof(argv[3])*PI;    //planetesimal angle around the planet (0-pi)
    double f = atof(argv[4])*2*PI;      //intiial position of planet in orbit.
    
    int seed = 10;
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->ri_hermes.solar_switch_factor = 20.;         //X*radius
    r->ri_hermes.adaptive_hill_switch_factor = 0; 
    r->testparticle_type = 1;
    r->heartbeat	= heartbeat;
    //r->usleep = 300;
    
    double afac = 1;
    double mfac = 1;
    double mpfac = 1;
    
    //which test?
    double tmax;
    if(test_dt){//testing dt
        tmax = 1e5;
        //tmax = 10 * pow(afac,1.5);
        r->ri_hermes.hill_switch_factor = 6;        //units of Hill radii
        r->dt = HSF_or_dt * 6.28319;
    } else {//testing HSR
        //tmax = 10 * pow(afac,1.5);
        tmax = 7;
        r->ri_hermes.hill_switch_factor = HSF_or_dt;
        r->dt = 0.001;
    }
    
    //collision
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;     //switch to track the energy from collisions/ejections
    r->collision_resolve_keep_sorted = 1;
    
    //boundary
    //r->boundary			= REB_BOUNDARY_OPEN;
    //double boxsize = 5;
    //reb_configure_box(r, boxsize, 1, 1, 1);
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
   
    srand(seed);
    int n_output = 50000;
    t_log_output = pow(tmax + 1, 1./(n_output - 1));
    t_output = r->dt;
    printf("tlogoutput=%f",t_log_output);
    
    //planet 1
    double a=1*afac;
    {
        double m=5e-5*mfac, e=0, inc = reb_random_normal(0.00001);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, f);
        p1.r = 1.6e-4;              //radius of particle is in AU!
        reb_add(r, p1);
    }
    
    r->N_active = r->N;
    double planetesimal_mass = 1e-8*mpfac;
    {//planetesimal
        double rr = 0.005;
        //double vy = 1.1*pow(2*r->G*r->particles[1].m/rr,0.5)*sin(theta);
        //double vx = pow(2*r->G*r->particles[1].m/rr,0.5)*cos(theta);
        struct reb_particle pt = {0};
        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a + rr*cos(theta), 0, 0, 0, 0, f+rr*sin(theta));
        pt.vy += r->particles[1].vy/6.5;
        pt.vx += r->particles[1].vx/6.5;
        pt.r = 4e-5;
        reb_add(r, pt);
    }
    //printf("\np=%f, pl=%f\n\n", sqrt(pow(r->particles[1].x,2)+pow(r->particles[1].y,2)), sqrt(pow(r->particles[2].x,2)+pow(r->particles[2].y,2)));

    reb_move_to_com(r);
    /*
    {//orbiting around body 1 satellite
        double x=0.01;
        struct reb_particle pt = {0};
        pt = reb_tools_orbit_to_particle(r->G, r->particles[1], planetesimal_mass, x, 0, 0, 0, 0, 0);    //works well with m2=5e-4
        //pt.vy += 0.1*r->particles[1].vy;
        //pt.vx -= 0.1*r->particles[1].vx;
        pt.r = 4e-5;
        reb_add(r, pt);
    }*/
    
    /*
    {//orbiting around body 1 satellite
        double x=-0.05;
        struct reb_particle pt = {0};
        pt = reb_tools_orbit_to_particle(r->G, star, r->particles[1].m, r->particles[1].x, 0, 0, 0, 0, x);
        pt.m = planetesimal_mass;
        pt.r = 4e-5;
        reb_add(r, pt);
    }*/
    //
    
    //energy
    E0 = reb_tools_energy(r);
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,argv[2]); strcat(syss,"*");
    system(syss);
    
    t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    N_prev = r->N;
    argv4= argv[2];
    
    //malloc stuff for keeping track of particles entering/leaving mini
    in_mini = calloc(sizeof(double),r->N);
    
    //Integrate!
    reb_integrate(r, tmax);
    
    //elapsed time stuff
    time_t t_fini = time(NULL);
    struct tm *tmp2 = gmtime(&t_fini);
    double time = t_fini - t_ini;
    char timeout[200] = {0};
    strcat(timeout,argv[2]); strcat(timeout,"_elapsedtime"); strcat(timeout,".txt");
    FILE* outt = fopen(timeout,"w");
    fprintf(outt,"\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    fprintf(outt,"System Parameters: dt=%f,tmax=%f,HSR=%f,N_active=%d. \n",r->dt,tmax,r->ri_hermes.hill_switch_factor,r->N_active);
    fclose(outt);
    printf("\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    
}

void heartbeat(struct reb_simulation* r){
    if(r->t > t_output){//log output
        t_output = r->t*t_log_output;
        
        double E = reb_tools_energy(r);
        double dE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
        
        time_t t_curr = time(NULL);
        struct tm *tmp2 = gmtime(&t_curr);
        double time = t_curr - t_ini;
        
        for(int i=0;i<r->N;i++){
            if(in_mini[i] != r->ri_hermes.is_in_mini[i]){
                in_mini[i] = r->ri_hermes.is_in_mini[i];
                N_CE++;
            }
        }
        
        FILE *append;
        append = fopen(output_name, "a");
        fprintf(append, "%e,%e,%d,%d,%.1f,%d\n",r->t,dE,r->N,r->ri_hermes.mini->N);
        fclose(append);
    }
    
    //record collisions in mini
    if(r->N < N_prev){
        char removed[200] = {0}; strcat(removed,argv4); strcat(removed,"_removed"); strcat(removed,".txt");
        FILE* append = fopen(removed,"a");
        fprintf(append,"Collision,%.5f\n",r->t);
        fclose(append);
        
        N_prev = r->N;
    }
    
}

/*
 //planetesimals
 double amin = 0.95, amax = 1.05;        //for planetesimal disk
 double powerlaw = 1;
 while(r->N<N_planetesimals + r->N_active){
 struct reb_particle pt = {0};
 double a	= reb_random_powerlaw(amin,amax,powerlaw);
 double phi 	= reb_random_uniform(0,2.*M_PI);
 double inc = reb_random_normal(0.0001);
 double Omega = reb_random_uniform(0,2.*M_PI);
 double apsis = reb_random_uniform(0,2.*M_PI);
 pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, 0., inc, Omega, apsis,phi);
 pt.r 		= 4e-5;
 reb_add(r, pt);
 }*/

    /*
    if(warning_message == 0 && r->t > r->dt){
        //record max velocity
        struct reb_particle* particles = r->particles;
        double min_dt_enc2 = INFINITY;
        double hill_switch_factor2 = r->ri_hermes.hill_switch_factor*r->ri_hermes.hill_switch_factor;
        struct reb_particle p0 = particles[0];
        for(int i=1;i<r->N_active;i++){
            struct reb_particle pi = particles[i];
            const double dxi = p0.x - pi.x;
            const double dyi = p0.y - pi.y;
            const double dzi = p0.z - pi.z;
            const double r0i2 = dxi*dxi + dyi*dyi + dzi*dzi;
            const double mi = pi.m/(p0.m*3.);
            double rhi = pow(mi*mi*r0i2*r0i2*r0i2,1./6.);
            for(int j=i+1;j<r->N;j++){
                struct reb_particle pj = particles[j];
                
                const double dxj = p0.x - pj.x;
                const double dyj = p0.y - pj.y;
                const double dzj = p0.z - pj.z;
                const double r0j2 = dxj*dxj + dyj*dyj + dzj*dzj;
                const double mj = pj.m/(p0.m*3.);
                double rhj = pow(mj*mj*r0j2*r0j2*r0j2,1./6.);
                const double rh_sum = rhi+rhj;
                const double rh_sum2 = rh_sum*rh_sum;
                
                const double dx = pi.x - pj.x;
                const double dy = pi.y - pj.y;
                const double dz = pi.z - pj.z;
                const double rij2 = dx*dx + dy*dy + dz*dz;
                if(rij2 < hill_switch_factor2*rh_sum2){
                    const double dvx = pi.vx - pj.vx;
                    const double dvy = pi.vy - pj.vy;
                    const double dvz = pi.vz - pj.vz;
                    const double vij2 = dvx*dvx + dvy*dvy + dvz*dvz;
                    const double dt_enc2 = hill_switch_factor2*rh_sum2/vij2;
                    min_dt_enc2 = MIN(min_dt_enc2,dt_enc2);
                    if (warning_message==0 && min_dt_enc2 < 16.*r->dt*r->dt){
                        warning_message = 1;
                        char warning[200] = {0}; strcat(warning,argv4); strcat(warning,"_warning"); strcat(warning,".txt");
                        FILE* append = fopen(warning,"a");
                        fprintf(append,"The timestep is likely too large. Close encounters might be missed. Decrease the timestep or increase the switching radius. This warning will appear only once.\n");
                        fprintf(append,"t=%f: min_dt_enc=%f < 4*dt=%f\n",r->t,sqrt(min_dt_enc2),4*r->dt);
                        fclose(append);
                    }
                }
            }
        }
    }
    */
    /*
    //ejections
    const double ED2 = 100;
    struct reb_particle* global = r->particles;
    struct reb_particle p0 = global[0];
    for(int i=1;i<r->N;i++){
        const double dx = global[i].x - p0.x;
        const double dy = global[i].y - p0.y;
        const double dz = global[i].z - p0.z;
        if(dx*dx+dy*dy+dz*dz > ED2){
            const double Ei = reb_tools_energy(r);
            reb_remove(r,i,1);
            const double Ef = reb_tools_energy(r);
            r->energy_offset += Ei - Ef;
            
            char removed[200] = {0}; strcat(removed,argv4); strcat(removed,"_removed"); strcat(removed,".txt");
            FILE* append = fopen(removed,"a");
            fprintf(append,"Ejection,%.5f\n",r->t);
            fclose(append);
            
            N_prev = r->N;
        }
    }*/

/*
 {//planetesimal
 int highHSR = 1;
 int same_orbit = 1; //orbiting on same path as planet but slightly ahead
 double rr, term1, term2, f;
 if(highHSR){
 if(same_orbit){
 rr = 0;
 term1 = 0;
 term2 = r->particles[1].vy/25.;
 f = 0.1;
 }else{
 rr = 0.07;
 term1 = 0;
 term2 = 0;
 f=0;
 }
 } else {
 if(same_orbit){
 rr = 0;
 term1 = 0;
 term2 = r->particles[1].vy/3.1;
 f = 0.001;
 } else {
 rr = 0.001;
 term1 = pow(2*r->G*r->particles[1].m/rr,0.5);
 term2 = r->particles[1].vy/60;
 f=0;
 }
 }
 //double rr = .1;
 struct reb_particle pt = {0};
 pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, r->particles[1].x - rr*r->particles[1].x, 0, 0, 0, 0, f);
 pt.vy += term1 + term2;
 pt.r = 4e-5;
 pt.hash = r->N;
 reb_add(r, pt);
 }*/
