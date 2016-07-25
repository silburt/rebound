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

void heartbeat(struct reb_simulation* r);

double E0, t_output, t_log_output, xyz_t = 0;
int xyz_counter = 0, numdt = 20;

//temp
int output_xyz = 0; //switch to 0 for no outputs
time_t t_ini;
int N_prev;
char* argv4;
char output_name[100] = {0};
int *in_mini;
int N_CE = 0;
int L_CE = 0;
int comr0 = 0;

//swifter/mercury compare
void output_to_mercury_swifter(struct reb_simulation* r, double HSR, double tmax, int n_output);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    int seed = atoi(argv[2]);
    strcat(output_name,argv[3]); strcat(output_name,".txt"); argv4=argv[3];
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->ri_hermes.radius_switch_factor = 20.;         //X*radius
    r->testparticle_type = 1;
    r->heartbeat	= heartbeat;
    r->ri_hermes.hill_switch_factor = 6;        //Hill radii
    r->dt = 0.015;
    double tmax = atof(argv[1]);
    
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->track_energy_offset = 1;     //switch to track the energy from collisions/ejections
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
   
    srand(seed);
    int n_output = 50000;
    t_log_output = pow(tmax + 1, 1./(n_output - 1));
    t_output = r->dt;
    printf("tlogoutput=%f\n",t_log_output);
    
    //planet 1
    double m_planet = 5e-5;
    {
        double a=1, m=m_planet, e=0.01, inc=reb_random_normal(0.00001);
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p2.r = 2e-4;
        p2.hash = r->N;
        reb_add(r, p2);
    }
    double rh = pow(m_planet/3.,1./3.);
    
    r->N_active = r->N;
    
    //planetesimals
    //double total_planetesimal_mass = 500e-8;
    //double planetesimal_mass = total_planetesimal_mass/N_planetesimals;
    int N_planetesimals = 200;
    double planetesimal_mass = 2.5e-7;
    double amin = 0.98, amax = 1.02;
    double powerlaw = 0.5;
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(amin,amax,powerlaw);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        double inc = reb_random_normal(0.00001);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, 0., inc, Omega, apsis,phi);
        /*
        double dx = pt.x - r->particles[1].x;
        double dy = pt.y - r->particles[1].y;
        double dz = pt.z - r->particles[1].z;
        while (pow(dx*dx+dy*dy+dz*dz,0.5) < rh*r->ri_hermes.hill_switch_factor){
            a	= reb_random_powerlaw(amin,amax,powerlaw);
            phi 	= reb_random_uniform(0,2.*M_PI);
            inc = reb_random_normal(0.00001);
            Omega = reb_random_uniform(0,2.*M_PI);
            apsis = reb_random_uniform(0,2.*M_PI);
            pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, 0., inc, Omega, apsis,phi);
            dx = pt.x - r->particles[1].x;
            dy = pt.y - r->particles[1].y;
            dz = pt.z - r->particles[1].z;
        }*/
		pt.r 		= 4e-5;
        pt.hash = r->N;
		reb_add(r, pt);
    }
    
    //com
    reb_move_to_com(r);
    
    //energy
    E0 = reb_tools_energy(r);
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,argv4); strcat(syss,"*");
    system(syss);
    
    t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    N_prev = r->N;
    
    //in_mini
    in_mini = calloc(sizeof(int),r->N);
    
    //Integrate!
    reb_integrate(r, tmax);
    
    //elapsed time stuff
    time_t t_fini = time(NULL);
    struct tm *tmp2 = gmtime(&t_fini);
    double time = t_fini - t_ini;
    char timeout[200] = {0};
    strcat(timeout,argv4); strcat(timeout,"_elapsedtime"); strcat(timeout,".txt");
    FILE* outt = fopen(timeout,"w");
    fprintf(outt,"\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    fprintf(outt,"System Parameters: dt=%f,tmax=%f,HSR=%f,N_planetesimals=%d,N_active=%d. \n",r->dt,tmax,r->ri_hermes.hill_switch_factor,N_planetesimals,r->N_active);
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
        
        //
        //counting close encounters + calculate angle
        int count_CE = 1;
        if(count_CE){
            struct reb_simulation* mini = r->ri_hermes.mini;
    
            int tempL_CE = L_CE;
            for(int i=mini->N_active;i<mini->N;i++){//check if entered
                int mini_id = mini->particles[i].hash;
                int found_in_mini = 0;
                for(int j=0;j<tempL_CE;j++)
                    if(in_mini[j] == mini_id){//already in in_mini?
                        found_in_mini = 1;
                    }
                if(found_in_mini == 0){//particle must have just entered CE region
                    in_mini[L_CE] = mini_id;
                    L_CE++;
                    N_CE++;
                }
            }
            
            int tempL_CE2 = L_CE;
            for(int i=0;i<tempL_CE2;i++){//check if left mini
                int id = in_mini[i];
                int found_in_mini = 0;
                for(int j=mini->N_active;j<mini->N;j++){
                    if(mini->particles[j].hash == id){
                        found_in_mini = 1;
                    }
                }
                if(found_in_mini == 0){//couldn't find in mini, must have left
                    L_CE--;
                    for(int k=i;k<tempL_CE2-1;k++) in_mini[k] = in_mini[k+1];
                    in_mini[tempL_CE2] = 0;
                }
            }
        }
        
        FILE *append;
        append = fopen(output_name, "a");
        fprintf(append, "%.16f,%.16f,%d,%d,%.1f,%d,%e\n",r->t,dE,r->N,r->ri_hermes.mini->N,time,N_CE,fabs(E-E0));
        fclose(append);
    }
    
    //record collisions in mini
    if(r->N < N_prev){
        char removed[200] = {0}; strcat(removed,argv4); strcat(removed,"_removed.txt");
        FILE* append = fopen(removed,"a");
        fprintf(append,"Collision,%.5f\n",r->t);
        fclose(append);
        
        N_prev = r->N;
    }
    
    //ejections
    {
        struct reb_particle* global = r->particles;
        const double ED2 = 25;
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
        }
    }
}