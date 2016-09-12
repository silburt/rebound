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
double c_angle(struct reb_simulation* r, double r_ps, int i, int incoming);
double unaltered_energy(const struct reb_simulation* const r);

double E0, t_output, t_log_output, xyz_t = 0;
int xyz_counter = 0, numdt = 20;
char* mercury_dir; char* swifter_dir;

//temp
int output_xyz = 0; //switch to 0 for no outputs
time_t t_ini;
int N_prev;
char* name;
char output_name[100] = {0};
int *in_mini;
int N_CE = 0;
int L_CE = 0;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    double tmax = atof(argv[1]);
    int N_planetesimals = atoi(argv[2]);
    int rad_fac = atof(argv[3]);
    int seed = atoi(argv[4]);
    strcat(output_name,argv[5]); strcat(output_name,".txt"); name=argv[5];
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->ri_hermes.solar_switch_factor = 20.;         //X*radius
    r->testparticle_type = 1;
    r->heartbeat	= heartbeat;
    r->ri_hermes.hill_switch_factor = 3;        //Hill radii
    r->ri_hermes.adaptive_hill_switch_factor = 0;
    r->dt = 0.01;
    
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;     //switch to track the energy from collisions/ejections
    r->collision_resolve_keep_sorted = 1;
    
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
    
    //planet
    {
        double a=1, m=5e-5, e=0.01, inc=reb_random_normal(0.00001);
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p2.r = 2e-4*rad_fac;
        p2.hash = r->N;
        reb_add(r, p2);
    }
    
    r->N_active = r->N;
    
    //planetesimals
    //double total_planetesimal_mass = 500e-8;
    //double planetesimal_mass = total_planetesimal_mass/N_planetesimals;
    double planetesimal_mass = 1e-8;
    //double amin = 0.4, amax = 0.6;        //for planetesimal disk
    double amin = 0.8, amax = 1.2;
    double powerlaw = 0.5;
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(amin,amax,powerlaw);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        double inc = reb_random_normal(0.00001);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, 0., inc, Omega, apsis,phi);
		pt.r 		= 4e-5;
        pt.hash = r->N;
		reb_add(r, pt);
    }
    
    //com
    reb_move_to_com(r);
    
    //energy
    E0 = reb_tools_energy(r);
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,argv[5]); strcat(syss,"*");
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
    strcat(timeout,argv[5]); strcat(timeout,"_elapsedtime"); strcat(timeout,".txt");
    FILE* outt = fopen(timeout,"w");
    fprintf(outt,"\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    fprintf(outt,"System Parameters: dt=%f,tmax=%f,HSR=%f,N_planetesimals=%d,N_active=%d. \n",r->dt,tmax,r->ri_hermes.hill_switch_factor,N_planetesimals,r->N_active);
    fclose(outt);
    printf("\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    
}

void heartbeat(struct reb_simulation* r){
    if(r->t > t_output){//log output
        t_output = r->t*t_log_output;
        
        //collision corrected energy
        double E = reb_tools_energy(r);
        double dE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
        
        //unaltered energy
        double Eraw = unaltered_energy(r);
        double dEraw = fabs((Eraw-E0)/E0);
        
        time_t t_curr = time(NULL);
        struct tm *tmp2 = gmtime(&t_curr);
        double time = t_curr - t_ini;
        
        /*
        //
        //counting close encounters + calculate angle
        int count_CE = 0;
        if(count_CE){
            struct reb_simulation* mini = r->ri_hermes.mini;
            double r_ps = 0, angle = 0;
            char ia[200] = {0}; strcat(ia,name); strcat(ia,"_inangle.txt");
            char oa[200] = {0}; strcat(oa,name); strcat(oa,"_outangle.txt");
            struct reb_particle* global = r->particles;
            double dx = global[0].x - global[1].x;
            double dy = global[0].y - global[1].y;
            r_ps = sqrt(dx*dx + dy*dy);  //planet-star distance
            
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
                    angle = c_angle(r, r_ps, mini_id,1);
                    //printf("\nParticle entered, angle=%f, n_CE=%d,mini_id=%d\n",angle*57.29578,N_CE,mini_id);
                    N_CE++;
                    FILE* aa = fopen(ia, "a");
                    fprintf(aa,"%f,%f,%d,%d\n",r->t,angle,mini_id,N_CE);
                    fclose(aa);
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
                    angle = c_angle(r, r_ps, id,0);
                    //printf("\nParticle left, angle=%f, n_CE=%d, id=%d\n",angle*57.29578,N_CE,id);
                    FILE* aa = fopen(oa, "a");
                    fprintf(aa,"%f,%f,%d,%d\n",r->t,angle,id,N_CE);
                    fclose(aa);
                }
            }
        }*/
        
        FILE *append;
        append = fopen(output_name, "a");
        fprintf(append, "%f,%e,%e,%e,%e,%d,%d,%.1f,%d\n",r->t,dE,dEraw,r->energy_offset/E0,r->energy_offset,r->N,r->ri_hermes.mini->N,time,N_CE);
        fclose(append);
    }
    
    //record collisions in mini
    if(r->N < N_prev){
        char removed[200] = {0}; strcat(removed,name); strcat(removed,"_removed.txt");
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
                
                char removed[200] = {0}; strcat(removed,name); strcat(removed,"_removed"); strcat(removed,".txt");
                FILE* append = fopen(removed,"a");
                fprintf(append,"Ejection,%.5f\n",r->t);
                fclose(append);
                
                N_prev = r->N;
            }
        }
    }
}

//works only for one planet!
double c_angle(struct reb_simulation* r, double r_ps, int id, int incoming){
    struct reb_particle* p = r->particles;
    for(int i=r->N_active;i<r->N;i++){//get correct index
        if(p[i].hash == id){
            //atan2, sensitive to 2pi
            double x1 = p[1].x - p[0].x;
            double y1 = p[1].y - p[0].y;
            double phi1 = atan2(y1,x1);
            double x2 = p[i].x - p[1].x;
            double y2 = p[i].y - p[1].y;
            double phi2 = atan2(y2,x2);
            double phi;
            if(phi1 < M_PI) phi = (M_PI - phi1) + phi2; else phi = phi2 - (phi1 - M_PI);
            while(phi > 2*M_PI) phi -= 2*M_PI;
            while(phi < 0) phi += 2*M_PI;
        
            return phi;
        }
    }
    return -1;
}

double unaltered_energy(const struct reb_simulation* const r){
    const int N = r->N;
    const int N_var = r->N_var;
    const int _N_active = ((r->N_active==-1)?N:r->N_active) - N_var;
    const struct reb_particle* restrict const particles = r->particles;
    double e_kin = 0.;
    double e_pot = 0.;
    int N_interact = (r->testparticle_type==0)?_N_active:(N-N_var);
    for (int i=0;i<N_interact;i++){
        struct reb_particle pi = particles[i];
        e_kin += 0.5 * pi.m * (pi.vx*pi.vx + pi.vy*pi.vy + pi.vz*pi.vz);
    }
    for (int i=0;i<_N_active;i++){
        struct reb_particle pi = particles[i];
        for (int j=i+1;j<N_interact;j++){
            struct reb_particle pj = particles[j];
            double dx = pi.x - pj.x;
            double dy = pi.y - pj.y;
            double dz = pi.z - pj.z;
            e_pot -= r->G*pj.m*pi.m/sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
    
    return e_kin + e_pot;
}