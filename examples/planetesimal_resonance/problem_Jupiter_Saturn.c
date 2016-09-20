/**
 * Planets in resonance embedded in a Planetesimal Disk
 *
 * In this example Jupiter and Saturn are migrated into resonance with planetesimals. The option for
 * a warm or cold start is possible. If warm start, there is one simulation with the planets 
 * migrating simultaneously with the planetesimal disk already added. If cold start, I've automated
 * the code such that first the planets migrate into resonance, the simulation is stopped, 
 * saved as a binary, then planetesimals are added, then the simulation starts again.
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
void migration_forces(struct reb_simulation* r);
void calc_resonant_angles(struct reb_simulation* r, FILE* f);
double calc_a(struct reb_simulation* r, int index);
void eia_snapshot(struct reb_simulation* r, char* time);

double E0;
char output_name[200] = {0};
time_t t_ini;
double* tau_a; 	/**< Migration timescale in years for all particles */
double* tau_e; 	/**< Eccentricity damping timescale in years for all particles */
double mig_time;
int rescale_energy;
int N_prev;     //keep track of ejected particles

//eia snapshot
char eia_out[100] = {0};
double output_time;
double output_inc;

//binary save
char binary_out[100] = {0};

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
//**********************Parameter List**********************
    int warm_start = 1;                 //if 0, cold start, if 1, warm start
    int N_planetesimals = atoi(argv[1]);
    int seed = atoi(argv[2]);
    strcat(output_name,argv[3]);        //name
    
    //Migration parameters
    mig_time = 10000;
    double mig_rate = 5e4;              //mini Jupiters = 2e4, Neptunes = 5e4
    double K = 10.0;                   //Lee & Peale (2002) K.
    rescale_energy = 0;                 //after migration is done, rescale energy
    
    //Planetesimal disk parameters
    double m_earth = 3e-6;
    double total_disk_mass = m_earth*50;
    double a_extend = 0.2;              //extension of planetesimal disk beyond planet orbits (AU)
    double planetesimal_mass = total_disk_mass/N_planetesimals;
    double powerlaw = 1;
    double N_eia_binary_outputs = 7;    //number of eia_snapshot and binary output intervals
//**********************************************************
    
	// Standard Simulation Setup
    r->integrator	= REB_INTEGRATOR_HERMES;
    r->heartbeat	= heartbeat;
    r->additional_forces = migration_forces;
    r->force_is_velocity_dependent = 1;
    r->testparticle_type = 1;
    
    //Important parameters
    r->ri_hermes.hill_switch_factor = 3;
    r->ri_hermes.solar_switch_factor = 20.;
    r->dt = 0.01;
    double tmax = 1e6;
    
    // Collisions
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
    // Boundaries
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 10;
    reb_configure_box(r,boxsize,1,1,1);
    
    srand(seed);
    output_inc = tmax/N_eia_binary_outputs;
    output_time = output_inc;
    
	// Star
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
    star.hash = 0;
	reb_add(r, star);
    
    //Jupiter
    {
        double a_real = 5.2;
        double a=1, m=0.0009543, e=0, inc=reb_random_normal(0.00001);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, reb_random_uniform(0,2*M_PI));
        p.r = 0.00046732617;
        p.hash = r->N;
        reb_add(r, p);
    }
    
    //Saturn
    {
        double a_real = 9.54;
        double a=1.64, m=0.0002857, e=0, inc=reb_random_normal(0.00001);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, reb_random_uniform(0,2*M_PI));
        p.r = 0.000389256877;
        p.hash = r->N;
        reb_add(r, p);
    }
    
    r->N_active = r->N;
    double amin=calc_a(r,1)-a_extend, amax=calc_a(r,2)+a_extend;    //planetesimal disk min/max
    
    //migration arrays
    tau_a = calloc(sizeof(double),r->N_active);
    tau_e = calloc(sizeof(double),r->N_active);
    tau_a[2] = 2.*M_PI*mig_rate;
    tau_e[2] = 2.*M_PI*mig_rate/K;
    
    //naming
    if(warm_start){strcat(output_name,"_warm");} else {strcat(output_name,"_cold");}
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,output_name); strcat(syss,"*");
    char info[200] = {0}; strcat(info,output_name);
    strcat(binary_out,output_name); strcat(binary_out,"_t=");
    strcat(eia_out, output_name); strcat(eia_out, "_eiasnapshot_t=");
    system(syss);
    strcat(output_name,".txt");
    
    //info output
    strcat(info,"_info.txt");
    FILE* out1 = fopen(info,"w");
    fprintf(out1, "Simulation Details:\n");
    int coll_on = 0; if(r->collision_resolve == reb_collision_resolve_merge) coll_on =1;
    fprintf(out1, "\nSetup Parmaeters:\nHSF=%.2f, RSF=%.1f, dt=%e, tmax=%e, collisions_on=%d\n",r->ri_hermes.hill_switch_factor,r->ri_hermes.solar_switch_factor,r->dt,tmax,coll_on);
    fprintf(out1, "\nPlanet(s):\n");
    for(int i=1;i<r->N_active;i++){
        struct reb_particle p = r->particles[i];
        fprintf(out1,"Planet %d: m=%e, r=%e, a=%e\n",i,p.m,p.r,calc_a(r,i));
    }
    fprintf(out1, "\nPlanetesimal Disk:\nNumber of planetesimals=%d\ntotal mass of planetesimal disk=%e\nmass of each planetesimal=%e\nsemi-major axis limits of planetesimal disk: a_min=%f, amax_pl=%f\npowerlaw of planetesimal disk=%.2f\n",N_planetesimals,total_disk_mass,planetesimal_mass,amin,amax,powerlaw);
    fprintf(out1, "\nMigration parameters:\nMigration Time=%f, Migration Rate=%f, K=%f",mig_time, mig_rate, K);
    fclose(out1);
    
    //run Simulations
    if(warm_start){
        // Generate Planetesimal Disk
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
        
        N_prev = r->N;
        reb_move_to_com(r);
        E0 = reb_tools_energy(r);
        
        //initial snapshot
        {char out_time[10] = {0}; sprintf(out_time,"%.0f",r->t); eia_snapshot(r, out_time);}
        
        // Integrate!
        reb_integrate(r, tmax);
    } else {
        reb_move_to_com(r);
        E0 = reb_tools_energy(r);
        
        reb_integrate(r, mig_time);
        
        //reset function pointers
        r->additional_forces = NULL;
        r->force_is_velocity_dependent = 0;
        
        // Generate Planetesimal Disk
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
        
        N_prev = r->N;
        reb_move_to_com(r);
        E0 = reb_tools_energy(r);
        rescale_energy = 1;
        
        //initial snapshot
        {char out_time[10] = {0}; sprintf(out_time,"%.0f",r->t); eia_snapshot(r, out_time);}
        
        // Integrate!
        reb_integrate(r, tmax);
    }
    
    //post integration
    {
        char out_time[10] = {0}; sprintf(out_time,"%.0f",r->t);
        eia_snapshot(r, out_time);
        char out[200] = {0}; strcat(out, binary_out); strcat(out, out_time); strcat(out, ".bin");
        reb_output_binary(r, out);
        printf("\nSimulation complete. Saved to binary \n\n");
    }
    
    free(tau_e); free(tau_a);

}

double tout = 0;
void heartbeat(struct reb_simulation* r){
    //output values to file
    if (tout <r->t){
        tout += 25;
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
    
    //output binary and eia snapshot
    if(output_time < r->t){
        char out_time[10] = {0}; sprintf(out_time,"%.0f",r->t);
        char out[200] = {0}; strcat(out, binary_out); strcat(out, out_time); strcat(out, ".bin");
        reb_output_binary(r, out);
        output_time += output_inc;
        eia_snapshot(r,out_time);
    }
    
    //Center if ejected particle
    if(r->N < N_prev){
        N_prev = r->N;
        double E = reb_tools_energy(r);
        reb_move_to_com(r);
        r->energy_offset += E - reb_tools_energy(r);
    }
    
    //output on screen
    if (reb_output_check(r, 100.*r->dt)){
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("%e",relE);
    }
    
    //rescale energy, also output post-migration eia snapshot
    if(rescale_energy == 0 && r->t >= mig_time){
        printf("\n **Migration is done, rescaling energy, outputting eia snapshot and binary**\n");
        rescale_energy = 1;
        E0 = reb_tools_energy(r);
        char out_time[10] = {0}; sprintf(out_time,"%.0f",r->t);
        eia_snapshot(r,out_time);
        char out[200] = {0}; strcat(out, binary_out); strcat(out, out_time); strcat(out, ".bin");
        reb_output_binary(r, out);
    }
}

void calc_resonant_angles(struct reb_simulation* r, FILE* f){
    double m[3] = {0};
    double e[3] = {0}; //Hardcoding, probably should change in the future.
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
    
    fprintf(f,"%e,%e,%f,%f,%f,%f,%f,%f,%f\n",m[1],m[2],a[1],e[1],a[2],e[2],phi,phi2,phi3);
    
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

void migration_forces(struct reb_simulation* r){
    if(r->t < mig_time){
        const double G = r->G;
        struct reb_particle* const particles = r->particles;
        struct reb_particle com = reb_get_com(r);
        for(int i=1;i<r->N_active;i++){
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
}

void eia_snapshot(struct reb_simulation* r, char* time){
    //name
    char dist[200] = {0}; strcat(dist,eia_out); strcat(dist,time); strcat(dist,".txt");
    FILE* append = fopen(dist,"a");
    
    //output
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = reb_get_com(r);
    double t = r->t;
    struct reb_particle p0 = particles[0];
    for(int i=1;i<r->N;i++){
        struct reb_particle p = particles[i];
        const double mu = r->G*(com.m + p.m);
        const double dvx = p.vx-com.vx;
        const double dvy = p.vy-com.vy;
        const double dvz = p.vz-com.vz;
        const double dx = p.x-com.x;
        const double dy = p.y-com.y;
        const double dz = p.z-com.z;
        const double hx = (dy*dvz - dz*dvy);
        const double hy = (dz*dvx - dx*dvz);
        const double hz = (dx*dvy - dy*dvx);
        const double h = sqrt(hx*hx + hy*hy + hz*hz);
        
        const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
        const double d = sqrt( dx*dx + dy*dy + dz*dz );    //distance
        const double dinv = 1./d;
        const double muinv = 1./mu;
        const double vr = (dx*dvx + dy*dvy + dz*dvz)*dinv;
        const double term1 = v2-mu*dinv;
        const double term2 = d*vr;
        const double ex = muinv*( term1*dx - term2*dvx );
        const double ey = muinv*( term1*dy - term2*dvy );
        const double ez = muinv*( term1*dz - term2*dvz );
        const double e = sqrt(ex*ex + ey*ey + ez*ez);   // eccentricity
        const double a = -mu/( v2 - 2.*mu*dinv );
        const double inc = acos(hz/h);
        const double rdist = sqrt((p.x-p0.x)*(p.x-p0.x)+(p.y-p0.y)*(p.y-p0.y)+(p.z-p0.z)*(p.z-p0.z));
        fprintf(append,"%f,%u,%f,%f,%f,%f,%e\n",t,p.hash,a,e,inc,rdist,p.m);
    }
    
    fclose(append);
}
