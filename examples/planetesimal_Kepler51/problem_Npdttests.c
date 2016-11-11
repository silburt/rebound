/**
 * Planets in resonance embedded in a Planetesimal Disk
 *
 * In this example the Kepler-58 planets are migrated into resonance with planetesimals. The option for
 * a warm or cold start is possible. If warm start, there is one simulation with the planets 
 * migrating simultaneously with the planetesimal disk already added. If cold start, I've automated
 * the code such that first the planets migrate into resonance, the simulation is stopped, 
 * saved as a binary, then planetesimals are added, then the simulation starts again.
 *
 * A bit about Kepler-51:
 * b: P = 45.15, a = 0.2514, M = 2.1M_Earth, e=0.04, r=7.1R_Earth
 * c: P = 85.31, a = 0.384, M = 4.0M_Earth, e=0.014, r=9.0R_Earth
 * d: P = 130.19, a = 0.509, M = 7.6M_Earth, e=0.008, r=9.7R_earth
 *
 * In this problem, I'm specifically wanting to track the period evolution vs. dt and Np. 
 * Assume a cold start for simplicity. Also, making the planetesimal disk limits much closer than
 * what I think Chatterjee is setting, and removing the inner planet, since this is just a 
 * convergence test.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "integrator_hermes.c"
#include <time.h>

void heartbeat(struct reb_simulation* r);
void migration_forces(struct reb_simulation* r);
void calc_resonant_angles(struct reb_simulation* r, FILE* f);
double calc_a(struct reb_simulation* r, int index);
void eia_snapshot(struct reb_simulation* r, char* time);

//General
double E0;
char output_name[200] = {0};
time_t t_ini;
int N_prev;     //keep track of ejected particles

//migration
double* tau_a; 	/**< Migration timescale in years for all particles */
double* tau_e; 	/**< Eccentricity damping timescale in years for all particles */
double mig_time, dispersal_time, mig_rate;
int rescale_energy;

//eia snapshot
char eia_out[100] = {0};
double output_time;
double output_inc;

//binary save
char binary_out[100] = {0};
int orbital_output_freq;

//resonance stuff
int inner, outer, order, jay;
double* mass; double* semi; double* ecc; double* omega; double* lambda;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    double m_earth = 0.000003003;       //convert Earth mass to Solar mass.
    double r_earth = 4.2587e-5;         //convert Earth radii to AU.
    
//**********************Parameter List**********************
    int warm_start = 0;     //if 0, cold start, if 1, warm start
    double dt = atof(argv[1])*2*M_PI;
    int N_planetesimals = atoi(argv[2]);
    int seed = atoi(argv[3]);
    strcat(output_name,argv[4]);        //name
    
    //Migration parameters
    mig_time = 35000;
    dispersal_time = 2000;              //gas is dispersing, mig_rate -> infinity.
    mig_rate = 5e5*2*M_PI;              //AU/(yr/2pi)
    double K = 100.0;                   //Lee & Peale (2002) K.
    double da = 0.00;                   //fractional offset of outer planet from current position
    
    //Planetesimal disk parameters
    double total_disk_mass=2.7*m_earth; //2.7m_earth = 20% mass of planets, 1.16m_earth = 10% of outer planets
    double alpha = 0;                   //reb_powerlaw uses x^alpha-http://mathworld.wolfram.com/RandomNumber.html
    double N_eia_binary_outputs = 7;    //number of eia_snapshot and binary output intervals
//**********************************************************
    
    //resonance parameters - shouldn't have to change
    inner = 1;                          //inner planet in resonance (0=star)
    outer = 2;                          //outer planet in resonance
    order = 1;                          //order of the resonance
    jay = 3;                            //j value of resonance (i.e. j*lambda_1 - (j-1)*lambda_2 - omega_1)
    
	// Standard Simulation Setup
    r->integrator	= REB_INTEGRATOR_HERMES;
    r->heartbeat	= heartbeat;
    r->additional_forces = migration_forces;
    r->force_is_velocity_dependent = 1;
    r->testparticle_type = 1;
    
    //****Testing**** simulation params
    r->ri_hermes.adaptive_hill_switch_factor = 0;
    
    //Important parameters
    r->ri_hermes.hill_switch_factor = 6;
    r->ri_hermes.solar_switch_factor = 20.;
    r->dt = dt;
    double tmax = 1e5*2*M_PI + mig_time + dispersal_time;
    
    // Collisions
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
    // Boundaries
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 1.5;
    reb_configure_box(r,boxsize,2,2,1);
    
    //output stuff
    output_inc = tmax/N_eia_binary_outputs;
    output_time = output_inc;
    orbital_output_freq = 25;       //how often to output orbital parameters to *.txt
//**********************************************************
    
    srand(seed);
    
	// Star
	struct reb_particle star = {0};
	star.m 		= 1.04;
    star.r		= 0.005;            // Radius of particle is in AU!
    star.hash = 0;
	reb_add(r, star);
    
    /*
    //b
    {
        double a_real = 0.2514; //<-true e, but start with 0!
        double a=a_real, m=2.1*m_earth, e=0, inc=reb_random_normal(0.0001);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, reb_random_uniform(0,2*M_PI));
        p.r = 7.1*r_earth;
        p.hash = r->N;
        reb_add(r, p);
    }*/
    
    //c
    {
        double a_real = 0.384; //<-true e, but start with 0!
        double a=a_real, m=4.0*m_earth, e=0, inc=reb_random_normal(0.0001);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, reb_random_uniform(0,2*M_PI));
        p.r = 9.0*r_earth;
        p.hash = r->N;
        reb_add(r, p);
    }
    
    //d
    {
        double a_real = 0.509;
        double a=pow(pow(calc_a(r,inner),1.5)*1.51,2./3.), m=7.6*m_earth, e=0, inc=reb_random_normal(0.0001);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, reb_random_uniform(0,2*M_PI));
        p.r = 9.7*r_earth;
        p.hash = r->N;
        reb_add(r, p);
    }
    
    r->N_active = r->N;
    //following Chatterjee, inner/outer edges are 0.01AU inside 3:1 period ratio for inner/outer planet
    double amin=pow(pow(calc_a(r,1),1.5)/1.5,2./3.) + 0.01;
    double amax=pow(pow(calc_a(r,outer)-da*calc_a(r,outer),1.5)*1.5,2./3.) - 0.01; //account for da offset
    double planetesimal_mass = total_disk_mass/N_planetesimals;
    
    //resonance arrays
    mass = calloc(sizeof(double),r->N_active);
    semi = calloc(sizeof(double),r->N_active);
    ecc = calloc(sizeof(double),r->N_active);
    omega = calloc(sizeof(double),r->N_active);
    lambda = calloc(sizeof(double),r->N_active);
    
    //naming
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
    fprintf(out1, "\nSetup Parmaeters:\nHSF=%.2f, RSF=%.1f, dt=%e, tmax=%e, collisions_on=%d, AdaptiveHSF=%d\n",r->ri_hermes.hill_switch_factor,r->ri_hermes.solar_switch_factor,r->dt,tmax,coll_on,r->ri_hermes.adaptive_hill_switch_factor);
    fprintf(out1, "\nPlanet(s):\n");
    for(int i=1;i<r->N_active;i++){
        struct reb_particle p = r->particles[i];
        fprintf(out1,"Planet %d: m=%e, r=%e, a=%e\n",i,p.m,p.r,calc_a(r,i));
    }
    fprintf(out1, "\nPlanetesimal Disk:\nNumber of planetesimals=%d\ntotal mass of planetesimal disk=%e\nmass of each planetesimal=%e\nsemi-major axis limits of planetesimal disk: a_min=%f, amax_pl=%f\npowerlaw of planetesimal disk=%.2f\n",N_planetesimals,total_disk_mass,planetesimal_mass,amin,amax,alpha);
    fprintf(out1, "\nMigration parameters:\nMigration Time=%f, Migration Rate=%f, K=%f, Dispersal_time=%f",mig_time, mig_rate, K, dispersal_time);
    fclose(out1);
    
    //run Simulations
    if(warm_start){//*****************************WARM START*****************************
        // Generate Planetesimal Disk
        while(r->N<(N_planetesimals + r->N_active)){
            struct reb_particle pt = {0};
            double a    = reb_random_powerlaw(amin,amax,alpha);
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
        
        //damping (migration) arrays
        tau_a = calloc(sizeof(double),r->N);
        tau_e = calloc(sizeof(double),r->N);
        //tau_a[inner] = mig_rate;
        //tau_e[inner] = tau_a[inner]/K;
        tau_a[outer] = mig_rate;
        tau_e[outer] = tau_a[outer]/K;
        for(int i=r->N_active;i<r->N;i++){tau_a[i] = 2*tau_a[outer]; tau_e[i] = 2*tau_e[i]/K;}   //planetesimal damping
        
        N_prev = r->N;
        reb_move_to_com(r);
        E0 = reb_tools_energy(r);
        
        //initial snapshot
        {char out_time[10] = {0}; sprintf(out_time,"%.0f",r->t); eia_snapshot(r, out_time);}
        
        //timing
        t_ini = time(NULL);
        struct tm *tmp = gmtime(&t_ini);
        
        // Integrate!
        reb_integrate(r, tmax);
    } else {//*****************************COLD START*****************************
        reb_move_to_com(r);
        E0 = reb_tools_energy(r);
        
        tau_a = calloc(sizeof(double),r->N);
        tau_e = calloc(sizeof(double),r->N);
        tau_a[outer] = mig_rate;
        tau_e[outer] = tau_a[outer]/K;
        
        reb_integrate(r, mig_time+dispersal_time);
        reb_output_binary(r, binary_out);
        r->integrator	= REB_INTEGRATOR_HERMES;
        r->heartbeat	= heartbeat;
        
        // Generate Planetesimal Disk
        while(r->N<(N_planetesimals + r->N_active)){
            struct reb_particle pt = {0};
            double a    = reb_random_powerlaw(amin,amax,alpha);
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
        
        //initial snapshot - for cold start, it's convenient to label the first output as t=0...
        {char out_time[10] = {0}; sprintf(out_time,"%.0f",r->t-mig_time-dispersal_time); eia_snapshot(r, out_time);}
        
        //timing
        t_ini = time(NULL);
        struct tm *tmp = gmtime(&t_ini);
        
        // Integrate!
        reb_integrate(r, tmax);
    }//*****************************POST INTEGRATION*****************************
    
    {
        char out_time[10] = {0}; sprintf(out_time,"%.0f",r->t);
        eia_snapshot(r, out_time);
        char out[200] = {0}; strcat(out, binary_out); strcat(out, out_time); strcat(out, ".bin");
        reb_output_binary(r, out);
        
        //timing
        time_t t_fini = time(NULL);
        struct tm *tmp2 = gmtime(&t_fini);
        double time = t_fini - t_ini;
        FILE* out1 = fopen(info,"a");
        fprintf(out1, "\n\n\nElapsed simulation time(s)= %f",time);
        fclose(out1);
        printf("\nSimulation complete. Saved to binary \n\n");
    }
    
    free(tau_e); free(tau_a); free(mass); free(semi); free(ecc); free(omega); free(lambda);

}

double tout = 0;
void heartbeat(struct reb_simulation* r){
    //output values to file
    if (tout <r->t){
        tout += orbital_output_freq;
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
    if(rescale_energy == 0 && r->t >= (mig_time + dispersal_time)){
        printf("\n **Migration is done, rescaling energy, outputting eia snapshot and binary**\n");
        rescale_energy = 1;
        E0 = reb_tools_energy(r);
        char out_time[10] = {0}; sprintf(out_time,"%.0f",r->t);
        eia_snapshot(r,out_time);
        char out[200] = {0}; strcat(out, binary_out); strcat(out, out_time); strcat(out, ".bin");
        reb_output_binary(r, out);
        //orbital_output_freq = 100;  //make txt outputs much longer now.
    }
}

void calc_resonant_angles(struct reb_simulation* r, FILE* f){
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
        mass[i] = p.m;
        ecc[i] = sqrt( ex*ex + ey*ey + ez*ez );   // eccentricity
        semi[i] = -mu/(v*v - 2.*mu/d);
        const double rdote = dx*ex + dy*ey + dz*ez;
        const double cosf = rdote/(ecc[i]*d);
        
        omega[i] = atan2(ey,ex);
        if(ey < 0.) omega[i] += 2*M_PI;
        double cosE = (semi[i] - d)/(semi[i]*ecc[i]);
        double E;
        if(cosf > 1. || cosf < -1.){
            E = M_PI - M_PI*cosE;
        } else {
            E = acos(cosE);
        }
        if(vr < 0.) E = 2.*M_PI - E;
        double MA = E - ecc[i]*sin(E);
        lambda[i] = MA + omega[i];
    }
    double phi = 0, phi2 = 0, phi3 = 0;
    phi = jay*lambda[outer] - (jay-order)*lambda[inner] - omega[inner];
    phi2 = jay*lambda[outer] - (jay-order)*lambda[inner] - omega[outer];
    phi3 = omega[inner] - omega[outer];
    while(phi >= 2*M_PI) phi -= 2*M_PI;
    while(phi < 0.) phi += 2*M_PI;
    while(phi2 >= 2*M_PI) phi2 -= 2*M_PI;
    while(phi2 < 0.) phi2 += 2*M_PI;
    while(phi3 >= 2*M_PI) phi3 -= 2*M_PI;
    while(phi3 < 0.) phi3 += 2*M_PI;
    
    for(int i=1;i<r->N_active;i++) fprintf(f,"%e,%f,%f,",mass[i],semi[i],ecc[i]);
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

void migration_forces(struct reb_simulation* r){
    if(r->t < (mig_time + dispersal_time)){
        if(r->t > mig_time){
            double fac = pow(1e7/mig_rate + 1, 1./(dispersal_time/r->dt - 1));
            tau_a[outer] *= fac;
            tau_e[outer] *= fac;
        }
        const double G = r->G;
        struct reb_particle* const particles = r->particles;
        struct reb_particle com = reb_get_com(r);
        for(int i=1;i<r->N;i++){
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
                    const double a = -mu/( v*v - 2.*mu/r );             // semi major axis
                    const double prefac1 = 1./(1.-e*e) /tau_e[i]/1.5;
                    const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))  /tau_e[i]/1.5;
                    p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
                    p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
                    p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;
                }
            }
            com = reb_get_com_of_pair(com,particles[i]);
        }
        reb_move_to_com(r);
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
