/**
 * Forward Migration of Jupiter
 *
 * This example integrates Jupiter's forward migration with an interior planetesimal disk.
 * It migrates forward in the disk due to gas, and the planetesimals it receives in its wake 
 * increases its eccentricity. After Saturn and Jupiter encounter and pass through 
 * the 2:1 MMR, Saturn is at a good eccentricity. Let's also check if Earth and Mars 
 * are okay too as Jupiter passes through secular resonances. The paper tells us the 
 * rates it has to be!
 *
 * Jupiter is at the inner edge of the planetesimal disk, and the goal is for Jupiter to send planetesimals towards Earth and Mars. 
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
void migration_forces(struct reb_simulation* r);
void calc_resonant_angles(struct reb_simulation* r, FILE* f);

double calc_a(struct reb_simulation* r, int index);

double E0;
char output_name[100] = {0};
time_t t_ini;
double* tau_a; double* tau_e; double* tau_m;    /**< Migration, e-damp, and mass gain timescales */
double* omega; double* lambda; double* a; double* e; double* m;
double m_earth = 3e-6;
int iJ;         /**< The index of Jupiter */
int N_prev;     //keep track of ejected particles
double tmax;

//double m_mars = 0.107*m_earth, m_jupiter = 317.8*m_earth, m_saturn = 95.16*m_earth;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    //Parameters
    double a_superearth = atof(argv[1]);
    double K = atof(argv[2]);       //Lee & Peale (2002) K.
    strcat(output_name,argv[3]);
    
    //other parameters
    double growthfac = 1; if(growthfac < 1)growthfac = 1;
    iJ = 4;     //index of Jupiter, this is important!
    
    r->integrator	= REB_INTEGRATOR_HERMES;
    
	// Simulation Setup
    r->heartbeat	= heartbeat;
    r->additional_forces = migration_forces;
    r->force_is_velocity_dependent = 1;
    r->testparticle_type = 1;
    r->dt = 0.125;
    tmax = 6e6;
    
    // Collisions
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
    // Boundaries
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 35;
    reb_configure_box(r,boxsize,2,2,1);
    
    int seed = 11;
    srand(seed);
    
	// Star
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
    star.hash = 0;
	reb_add(r, star);
    
    //Earth
    double m_earth = 3e-6, r_earth=4.26349651e-5, a_uranus=19.18;
    {
        double a=1, m=m_earth, e=0, inc = reb_random_normal(0.00001);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, reb_random_uniform(0,2*M_PI));
        p1.r = r_earth;
        reb_add(r, p1);
    }
    
    //Mars
    {
        double a=1.5273, m=0.107*m_earth, e=0, inc = reb_random_normal(0.00001);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, reb_random_uniform(0,2*M_PI));
        p1.r = 2.27075425e-5;
        reb_add(r, p1);
    }

    //Super Earth - tends to get ejected
    {
        //double a_default = 3.7;
        double a=a_superearth, m=5*m_earth, e=0, inc = reb_random_normal(0.00001);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, reb_random_uniform(0,2*M_PI));
        p1.r = 2.27075425e-5;
        reb_add(r, p1);
    }
    
    //Jupiter
    {
        double a_real = 5.2, m_real = 0.0009543;
        double a=a_uranus-10, m=0.0009543/growthfac, e=0, inc = reb_random_normal(0.00001);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, reb_random_uniform(0,2*M_PI));
        p1.r = 0.00046732617;
        reb_add(r, p1);
    }
    
    //Saturn
    {
        double a_real = 9.54;
        double a=a_uranus-6, m=0.0002857/growthfac, e=0, inc=reb_random_normal(0.00001);
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, reb_random_uniform(0,2*M_PI));
        p2.r = 0.000389256877;
        reb_add(r, p2);
    }
    
    //Uranus
    {
        double a=a_uranus-2, m=0.00004511, e=0.0, inc=reb_random_normal(0.00001);
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, reb_random_uniform(0,2*M_PI));
        p2.r = 0.000164939;
        reb_add(r, p2);
    }
    
    //Neptune
    {
        double a_real = 30.06;
        double a=a_uranus+5, m=0.0000511, e=0.0, inc=reb_random_normal(0.00001);
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, reb_random_uniform(0,2*M_PI));
        p2.r = 0.000159835;
        reb_add(r, p2);
    }
    
    r->N_active = r->N;
    
    //Planetesimal disk parameters
    int N_planetesimals = 0;
    double planetesimal_mass = 0.012*m_earth;    //moon sized bodies
    double amin = calc_a(r, 2), amax = calc_a(r, 3) + 0.5;
    double powerlaw = 1;
    
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
        pt.r 		= 0.27*r_earth; //moon sized body
        pt.hash = r->N;
        reb_add(r, pt);
    }
    
    //Migration Stuff
    //parameters
    double mig_rate = 2.*M_PI * 5e5;
    //arrays
    tau_a = calloc(sizeof(double),r->N);
    tau_e = calloc(sizeof(double),r->N);
    tau_m = calloc(sizeof(double),r->N);                //mass doubling every mig_rate years
    tau_a[iJ] = mig_rate;                               //Jupiter  - mass growth rate too!
    tau_e[iJ] = mig_rate/K;
    tau_m[iJ] = pow(growthfac, (r->dt/tau_a[iJ]));
    tau_a[iJ+1] = 1.5*mig_rate;                        //Saturn - mass growth rate too!
    tau_e[iJ+1] = mig_rate/K;
    tau_m[iJ+1] = pow(growthfac, (r->dt/tau_a[iJ]));
    tau_e[iJ+2] = mig_rate/K;                           //Uranus - eccentricity damping only
    tau_a[iJ+3] = -5*mig_rate;                         //Neptune
    tau_e[iJ+3] = 5*mig_rate/K;
 
    //extra arrays
    a = calloc(sizeof(double),r->N);
    e = calloc(sizeof(double),r->N);
    m = calloc(sizeof(double),r->N);
    omega = calloc(sizeof(double),r->N);
    lambda = calloc(sizeof(double),r->N);
    
    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    N_prev = r->N;
    
    //naming
    if(growthfac > 1) strcat(output_name,"_massgrowth");
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,output_name); strcat(syss,"*");
    char info[150] = {0}; strcat(info,output_name);
    char binary[150] = {0}; strcat(binary,output_name); strcat(binary,".bin");
    system(syss);
    strcat(output_name,".txt");
    
    //info output
    strcat(info,"_info.txt");
    FILE* out1 = fopen(info,"w");
    fprintf(out1, "Simulation Details:\n");
    fprintf(out1, "\nSetup Parmaeters:\nHSF=%.2f, SSF=%.1f, dt=%e, tmax=%e d\n",r->ri_hermes.hill_switch_factor,r->ri_hermes.solar_switch_factor,r->dt,tmax);
    fprintf(out1, "\nPlanet(s):\n");
    for(int i=1;i<r->N_active;i++){
        struct reb_particle p = r->particles[i];
        fprintf(out1,"Planet %d: m=%e, Rp=%e, a=%e\n",i,p.m,p.r,calc_a(r,i));
    }
    fprintf(out1, "\nPlanetesimal Disk:\nNumber of planetesimals=%d\ntotal mass of planetesimal disk=%e\nmass of each planetesimal=%e\nsemi-major axis limits of planetesimal disk: a_min=%f, amax_pl=%f\npowerlaw of planetesimal=%.2f\n",N_planetesimals,planetesimal_mass*N_planetesimals,planetesimal_mass,amin,amax,powerlaw);
    int tau_e_pl = 0; if(tau_e[r->N_active]!=0) tau_e_pl=1;
    fprintf(out1, "\nMigration parameters:\n Migration Rate=2pi*%f, K=%f, Jup/Sat_mass_growth=%d\n", mig_rate/(2*M_PI), K, growthfac>1?1:0);
    fclose(out1);
    
    // Integrate!
    reb_integrate(r, tmax);
    
    reb_output_binary(r, binary);
    printf("\nSimulation complete. Saved to binary \n\n");
}

double tout = 0;
void heartbeat(struct reb_simulation* r){
    if (tout <r->t){
        tout += 200;
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        int N_mini = 0;
        if (r->ri_hermes.mini_active){
            N_mini = r->ri_hermes.mini->N;
        }
        
        FILE* f = fopen(output_name, "a");
        fprintf(f,"%e,%e,%d,%d",r->t,relE,r->N,N_mini);
        calc_resonant_angles(r,f);
        fclose(f);
    }
    
    if (reb_output_check(r, 100.*r->dt)){
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("%e",relE);
    }
    
    //Center if ejected particle
    if(r->N < N_prev){
        r->t = tmax;
    }
    
    //mass growth
    for(int i=0;i<r->N_active;i++){
        struct reb_particle* p = &(r->particles[i]);
        if (tau_m[i]!=0){       //mass doubling every mig_rate years
            p->m *= tau_m[i];
        }
    }

}

void calc_resonant_angles(struct reb_simulation* r, FILE* f){
    struct reb_particle com = reb_get_com(r);
    for(int i=1;i<r->N;i++){
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
        m[i] = p.m;
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
    int ip = 3, op = 4; //inner and outer resonant planets - Juptier and Saturn
    
    double phi = 0, phi2 = 0, phi3 = 0;
    phi = 2.*lambda[op] - lambda[ip] - omega[ip];
    phi2 = 2.*lambda[op] - lambda[ip] - omega[op];
    phi3 = omega[ip] - omega[op];
    while(phi >= 2*M_PI) phi -= 2*M_PI;
    while(phi < 0.) phi += 2*M_PI;
    while(phi2 >= 2*M_PI) phi2 -= 2*M_PI;
    while(phi2 < 0.) phi2 += 2*M_PI;
    while(phi3 >= 2*M_PI) phi3 -= 2*M_PI;
    while(phi3 < 0.) phi3 += 2*M_PI;
    
    for(int i=1;i<r->N_active;i++) fprintf(f,",%e,%e,%e",a[i],e[i],m[i]);
    fprintf(f,",%f,%f,%f\n",phi,phi2,phi3);
    
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
    double ajup_curr = calc_a(r,iJ);
    double ajup_real = 5.2;
    if(ajup_curr <= ajup_real){
        for(int i=0;i<r->N;i++){
            tau_a[i] = 0;
            tau_e[i] = 0;
            tau_m[i] = 0;
        }
    }
    
    if(tau_a[iJ]>0){
        const double G = r->G;
        const int N = r->N;
        struct reb_particle* const particles = r->particles;
        struct reb_particle com = reb_get_com(r);
        for(int i=1;i<N;i++){
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