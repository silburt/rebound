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
struct reb_particle ari_get_com(struct reb_simulation* r, int N_choice);
double calc_a(struct reb_simulation* r, int index);

double E0, tmax = 0;
char* mercury_dir; char* swifter_dir;

//temp
int N_prev;
char* argv4;
char output_name[100] = {0};
int *in_mini;
int N_CE = 0;
int L_CE = 0;

//swifter/mercury compare
void output_to_mercury_swifter(struct reb_simulation* r, double HSR, double tmax);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    tmax = atof(argv[1]);
    int N_planetesimals = atoi(argv[2]);
    int seed = atoi(argv[3]);
    strcat(output_name,argv[4]); strcat(output_name,".txt"); argv4=argv[4];
    
    int mercury_swifter_comp = 1;   //if set to 1, need argv[5] and argv[6]
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->ri_hermes.solar_switch_factor = 20.;         //X*radius
    r->testparticle_type = 1;
    r->heartbeat	= heartbeat;
    r->ri_hermes.hill_switch_factor = 3;        //Hill radii
    r->ri_hermes.adaptive_hill_switch_factor = 1;
    r->dt = 0.05;
    
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
    
    //planet
    {
        double m_neptune = 5e-5, r_neptune = 1.6e-4;
        double m_earth = 3e-6, r_earth = 0.00004258689;
        double a=1, m=m_neptune, e=0.01, inc=reb_random_normal(0.00001);
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p2.r = r_neptune;
        p2.hash = r->N;
        reb_add(r, p2);
    }
    
    r->N_active = r->N;
    
    //planetesimals
    double planetesimal_mass = 1e-8;
    double amin = 0.85, amax = 1.15;
    double powerlaw = 0;
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
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,argv[4]); strcat(syss,"*");
    system(syss);
    
    //swifter/mercury compare
    if(mercury_swifter_comp){
        mercury_dir = argv[5];
        swifter_dir = argv[6];
        output_to_mercury_swifter(r, r->ri_hermes.hill_switch_factor, tmax);
    }
    
    N_prev = r->N;
    argv4= argv[4];
    
    //in_mini
    in_mini = calloc(sizeof(int),r->N);
    
    //Integrate!
    reb_integrate(r, tmax);
    
}

void heartbeat(struct reb_simulation* r){
    if(tmax - r->t < 4*r->dt){//ending output
        
        double E = reb_tools_energy(r);
        double dE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
        
        int mini_N = 0;
        int mini_active = 0;
        if(r->integrator==REB_INTEGRATOR_HERMES){ mini_N =r->ri_hermes.mini->N; mini_active=r->ri_hermes.mini_active;}
        
        FILE *append;
        append = fopen(output_name, "a");
        fprintf(append, "%f,%e,%d,%d,%d,%f,%e\n",r->t,dE,r->N,mini_N,mini_active,calc_a(r,1),r->ri_hermes.hill_switch_factor);
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

void output_to_mercury_swifter(struct reb_simulation* r, double HSR, double tmax){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle p0 = particles[0];
    int N = r->N;
    int N_active = r->N_active;
    
    char s1[200]={0}; strcat(s1,swifter_dir); strcat(s1,"swifter_pl.in");
    char s2[200]={0}; strcat(s2,swifter_dir); strcat(s2,"param.in");
    char m1[200]={0}; strcat(m1,mercury_dir); strcat(m1,"mercury_big.in");
    char m2[200]={0}; strcat(m2,mercury_dir); strcat(m2,"mercury_small.in");
    char m3[200]={0}; strcat(m3,mercury_dir); strcat(m3,"mercury_param.in");
    char m4[200]={0}; strcat(m4,mercury_dir); strcat(m4,"mercury.inc");
    
    printf("\ns1=%s\n",s1);
    
    //Need Hill radii for swifter too.
    FILE* swifter = fopen(s1,"w");
    FILE* swifterparams = fopen(s2,"w");
    FILE* mercuryb = fopen(m1,"w");
    FILE* mercurys = fopen(m2,"w");
    FILE* mercuryparams = fopen(m3,"w");
    FILE* mercuryinc = fopen(m4,"w");
    
    //conversion options - swifter
    int alt_units = 0;
    double mass_conv = 1, vel_conv = 1, time_conv = 1;
    if(alt_units == 1){
        mass_conv = 2.959139768995959e-04;  //solar masses to this unit
        vel_conv = 0.017202424;             //converts [v] = AU/(yr/2pi) -> AU/day
        time_conv = 58.09155423;            //converts [yr/2pi] -> days
    }
    
    //swifter initial - Nbodies and sun:
    fprintf(swifter," %d\n",N);
    fprintf(swifter," 1 %.16f\n",particles[0].m*mass_conv);
    fprintf(swifter," .0 .0 .0\n");
    fprintf(swifter," .0 .0 .0\n");
    
    //SWIFTER - heliocentric coords
    for(int i=1;i<N;i++){
        struct reb_particle p = particles[i];
        double m = p.m*mass_conv;
        double rr = sqrt((p.x-p0.x)*(p.x-p0.x) + (p.y-p0.y)*(p.y-p0.y) + (p.z-p0.z)*(p.z-p0.z));
        fprintf(swifter," %d %.16f %f\n",i+1,m,rr*pow(p.m/(particles[0].m*3.),2./3.));
        fprintf(swifter," %f\n",p.r);
        fprintf(swifter," %.16f %.16f %.16f\n",p.x - p0.x, p.y - p0.y, p.z - p0.z);
        fprintf(swifter," %.16f %.16f %.16f\n",(p.vx - p0.vx)*vel_conv,(p.vy - p0.vy)*vel_conv,(p.vz - p0.vz)*vel_conv);
    }
    
    //SWIFTER - Other params (time, dt, etc.)
    fprintf(swifterparams,"! \n");
    fprintf(swifterparams,"! Parameter file for Swifter, with N=%d total bodies. \n",r->N);
    fprintf(swifterparams,"! \n! \n");
    fprintf(swifterparams,"T0             0.0E0 \n");
    fprintf(swifterparams,"TSTOP          %e        !In units where G=1\n",tmax);
    fprintf(swifterparams,"DT             %e        !In units where G=1\n",r->dt);
    fprintf(swifterparams,"PL_IN          swifter_pl.in\n");
    fprintf(swifterparams,"!TP_IN         tp.in     !Commented out for now, no test par\n");
    fprintf(swifterparams,"IN_TYPE        ASCII\n");
    fprintf(swifterparams,"ISTEP_OUT      10000        !# timesteps between outputs \n");
    fprintf(swifterparams,"BIN_OUT        out.dat\n");
    fprintf(swifterparams,"OUT_TYPE       REAL8\n");
    fprintf(swifterparams,"OUT_FORM       XV\n");
    fprintf(swifterparams,"OUT_STAT       NEW\n");
    fprintf(swifterparams,"ISTEP_DUMP     10000     !Dump parameters (incase of crash)\n");
    fprintf(swifterparams,"J2             0.0E0\n");
    fprintf(swifterparams,"J4             0.0E0\n");
    fprintf(swifterparams,"CHK_CLOSE      yes\n");
    fprintf(swifterparams,"CHK_RMIN       -1.0\n");
    fprintf(swifterparams,"CHK_RMAX       100.0\n");
    fprintf(swifterparams,"CHK_EJECT      -1.0\n");
    fprintf(swifterparams,"CHK_QMIN       -1.0\n");
    fprintf(swifterparams,"!CHK_QMIN_COORD HELIO\n");
    fprintf(swifterparams,"!CHK_QMIN_RANGE 1.0 1000.0\n");
    fprintf(swifterparams,"ENC_OUT        enc.dat\n");
    fprintf(swifterparams,"EXTRA_FORCE    no\n");
    fprintf(swifterparams,"BIG_DISCARD    yes\n");
    fprintf(swifterparams,"RHILL_PRESENT  yes\n");
    
    //mercury initial:
    //double day_zero = 2451179.5;
    double day_zero = 0;
    fprintf(mercuryb,")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n");
    fprintf(mercuryb,") Lines beginning with `)' are ignored.\n");
    fprintf(mercuryb,")---------------------------------------------------------------------\n");
    fprintf(mercuryb," style (Cartesian, Asteroidal, Cometary) = Cartesian\n");
    fprintf(mercuryb," epoch (in days) = %f\n",day_zero);
    fprintf(mercuryb,")---------------------------------------------------------------------\n");
    fprintf(mercurys,")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n");
    fprintf(mercurys,") Lines beginning with `)' are ignored.\n");
    fprintf(mercurys,")---------------------------------------------------------------------\n");
    fprintf(mercurys," style (Cartesian, Asteroidal, Cometary) = Cartesian\n");
    fprintf(mercurys,")---------------------------------------------------------------------\n");
    
    //MERCURY - heliocentric coords
    //massive planets
    double AU_d = 0.01720242383; //converts [v] = AU/(yr/2pi) -> AU/day
    for(int i=1;i<N_active;i++){
        struct reb_particle p = particles[i];
        double rr = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        fprintf(mercuryb," BD%d        m=%.16f r=%f\n",i,p.m,HSR*rr*pow(p.m/(particles[0].m*3.),2./3.));
        fprintf(mercuryb," %.16f %.16f %.16f\n",p.x - p0.x,p.y - p0.y,p.z - p0.z);
        fprintf(mercuryb," %.16f %.16f %.16f\n",(p.vx - p0.vx)*AU_d,(p.vy - p0.vy)*AU_d,(p.vz - p0.vz)*AU_d);   //AU/day
        fprintf(mercuryb," 0. 0. 0.\n");
    }
    //mini bodies
    for(int i=N_active;i<N;i++){
        struct reb_particle p = particles[i];
        double rr = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        fprintf(mercurys," BD%d        m=%.16f r=%f\n",i,p.m,HSR*rr*pow(p.m/(particles[0].m*3.),2./3.));
        fprintf(mercurys," %.16f %.16f %.16f\n",p.x - p0.x,p.y - p0.y,p.z - p0.z);     //AU, heliocentric
        fprintf(mercurys," %.16f %.16f %.16f\n",(p.vx - p0.vx)*AU_d,(p.vy - p0.vy)*AU_d,(p.vz - p0.vz)*AU_d);   //AU/day
        fprintf(mercurys," 0. 0. 0.\n");
    }
    
    //Mercury param file
    double yr2day = 1.0/AU_d; //yr/2pi -> day
    fprintf(mercuryparams,")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n");
    fprintf(mercuryparams,") Lines beginning with `)' are ignored.\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams,") Important integration parameters:\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams," algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = hyb\n");
    fprintf(mercuryparams," start time (days)= %f\n",day_zero);
    fprintf(mercuryparams," stop time (days) =%.1f\n",tmax*yr2day + day_zero);
    fprintf(mercuryparams," output interval (days) = %.1f\n",tmax*yr2day + day_zero);
    fprintf(mercuryparams," timestep (days) = %f\n",r->dt*yr2day);
    fprintf(mercuryparams," accuracy parameter=1.d-12\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams,") Integration options:\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams," stop integration after a close encounter = no\n");
    fprintf(mercuryparams," allow collisions to occur = yes\n");
    fprintf(mercuryparams," include collisional fragmentation = no\n");
    fprintf(mercuryparams," express time in days or years = years\n");
    fprintf(mercuryparams," express time relative to integration start time = yes\n");
    fprintf(mercuryparams," output precision = medium\n");
    fprintf(mercuryparams," < not used at present >\n");
    fprintf(mercuryparams," include relativity in integration= no\n");
    fprintf(mercuryparams," include user-defined force = no\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams,") These parameters do not need to be adjusted often:\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams," ejection distance (AU)= 100\n");
    fprintf(mercuryparams," radius of central body (AU) = 0.005\n");
    fprintf(mercuryparams," central mass (solar) = 1.0\n");
    fprintf(mercuryparams," central J2 = 0\n");
    fprintf(mercuryparams," central J4 = 0\n");
    fprintf(mercuryparams," central J6 = 0\n");
    fprintf(mercuryparams," < not used at present >\n");
    fprintf(mercuryparams," < not used at present >\n");
    fprintf(mercuryparams," Hybrid integrator changeover (Hill radii) = 3.\n");
    fprintf(mercuryparams," number of timesteps between data dumps = 10000\n");
    fprintf(mercuryparams," number of timesteps between periodic effects = 100\n");
    
    //Mercury.inc
    int cmax; if(N > 1e4) cmax = N/100; else cmax = 50;
    fprintf(mercuryinc,"c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
    fprintf(mercuryinc,"c\n");
    fprintf(mercuryinc,"c      MERCURY.INC    (ErikSoft   4 March 2001)\n");
    fprintf(mercuryinc,"c\n");
    fprintf(mercuryinc,"c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
    fprintf(mercuryinc,"c\n");
    fprintf(mercuryinc,"c Author: John E. Chambers\n");
    fprintf(mercuryinc,"c\n");
    fprintf(mercuryinc,"c Parameters that you may want to alter at some point:\n");
    fprintf(mercuryinc,"c\n");
    fprintf(mercuryinc,"c NMAX  = maximum number of bodies\n");
    fprintf(mercuryinc,"c CMAX  = maximum number of close-encounter minima monitored simultaneously\n");
    fprintf(mercuryinc,"c NMESS = maximum number of messages in message.in\n");
    fprintf(mercuryinc,"c HUGE  = an implausibly large number\n");
    fprintf(mercuryinc,"c NFILES = maximum number of files that can be open at the same time\n");
    fprintf(mercuryinc,"c\n");
    fprintf(mercuryinc,"      integer NMAX, CMAX, NMESS, NFILES\n");
    fprintf(mercuryinc,"      real*8 HUGE\n");
    fprintf(mercuryinc,"c\n");
    fprintf(mercuryinc,"      parameter (NMAX = %d)\n",N);
    fprintf(mercuryinc,"      parameter (CMAX = %d)\n",cmax);
    fprintf(mercuryinc,"      parameter (NMESS = 200)\n");
    fprintf(mercuryinc,"      parameter (HUGE = 9.9d29)\n");
    fprintf(mercuryinc,"      parameter (NFILES = 50)\n");
    fprintf(mercuryinc,"c\n");
    fprintf(mercuryinc,"c------------------------------------------------------------------------------\n");
    fprintf(mercuryinc,"c\n");
    fprintf(mercuryinc,"c Constants:\n");
    fprintf(mercuryinc,"c\n");
    fprintf(mercuryinc,"c DR = conversion factor from degrees to radians\n");
    fprintf(mercuryinc,"c K2 = Gaussian gravitational constant squared\n");
    fprintf(mercuryinc,"c AU = astronomical unit in cm\n");
    fprintf(mercuryinc,"c MSUN = mass of the Sun in g\n");
    fprintf(mercuryinc,"c\n");
    fprintf(mercuryinc,"      real*8 PI,TWOPI,PIBY2,DR,K2,AU,MSUN\n");
    fprintf(mercuryinc,"c\n");
    fprintf(mercuryinc,"      parameter (PI = 3.141592653589793d0)\n");
    fprintf(mercuryinc,"      parameter (TWOPI = PI * 2.d0)\n");
    fprintf(mercuryinc,"      parameter (PIBY2 = PI * .5d0)\n");
    fprintf(mercuryinc,"      parameter (DR = PI / 180.d0)\n");
    fprintf(mercuryinc,"      parameter (K2 = 2.959122082855911d-4)\n");
    fprintf(mercuryinc,"      parameter (AU = 1.4959787e13)\n");
    fprintf(mercuryinc,"      parameter (MSUN = 1.9891e33)\n");
    
    fclose(mercuryb);
    fclose(mercurys);
    fclose(swifter);
    fclose(swifterparams);
    fclose(mercuryparams);
    fclose(mercuryinc);
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
