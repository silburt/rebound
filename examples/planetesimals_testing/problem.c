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
double c_dcom(struct reb_simulation* r);
void c_orb(struct reb_simulation* r, int planet_index, double* a, double* e, double* w, double* sinf, double* cosf, double* ey);
void apsidal_precession(struct reb_simulation* r, int index, double dwdt);
void c_momentum(struct reb_simulation* r, double* dL_ang, double* dL_lin, double Lang0, double Llin0);
void nearmiss_hit_tests(struct reb_simulation* r);

double E0, LA0, LL0, t_output, t_log_output, xyz_t = 0;
int xyz_counter = 0, numdt = 20;
char* mercury_dir; char* swifter_dir;

double HSRbase;

//temp
int output_xyz = 0; //switch to 0 for no outputs
time_t t_ini;
int N_prev;
char* argv4;
char output_name[100] = {0};
int *in_mini;
int N_CE = 0;
int L_CE = 0;
struct reb_particle comr0;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    r->usleep = 1000;
    
    double tmax = INFINITY;
    //int N_planetesimals = 1000;
    int N_planetesimals = 0;
    double planet_radius_fac = atof(argv[1]);
    int seed = atoi(argv[2]);
    strcat(output_name,argv[3]); strcat(output_name,".txt"); argv4=argv[3];
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->ri_hermes.radius_switch_factor = 20.;         //X*radius
    r->testparticle_type = 1;
    r->heartbeat	= heartbeat;
    r->ri_hermes.hill_switch_factor = 6;        //Hill radii
    r->dt = 0.001;
    
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;     //switch to track the energy from collisions/ejections
    
    //For many ejected planetesimals this leads to error jumps.
    //r->boundary	= REB_BOUNDARY_OPEN;
    //const double boxsize = 10;
    //reb_configure_box(r,boxsize,1,1,1);
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
   
    srand(seed);
    t_log_output = 1.00048;
    t_output = r->dt;
    
    int planetesimal_nearmiss_hit_tests = 1;
    
    //planet 1
    {
        double e,a,m;
        if(planetesimal_nearmiss_hit_tests){a = 0.5, m = 1e-5; e=0.1;    //for collisions with planetesimals
        }else{a=0.5,m=5e-5,e=0;}
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, 0, 0, 0, 0);
        p.r = planet_radius_fac*1.6e-4;              //radius of particle is in AU!
        p.id = r->N;
        reb_add(r, p);
    }
    
    r->N_active = r->N;
    
    //planetesimals
    if(planetesimal_nearmiss_hit_tests) N_planetesimals = 0;
    double planetesimal_mass = 1e-8;
    double amin = 0.4, amax = 0.6;        //for planetesimal disk
    double powerlaw = 0.5;
    //double e_pp = 0.1;
    double e_pp = 0;
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(amin,amax,powerlaw);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        double inc = reb_random_normal(0.00001);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, e_pp, inc, Omega, apsis,phi);
		pt.r 		= 4e-5;
        pt.id = r->N;
        
        reb_add(r, pt);
    }
    
    if(planetesimal_nearmiss_hit_tests) nearmiss_hit_tests(r);
    
    //com
    reb_move_to_com(r);
    comr0 = reb_get_com(r);
    
    //energy & mom
    E0 = reb_tools_energy(r);
    double junk=0, junk2=0;
    c_momentum(r, &LA0, &LL0, junk, junk2);
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,argv[3]); strcat(syss,"*");
    system(syss);
    
    t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    N_prev = r->N;
    
    //in_mini
    in_mini = calloc(sizeof(int),r->N);
    
    //Integrate!
    reb_integrate(r, tmax);
    
}

void heartbeat(struct reb_simulation* r){
    //change hill switch radius
    
    if(r->t > t_output){//log output
        t_output = r->t*t_log_output;
        
        /*
        //update COM each iteration - testing
        int move_com = 0;
        if(move_com){
            double Eicom = reb_tools_energy(r);
            reb_move_to_com(r);
            double Efcom = reb_tools_energy(r);
            r->ri_hermes.com_dE += Eicom - Efcom;
        }*/
        
        //calc e, w
        int calc_orb = 0;
        double a=0, e=0, w=0, sinf=0, cosf=0, ey=0;
        if(calc_orb){
            int body_index = 2;
            c_orb(r,body_index,&a,&e,&w,&sinf,&cosf,&ey);
            
            int apsidal_precess = 0;            //apdidal precession
            double dwdt=0.1;
            if(apsidal_precess && r->N == 2) apsidal_precession(r, body_index, dwdt);
        }
        
        int calc_mom = 1;
        double dLA = 0, dLL = 0; //angular momentum, linear momentum
        if(calc_mom) c_momentum(r, &dLA, &dLL, LA0, LL0);
        
        int N_mini = 0;
        if(r->ri_hermes.mini_active) N_mini = r->ri_hermes.mini->N;
        
        double E = reb_tools_energy(r);// + r->ri_hermes.com_dE;
        //double dE = fabs((E-E0)/E0);
        double dE = (E-E0)/E0;
        reb_output_timing(r, 0);
        printf("    dE=%e,N_mini=%d",dE,N_mini);
        
        time_t t_curr = time(NULL);
        struct tm *tmp2 = gmtime(&t_curr);
        double time = t_curr - t_ini;
        
        FILE *append;
        append = fopen(output_name, "a");
        fprintf(append, "%.16f,%.16f,%d,%d,%.1f,%d,%e,%e,%e,%e,%e,%e\n",r->t,dE,r->N,N_mini,time,N_CE,E-E0,E,E0,dLA,dLL,c_dcom(r));
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
                reb_move_to_com(r);
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

//works only for one planet!
double c_angle(struct reb_simulation* r, double r_ps, int id, int incoming){
    struct reb_particle* p = r->particles;
    for(int i=r->N_active;i<r->N;i++){//get correct index
        if(p[i].id == id){
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

double c_dcom(struct reb_simulation* r){
    struct reb_particle com = reb_get_com(r);
    double dx = comr0.x - com.x;
    double dy = comr0.y - com.y;
    double dz = comr0.z - com.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

void c_orb(struct reb_simulation* r, int planet_index, double* a, double* e, double* w, double* sinf, double* cosf, double* eyy){
    double ee, ww;
    struct reb_particle* global = r->particles;
    struct reb_particle par = global[planet_index];
    struct reb_particle com = reb_get_com(r);
    const double m = par.m;
    const double mu = r->G*(com.m + m);
    const double muinv = 1./mu;
    
    const double dvx = par.vx-com.vx;
    const double dvy = par.vy-com.vy;
    const double dvz = par.vz-com.vz;
    const double dx = par.x-com.x;
    const double dy = par.y-com.y;
    const double dz = par.z-com.z;
    
    const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
    const double rr = sqrt ( dx*dx + dy*dy + dz*dz );
    const double vr = (dx*dvx + dy*dvy + dz*dvz)/rr;
    
    const double term1 = v*v-mu/rr;
    const double term2 = rr*vr;
    const double ex = muinv*( term1*dx - term2*dvx );
    const double ey = muinv*( term1*dy - term2*dvy );
    const double ez = muinv*( term1*dz - term2*dvz );
    double e2 = ex*ex + ey*ey + ez*ez ;
    ee = sqrt(e2); // eccentricity
    if(ey >= 0.) ww = acos(ex/ee); else ww = 2.*M_PI - acos(ex/ee);//assumes 0 inclination!!
    
    const double rdote = dx*ex + dy*ey + dz*ez;
    double cosff = rdote/(ee*rr);
    if(cosff >= 1.) cosff = 1.;
    if(cosff <= -1.) cosff = -1.;
    double sinff = sqrt(1. - cosff*cosff);
    if(vr < 0.) sinff *= -1.;
    
    double aa = rr*(1. + ee*cosff)/(1. - e2);
    
    *cosf = cosff;
    *sinf = sinff;
    *a = aa;
    *e = ee;
    *w = ww;
    *eyy = ey;
}

void apsidal_precession(struct reb_simulation* r, int index, double dwdt){
    double Ei = reb_tools_energy(r);
    struct reb_particle* p = &(r->particles[index]);
    struct reb_particle com = reb_get_com(r);
    double a, e, w, sinf, cosf, ey;
    c_orb(r,index,&a,&e,&w,&sinf,&cosf,&ey);
    w += dwdt*r->dt;
    
    const double mu = r->G*(com.m + p->m);
    double const cosw = cos(w);
    double sinw = sqrt(1. - cosw*cosw);
    if(ey < 0.) sinw *= -1;
    double const coswf = cosw*cosf - sinw*sinf;
    double const sinwf = sinw*cosf + sinf*cosw;
    
    const double r_new = a*(1. - e*e)/(1. + e*cosf);
    const double x_new = r_new*coswf + com.x;
    const double y_new = r_new*sinwf + com.y;
    double n = sqrt(mu/(a*a*a));
    
    const double term = n*a/sqrt(1.- e*e);
    const double rdot = term*e*sinf;
    const double rfdot = term*(1. + e*cosf);
    const double vx_new = rdot*coswf - rfdot*sinwf + com.vx;
    const double vy_new = rdot*sinwf + rfdot*coswf + com.vy;
    
    //Stop program if nan values being produced.
    if(x_new!=x_new || y_new!=y_new || vx_new!=vx_new ||vy_new!=vy_new){
        printf("\n\nProblem in apsidal precession calc. Exiting.\n");
        exit(0);
    }
    
    p->x = x_new;
    p->y = y_new;
    p->vx = vx_new;
    p->vy = vy_new;
    
    r->energy_offset += Ei - reb_tools_energy(r); //account for energy change because of this.
}

void c_momentum(struct reb_simulation* r, double* dL_ang, double* dL_lin, double Lang0, double Llin0){
    struct reb_particle com = reb_get_com(r);
    double Lang = 0, Llin = 0;
    for(int i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        double dx = p.x - com.x;
        double dy = p.y - com.y;
        double dz = p.z - com.z;
        double dvx = p.vx - com.vx;
        double dvy = p.vy - com.vy;
        double dvz = p.vz - com.vz;
        double hx = (dy*dvz - dz*dvy);
        double hy = (dz*dvx - dx*dvz);
        double hz = (dx*dvy - dy*dvx);
        Lang += p.m*sqrt(hx*hx + hy*hy + hz*hz);
        Llin += p.m*sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
    }
    
    if(r->t > 0){
        *dL_ang = (Lang - Lang0)/Lang0;
        *dL_lin = (Llin - Llin0)/Llin0;
    } else {
        *dL_ang = Lang; //iteration 0, calculating Lang0, Llin0
        *dL_lin = Llin;
    }
}

void nearmiss_hit_tests(struct reb_simulation* r){
    struct reb_particle star = r->particles[0];
    
    int hit = 1;
    double mass = 1e-9;
    {//planetesimal 1
        double f;
        if(hit)f = -0.94; else f = -0.91; //nearmiss/short encounter
        struct reb_particle pt = {0};
        pt = reb_tools_orbit_to_particle(r->G, star, mass, r->particles[1].x+0.1, 0.4, 0, 0, 0, f);
        pt.r = 4e-5;
        pt.id = r->N;
        reb_add(r, pt);
    }
    {//planetesimal 2
        double f;
        if(hit)f = 2.6748; else f = 2.67; //nearmiss
        struct reb_particle pt = {0};
        pt = reb_tools_orbit_to_particle(r->G, star, mass, r->particles[1].x+0.2, 0.4, 0, 0, 0, f);
        pt.r = 4e-5;
        pt.id = r->N;
        reb_add(r, pt);
    }
    
    {//planetesimal 3
        double f;
        if(hit)f = 1.828; else f = 1.83; //nearmiss
        struct reb_particle pt = {0};
        pt = reb_tools_orbit_to_particle(r->G, star, mass, r->particles[1].x+0.3, 0.3, 0, 0, M_PI, f);
        pt.r = 4e-5;
        pt.id = r->N;
        reb_add(r, pt);
    }
    
    {//planetesimal 4
        double f;
        if(hit)f= 1.9522; else f=1.951;
        struct reb_particle pt = {0};
        pt = reb_tools_orbit_to_particle(r->G, star, mass, r->particles[1].x+0.3, 0.8, 0, 0, M_PI, f);
        pt.r = 4e-5;
        pt.id = r->N;
        reb_add(r, pt);
    }
    
    /*
     {//planetesimal glide beside
     double a=0.4895, e=0.1, f=-0.27;  //long encounter
     struct reb_particle pt = {0};
     pt = reb_tools_orbit_to_particle(r->G, star, 1e-9, a, e, 0, 0, 0, f);
     pt.r = 4e-5;
     pt.id = r->N;
     reb_add(r, pt);
     }*/
}
