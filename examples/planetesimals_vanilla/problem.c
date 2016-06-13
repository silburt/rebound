//This is essentially the Kirsch example, using it as a vanilla example, need to vary various parameters and make sure the scaling relations of Armitage work.  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double calc_a(struct reb_simulation* r, int index);
double draw_ainv_powerlaw(double min, double max);
void c_momentum(struct reb_simulation* r, double* dL_ang, double* dL_lin, double Lang0, double Llin0);
double c_dcom(struct reb_simulation* r);
void eia_snapshot(struct reb_simulation* r, int output_number);

double E0, LA0, LL0, eia_snapshot_dt;
struct reb_particle comr0;
int N_prev, eia_snapshot_inc, eia_active;
char output_name[100] = {0};
char removed[200] = {0};
char* argv4;
double log_constant, tlog_output, lin_constant, tlin_output;
time_t t_ini;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    double m_earth = 0.000003003;
    //double m_neptune = 0.00005149;
    
    double powerlaw = atof(argv[1]);
    int seed = atoi(argv[2]);
    srand(seed);
    strcat(output_name,argv[3]);
    argv4 = argv[3];
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->ri_hermes.hill_switch_factor = 3;         //Hill radii
    r->ri_hermes.radius_switch_factor = 20.;          //X*radius
    r->testparticle_type = 1;
	r->heartbeat	= heartbeat;
    double tmax = 1e5 * 6.283;
    
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 11;
    reb_configure_box(r,boxsize,3,3,2);
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
    
    //does order of adding bodies matter? Yes, apparently it does says Hanno. Jacobi coords.
    
    //inner massive planet to scatter planetesimals out
    double a2=3, m2=5e-4, e2=0, inc2=reb_random_normal(0.00001);
    struct reb_particle p1 = {0};
    p1 = reb_tools_orbit_to_particle(r->G, star, m2, a2, e2, inc2, 0, 0, 0);
    p1.r = 0.000467;       //radius of Jupiter (AU)
    p1.hash = r->N;
    reb_add(r, p1);
    
    //planet 2
    double a1=5, m1=2.3*m_earth, e1=0, inc1=reb_random_normal(0.00001);
    struct reb_particle p2 = {0};
    p2 = reb_tools_orbit_to_particle(r->G, star, m1, a1, e1, inc1, 0, 0, 0);
    p2.r = 0.0000788215;       //radius of particle using 2g/cm^3 (AU)
    //p2.r = 5e-4;
    p2.hash = r->N;
    reb_add(r, p2);

    r->N_active = r->N;
    r->dt = pow(a2,1.5)/30;
    
    //planetesimals-what's a reasonable ini? Perfectly cold disk seems unlikely...
    //double planetesimal_mass = m1/600;     //each planetesimal = 1/600th of planet mass
    //int N_planetesimals = 230.*m_earth/planetesimal_mass;
    double total_disk_mass = m1*10;
    int N_planetesimals = 20000;
    double planetesimal_mass = total_disk_mass / N_planetesimals;
    //double amin = a1 - 2, amax = a1 + 2;          //planet in center of disk
    //double amin = a1 - 0.5, amax = a1 + 2.5;      //planet in asymmetric disk
    double amin = a1, amax = a1 + 3;                //planet at edge of disk
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(amin,amax,powerlaw);
        //double a = draw_ainv_powerlaw(amin,amax);
        double e = reb_random_rayleigh(0.005);   //rayleigh dist
        double inc = reb_random_rayleigh(0.005);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, r->testparticle_type?planetesimal_mass:0., a, e, inc, Omega, apsis, phi);
		pt.r 		= 0.00000934532;
        pt.hash = r->N;
		reb_add(r, pt);
    }
    
    int n_output = 25000;
    log_constant = pow(tmax + 1, 1./(n_output - 1));
    tlog_output = r->dt;
    lin_constant = tmax/n_output;
    tlin_output = r->dt;

    reb_move_to_com(r);
    comr0 = reb_get_com(r);
    double junk=0, junk2=0;
    c_momentum(r, &LA0, &LL0, junk, junk2);
    E0 = reb_tools_energy(r);
    N_prev = r->N;
    
    //eia
    eia_active = 1;
    eia_snapshot_inc = 0;
    eia_snapshot_dt = 5e3*6.283;
    
    //naming stuff
    char timeout[200] = {0};
    strcat(timeout,output_name);
    strcat(removed,output_name); strcat(removed,"_removed.txt");
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,output_name); strcat(syss,"*");
    system(syss);
    strcat(output_name,".txt");
    time_t t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    
    //Integrate!
    reb_integrate(r, tmax);
    
    time_t t_fini = time(NULL);
    struct tm *tmp2 = gmtime(&t_fini);
    double time = t_fini - t_ini;
    strcat(timeout,"_elapsedtime.txt");
    FILE* outt = fopen(timeout,"w");
    fprintf(outt,"\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    fclose(outt);
    printf("\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    
}

void heartbeat(struct reb_simulation* r){
    if(r->t > tlog_output || r->t > tlin_output){//log output or linear output!!
        if(r->t > tlog_output)tlog_output = r->t*log_constant; else tlin_output = r->t+lin_constant;
        
        double E = reb_tools_energy(r);
        double dE = (E-E0)/E0;
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
        double a1 = calc_a(r,1);
        double a2 = 0;
        if(r->N_active > 2) a2 = calc_a(r,2);
        
        //calculate actual distance of planet from star, not semi-major axis
        struct reb_particle* global = r->particles;
        double dx = global[0].x - global[1].x;
        double dy = global[0].y - global[1].y;
        double dz = global[0].z - global[1].z;
        double rdist1 = sqrt(dx*dx + dy*dy + dz*dz);
        double dx2 = global[0].x - global[2].x;
        double dy2 = global[0].y - global[2].y;
        double dz2 = global[0].z - global[2].z;
        double rdist2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
        
        int calc_mom = 1;
        double dLA = 0, dLL = 0; //angular momentum, linear momentum
        if(calc_mom) c_momentum(r, &dLA, &dLL, LA0, LL0);
        
        int N_mini = 0;
        if(r->integrator == REB_INTEGRATOR_HERMES) N_mini = r->ri_hermes.mini->N - r->ri_hermes.mini->N_active;
        
        time_t t_curr = time(NULL);
        struct tm *tmp2 = gmtime(&t_curr);
        double time = t_curr - t_ini;
        FILE *append;
        append = fopen(output_name, "a");
        fprintf(append, "%.16f,%.16f,%.16f,%d,%d,%.1f,%e,%e,%e,%f,%f,%f\n",r->t,dE,a1,r->N,N_mini,time,dLA,dLL,c_dcom(r),rdist1,rdist2,a2);
        fclose(append);
        
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
    }
    
    if(r->t > eia_snapshot_dt*eia_snapshot_inc && eia_active == 1){
        eia_snapshot(r,eia_snapshot_inc);
        if(eia_snapshot_inc == 0) printf("\n outputting eia_snapshots every %f years\n",eia_snapshot_dt);
        eia_snapshot_inc++;
    }
    
    //record collisions in mini
    if(r->N < N_prev){
        FILE* append = fopen(removed,"a");
        fprintf(append,"Collision,%.5f\n",r->t);
        fclose(append);
        
        N_prev = r->N;
    }
    
    /*
    //ejections
    {
        struct reb_particle* global = r->particles;
        const double ED2 = 100;
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
                
                //char removed[200] = {0}; strcat(removed,argv4); strcat(removed,"_removed"); strcat(removed,".txt");
                FILE* append = fopen(removed,"a");
                fprintf(append,"Ejection,%.5f\n",r->t);
                fclose(append);
                
                N_prev = r->N;
            }
        }
    }*/
}

double calc_a(struct reb_simulation* r, int index){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = reb_get_com(r);
    struct reb_particle p = particles[index]; //output planet only.
    const double mu = r->G*(com.m + p.m);
    const double dvx = p.vx-com.vx;
    const double dvy = p.vy-com.vy;
    const double dvz = p.vz-com.vz;
    const double dx = p.x-com.x;
    const double dy = p.y-com.y;
    const double dz = p.z-com.z;
    
    const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
    const double d = sqrt( dx*dx + dy*dy + dz*dz );    //distance
    const double dinv = 1./d;
    const double a = -mu/( v2 - 2.*mu*dinv );
    
    return a;
}

//returns value randomly drawn from P(x) = 1/x distribution
double draw_ainv_powerlaw(double min, double max){
    double y = reb_random_uniform(0., 1.);
    return exp(y*log(max/min) + log(min));
}

double c_dcom(struct reb_simulation* r){
    struct reb_particle com = reb_get_com(r);
    double dx = comr0.x - com.x;
    double dy = comr0.y - com.y;
    double dz = comr0.z - com.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
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

void eia_snapshot(struct reb_simulation* r, int output_number){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = reb_get_com(r);
    
    //name
    char outstr[15];
    sprintf(outstr, "%d", output_number);
    char dist[200] = {0}; strcat(dist,argv4); strcat(dist,"_ei"); strcat(dist,outstr); strcat(dist,".txt");
    FILE* append = fopen(dist,"a");
    
    double t = r->t;
    struct reb_particle p0 = particles[0];
    for(int i=1;i<r->N;i++){
        struct reb_particle p = particles[i]; //output planet only.
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
