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

double E0, t_output, t_log_output, xyz_t = 0;
int xyz_counter = 0, numdt = 20;
char* mercury_dir; char* swifter_dir;

//temp
int output_xyz = 0; //switch to 0 for no outputs
time_t t_ini;
int N_prev;
char* argv4;
char output_name[100] = {0};
int *in_mini;
int N_CE = 0;
int L_CE = 0;

//swifter/mercury compare
void output_to_mercury_swifter(struct reb_simulation* r, double HSR, double tmax, int n_output);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    double tmax = atof(argv[1]);
    int N_planetesimals = atoi(argv[2]);
    int seed = atoi(argv[3]);
    strcat(output_name,argv[4]); strcat(output_name,".txt"); argv4=argv[4];
    
    int mercury_swifter_comp = 1;   //if set to 1, need argv[5] and argv[6]
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->ri_hermes.radius_switch_factor = 20.;         //X*radius
    r->testparticle_type = 1;
    r->heartbeat	= heartbeat;
    r->ri_hermes.hill_switch_factor = 3;        //Hill radii
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
        p2.r = 2e-4;
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
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,argv[4]); strcat(syss,"*");
    system(syss);
    
    //swifter/mercury compare
    if(mercury_swifter_comp){
        mercury_dir = argv[5];
        swifter_dir = argv[6];
        output_to_mercury_swifter(r, r->ri_hermes.hill_switch_factor, tmax, n_output);
    }
    
    t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    N_prev = r->N;
    argv4= argv[4];
    
    //in_mini
    in_mini = calloc(sizeof(int),r->N);
    
    //Integrate!
    reb_integrate(r, tmax);
    
    //elapsed time stuff
    time_t t_fini = time(NULL);
    struct tm *tmp2 = gmtime(&t_fini);
    double time = t_fini - t_ini;
    char timeout[200] = {0};
    strcat(timeout,argv[4]); strcat(timeout,"_elapsedtime"); strcat(timeout,".txt");
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
        
        /*
        //
        //counting close encounters + calculate angle
        int count_CE = 0;
        if(count_CE){
            struct reb_simulation* mini = r->ri_hermes.mini;
            double r_ps = 0, angle = 0;
            char ia[200] = {0}; strcat(ia,argv4); strcat(ia,"_inangle.txt");
            char oa[200] = {0}; strcat(oa,argv4); strcat(oa,"_outangle.txt");
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
        fprintf(append, "%.16f,%.16f,%d,%d,%.1f,%d,%e,%e,%e\n",r->t,dE,r->N,r->ri_hermes.mini->N,time,N_CE,fabs(E-E0),E,E0);
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

void output_to_mercury_swifter(struct reb_simulation* r, double HSR, double tmax, int n_output){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle p0 = particles[0];
    int N = r->N;
    int N_active = r->N_active;
    
    char s1[200]={0}; strcat(s1,swifter_dir); strcat(s1,"swifter_pl.in");
    char s2[200]={0}; strcat(s2,swifter_dir); strcat(s2,"param.in");
    char m1[200]={0}; strcat(m1,mercury_dir); strcat(m1,"mercury_big.in");
    char m2[200]={0}; strcat(m2,mercury_dir); strcat(m2,"mercury_small.in");
    char m3[200]={0}; strcat(m3,mercury_dir); strcat(m3,"mercury_param.in");
    
    printf("\ns1=%s\n",s1);
    
    //Need Hill radii for swifter too.
    FILE* swifter = fopen(s1,"w");
    FILE* swifterparams = fopen(s2,"w");
    FILE* mercuryb = fopen(m1,"w");
    FILE* mercurys = fopen(m2,"w");
    FILE* mercuryparams = fopen(m3,"w");
    
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
    int output_rate = round(tmax/r->dt/n_output);
    if(output_rate == 0) output_rate = 1;
    fprintf(swifterparams,"! \n");
    fprintf(swifterparams,"! Parameter file for Swifter, with N=%d total bodies. \n",r->N);
    fprintf(swifterparams,"! \n! \n");
    fprintf(swifterparams,"T0             0.0E0 \n");
    fprintf(swifterparams,"TSTOP          %e        !In units where G=1\n",tmax);
    fprintf(swifterparams,"DT             %e        !In units where G=1\n",r->dt);
    fprintf(swifterparams,"PL_IN          swifter_pl.in\n");
    fprintf(swifterparams,"!TP_IN         tp.in     !Commented out for now, no test par\n");
    fprintf(swifterparams,"IN_TYPE        ASCII\n");
    fprintf(swifterparams,"ISTEP_OUT      %d        !# timesteps between outputs \n",output_rate);
    fprintf(swifterparams,"BIN_OUT        out.dat\n");
    fprintf(swifterparams,"OUT_TYPE       REAL8\n");
    fprintf(swifterparams,"OUT_FORM       XV\n");
    fprintf(swifterparams,"OUT_STAT       NEW\n");
    fprintf(swifterparams,"ISTEP_DUMP     %d     !Dump parameters (incase of crash)\n",output_rate*500);
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
    double day_zero = 2451179.5;
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
        fprintf(mercuryb," BODY%d      m=%.16f r=%f\n",i,p.m,HSR*rr*pow(p.m/(particles[0].m*3.),2./3.));
        fprintf(mercuryb," %.16f %.16f %.16f\n",p.x - p0.x,p.y - p0.y,p.z - p0.z);
        fprintf(mercuryb," %.16f %.16f %.16f\n",(p.vx - p0.vx)*AU_d,(p.vy - p0.vy)*AU_d,(p.vz - p0.vz)*AU_d);   //AU/day
        fprintf(mercuryb," 0. 0. 0.\n");
    }
    //mini bodies
    for(int i=N_active;i<N;i++){
        struct reb_particle p = particles[i];
        double rr = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        fprintf(mercurys," BODY%d      m=%.16f r=%f\n",i,p.m,HSR*rr*pow(p.m/(particles[0].m*3.),2./3.));
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
    fprintf(mercuryparams," output interval (days) = %.2fd0\n",(tmax/r->dt/n_output)*yr2day);
    fprintf(mercuryparams," timestep (days) = %f\n",r->dt*yr2day);
    fprintf(mercuryparams," accuracy parameter=1.d-12\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams,") Integration options:\n");
    fprintf(mercuryparams,")---------------------------------------------------------------------\n");
    fprintf(mercuryparams," stop integration after a close encounter = no\n");
    fprintf(mercuryparams," allow collisions to occur = yes\n");
    fprintf(mercuryparams," include collisional fragmentation = no\n");
    fprintf(mercuryparams," express time in days or years = years\n");
    fprintf(mercuryparams," express time relative to integration start time = no\n");
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
    fprintf(mercuryparams," number of timesteps between data dumps = 500\n");
    fprintf(mercuryparams," number of timesteps between periodic effects = 100\n");
    
    fclose(mercuryb);
    fclose(mercurys);
    fclose(swifter);
    fclose(swifterparams);
    fclose(mercuryparams);
}
