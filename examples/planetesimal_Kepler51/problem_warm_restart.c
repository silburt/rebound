/**
 * Problem Warm Restart
 *
 * This macro restarts a given warm_start simulation that had been integrated for the full 
 * simulation time, but I want to simulate it further.
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
double calc_a(struct reb_simulation* r, int index);
void calc_resonant_angles(struct reb_simulation* r, FILE* f);
void eia_snapshot(struct reb_simulation* r, char* time);

double E0, tout=0;
int N_prev;
char output_name[100] = {0};
time_t t_ini;

//resonance stuff
int inner, outer, order, jay;
double* mass; double* semi; double* ecc; double* omega; double* lambda;

//binary save
char binary_out[100] = {0};

//eia snapshot
char eia_out[300] = {0};
double output_time;
double output_inc;

int main(int argc, char* argv[]){
    char binary[300] = {0}; strcat(binary, argv[1]); strcat(binary,".bin");
    struct reb_simulation* r = reb_create_simulation_from_binary(binary);
    strcat(output_name,argv[1]); strcat(output_name,"_restart");
    
    // Standard Simulation Setup
    r->integrator	= REB_INTEGRATOR_HERMES;
    r->heartbeat	= heartbeat;
    r->testparticle_type = 1;
    
    //Important parameters
    r->ri_hermes.hill_switch_factor = 6;
    r->ri_hermes.solar_switch_factor = 20.;
    r->dt = 0.01;
    double tmax = 1e6*2*M_PI;
    
    //resonance parameters - shouldn't have to change
    inner = 2;                          //inner planet in resonance (0=star)
    outer = 3;                          //outer planet in resonance
    order = 1;                          //order of the resonance
    jay = 3;                            //j value of resonance (i.e. j*lambda_1 - (j-1)*lambda_2 - omega_1)
    
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
    double N_eia_binary_outputs = 7;
    output_inc = (tmax-r->t)/N_eia_binary_outputs;
    output_time = output_inc + r->t;
    
    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    N_prev = r->N;
    output_time = r->t + 1e5;
    tout = r->t;
    
    //resonance arrays
    mass = calloc(sizeof(double),r->N_active);
    semi = calloc(sizeof(double),r->N_active);
    ecc = calloc(sizeof(double),r->N_active);
    omega = calloc(sizeof(double),r->N_active);
    lambda = calloc(sizeof(double),r->N_active);
    
    //naming
    char binary_out[300]={0}; strcat(binary_out,output_name); strcat(binary_out,"_t=");
    char timeout[300]={0}; strcat(timeout, output_name);
    strcat(eia_out, output_name); strcat(eia_out, "_eiasnapshot_t=");
    strcat(output_name,".txt");
    
    // Integrate!
    time_t t_ini = time(NULL);
    struct tm *tmp1 = gmtime(&t_ini);
    reb_integrate(r, tmax);
    
    //post integration
    {
        char out_time[10] = {0}; sprintf(out_time,"%.0f",r->t);
        eia_snapshot(r, out_time);
        char out[200] = {0}; strcat(out, binary_out); strcat(out, out_time); strcat(out, ".bin");
        reb_output_binary(r, out);
        printf("\nSimulation complete. Saved to binary \n\n");
    }
    
    //time output
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
    //output values to file
    if (tout <r->t){
        tout *= 1.01;
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
