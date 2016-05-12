#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double E0, t_output, t_log_output;
char output_name[100] = {0};
int N_prev;
char* argv1;

int main(int argc, char* argv[]){
    strcat(output_name,argv[1]); strcat(output_name,".txt"); argv1=argv[1];
    char binary[200] = {0}; strcat(binary,argv[1]); strcat(binary,".bin");
    
    //I think the binary file is currently not carrying through hybarid
    struct reb_simulation* r = reb_create_simulation_from_binary(binary);
    double tmax = atof(argv[2]);
    E0 = atof(argv[3]);
    N_prev = r->N;
    
    int n_output = 10000;
    t_log_output = pow(tmax + 1, 1./(n_output - 1));
    t_output = r->t;
    printf("tlogoutput=%f",t_log_output);
    
    reb_integrate(r,tmax);
}

void heartbeat(struct reb_simulation* r){
    if(r->t > t_output){//log output
        t_output = r->t*t_log_output;
        
        double E = reb_tools_energy(r) + r->ri_hybarid.com_dE;
        double dE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
        
        FILE *append;
        append = fopen(output_name, "a");
        fprintf(append, "%.16f,%.16f,%d,%d,%.1f,0,0,%e,%e,%e,%e,%e,%e,%e\n",r->t,dE,r->N,r->ri_hybarid.mini->N,fabs(E-E0),E,E0,r->ri_hybarid.com_gp,r->ri_hybarid.com_gf,r->ri_hybarid.com_mp,r->ri_hybarid.com_mf,r->ri_hybarid.com_dE);
        fclose(append);
    }
    
    //record collisions in mini
    if(r->N < N_prev){
        char removed[200] = {0}; strcat(removed,argv1); strcat(removed,"_removed.txt");
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
                r->collisions_dE += Ei - Ef;
                
                char removed[200] = {0}; strcat(removed,argv1); strcat(removed,"_removed"); strcat(removed,".txt");
                FILE* append = fopen(removed,"a");
                fprintf(append,"Ejection,%.5f\n",r->t);
                fclose(append);
                
                N_prev = r->N;
            }
        }
    }

}

