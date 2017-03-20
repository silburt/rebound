/**
 * @file    transformations.c
 * @brief   Transformations back and forth between different coordinate systems.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details This file collects all the transformations used by the different integrators
 * between different coordinate systems.
 *
 * @section     LICENSE
 * Copyright (c) 2017 Hanno Rein, Dan Tamayo.
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "transformations.h"
#include "rebound.h"

/******************************
 * Jacobi */

void reb_transformations_calculate_jacobi_eta(const struct reb_particle* const ps, double* const eta, const int N){
    eta[0] = ps[0].m;
    for (unsigned int i=1;i<N;i++){
        eta[i] = eta[i-1] + ps[i].m;
    }
}

void reb_transformations_calculate_jacobi_masses(const struct reb_particle* const ps, double* const m_j, double* const eta, const int N){
    eta[0] = ps[0].m;
    for (unsigned int i=1;i<N;i++){
        eta[i] = eta[i-1] + ps[i].m;
        m_j[i] = ps[i].m*eta[i-1]/eta[i];
    }
    m_j[0] = eta[N-1];
}

void reb_transformations_inertial_to_jacobi_posvel(const struct reb_particle* const particles, struct reb_particle* const p_j, const double* const eta, const struct reb_particle* const p_mass, const int N){
    double s_x = eta[0] * particles[0].x;
    double s_y = eta[0] * particles[0].y;
    double s_z = eta[0] * particles[0].z;
    double s_vx = eta[0] * particles[0].vx;
    double s_vy = eta[0] * particles[0].vy;
    double s_vz = eta[0] * particles[0].vz;
    for (unsigned int i=1;i<N;i++){
        const double ei = 1./eta[i-1];
        const struct reb_particle pi = particles[i];
        const double pme = eta[i]*ei;
        p_j[i].x = pi.x - s_x*ei;
        p_j[i].y = pi.y - s_y*ei;
        p_j[i].z = pi.z - s_z*ei;
        p_j[i].vx = pi.vx - s_vx*ei;
        p_j[i].vy = pi.vy - s_vy*ei;
        p_j[i].vz = pi.vz - s_vz*ei;
        s_x  = s_x  * pme + p_mass[i].m*p_j[i].x ;
        s_y  = s_y  * pme + p_mass[i].m*p_j[i].y ;
        s_z  = s_z  * pme + p_mass[i].m*p_j[i].z ;
        s_vx = s_vx * pme + p_mass[i].m*p_j[i].vx;
        s_vy = s_vy * pme + p_mass[i].m*p_j[i].vy;
        s_vz = s_vz * pme + p_mass[i].m*p_j[i].vz;
    }
    const double Mtotal  = eta[N-1];
    const double Mtotali = 1./Mtotal;
    p_j[0].x = s_x * Mtotali;
    p_j[0].y = s_y * Mtotali;
    p_j[0].z = s_z * Mtotali;
    p_j[0].vx = s_vx * Mtotali;
    p_j[0].vy = s_vy * Mtotali;
    p_j[0].vz = s_vz * Mtotali;
}

void reb_transformations_inertial_to_jacobi_posvelacc(const struct reb_particle* const particles, struct reb_particle* const p_j, const double* const eta, const struct reb_particle* const p_mass, const int N){
    double s_x = eta[0] * particles[0].x;
    double s_y = eta[0] * particles[0].y;
    double s_z = eta[0] * particles[0].z;
    double s_vx = eta[0] * particles[0].vx;
    double s_vy = eta[0] * particles[0].vy;
    double s_vz = eta[0] * particles[0].vz;
    double s_ax = eta[0] * particles[0].ax;
    double s_ay = eta[0] * particles[0].ay;
    double s_az = eta[0] * particles[0].az;
    for (unsigned int i=1;i<N;i++){
        const double ei = 1./eta[i-1];
        const struct reb_particle pi = particles[i];
        const double pme = eta[i]*ei;
        p_j[i].x = pi.x - s_x*ei;
        p_j[i].y = pi.y - s_y*ei;
        p_j[i].z = pi.z - s_z*ei;
        p_j[i].vx = pi.vx - s_vx*ei;
        p_j[i].vy = pi.vy - s_vy*ei;
        p_j[i].vz = pi.vz - s_vz*ei;
        p_j[i].ax = pi.ax - s_ax*ei;
        p_j[i].ay = pi.ay - s_ay*ei;
        p_j[i].az = pi.az - s_az*ei;
        s_x  = s_x  * pme + p_mass[i].m*p_j[i].x ;
        s_y  = s_y  * pme + p_mass[i].m*p_j[i].y ;
        s_z  = s_z  * pme + p_mass[i].m*p_j[i].z ;
        s_vx = s_vx * pme + p_mass[i].m*p_j[i].vx;
        s_vy = s_vy * pme + p_mass[i].m*p_j[i].vy;
        s_vz = s_vz * pme + p_mass[i].m*p_j[i].vz;
        s_ax = s_ax * pme + p_mass[i].m*p_j[i].ax;
        s_ay = s_ay * pme + p_mass[i].m*p_j[i].ay;
        s_az = s_az * pme + p_mass[i].m*p_j[i].az;
    }
    const double Mtotal  = eta[N-1];
    const double Mtotali = 1./Mtotal;
    p_j[0].x = s_x * Mtotali;
    p_j[0].y = s_y * Mtotali;
    p_j[0].z = s_z * Mtotali;
    p_j[0].vx = s_vx * Mtotali;
    p_j[0].vy = s_vy * Mtotali;
    p_j[0].vz = s_vz * Mtotali;
    p_j[0].ax = s_ax * Mtotali;
    p_j[0].ay = s_ay * Mtotali;
    p_j[0].az = s_az * Mtotali;
}

void reb_transformations_inertial_to_jacobi_acc(const struct reb_particle* const particles, struct reb_particle* const p_j, const double* const eta, const struct reb_particle* const p_mass, const int N){
    double s_ax = eta[0] * particles[0].ax;
    double s_ay = eta[0] * particles[0].ay;
    double s_az = eta[0] * particles[0].az;
    for (unsigned int i=1;i<N;i++){
        const double ei = 1./eta[i-1];
        const struct reb_particle pi = particles[i];
        const double pme = eta[i]*ei;
        p_j[i].ax = pi.ax - s_ax*ei;
        p_j[i].ay = pi.ay - s_ay*ei;
        p_j[i].az = pi.az - s_az*ei;
        s_ax = s_ax * pme + p_mass[i].m*p_j[i].ax;
        s_ay = s_ay * pme + p_mass[i].m*p_j[i].ay;
        s_az = s_az * pme + p_mass[i].m*p_j[i].az;
    }
    const double Mtotal  = eta[N-1];
    const double Mtotali = 1./Mtotal;
    p_j[0].ax = s_ax * Mtotali;
    p_j[0].ay = s_ay * Mtotali;
    p_j[0].az = s_az * Mtotali;
}

void reb_transformations_jacobi_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_j, const double* const eta, const struct reb_particle* const p_mass, const int N){
    const double Mtotal  = eta[N-1];
    double s_x  = p_j[0].x  * Mtotal;
    double s_y  = p_j[0].y  * Mtotal;
    double s_z  = p_j[0].z  * Mtotal;
    double s_vx = p_j[0].vx * Mtotal;
    double s_vy = p_j[0].vy * Mtotal;
    double s_vz = p_j[0].vz * Mtotal;
    for (unsigned int i=N-1;i>0;i--){
        const struct reb_particle pji = p_j[i];
        const double ei = 1./eta[i];
        s_x  = (s_x  - p_mass[i].m * pji.x ) * ei;
        s_y  = (s_y  - p_mass[i].m * pji.y ) * ei;
        s_z  = (s_z  - p_mass[i].m * pji.z ) * ei;
        s_vx = (s_vx - p_mass[i].m * pji.vx) * ei;
        s_vy = (s_vy - p_mass[i].m * pji.vy) * ei;
        s_vz = (s_vz - p_mass[i].m * pji.vz) * ei;
        particles[i].x  = pji.x  + s_x ;
        particles[i].y  = pji.y  + s_y ;
        particles[i].z  = pji.z  + s_z ;
        particles[i].vx = pji.vx + s_vx;
        particles[i].vy = pji.vy + s_vy;
        particles[i].vz = pji.vz + s_vz;
        s_x  *= eta[i-1];
        s_y  *= eta[i-1];
        s_z  *= eta[i-1];
        s_vx *= eta[i-1];
        s_vy *= eta[i-1];
        s_vz *= eta[i-1];
    }
    const double mi = 1./eta[0];
    particles[0].x  = s_x  * mi;
    particles[0].y  = s_y  * mi;
    particles[0].z  = s_z  * mi;
    particles[0].vx = s_vx * mi;
    particles[0].vy = s_vy * mi;
    particles[0].vz = s_vz * mi;
}

void reb_transformations_jacobi_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_j, const double* const eta, const struct reb_particle* const p_mass, const int N){
    const double Mtotal  = eta[N-1];
    double s_x  = p_j[0].x  * Mtotal;
    double s_y  = p_j[0].y  * Mtotal;
    double s_z  = p_j[0].z  * Mtotal;
    for (unsigned int i=N-1;i>0;i--){
        const struct reb_particle pji = p_j[i];
        const double ei = 1./eta[i];
        s_x  = (s_x  - p_mass[i].m * pji.x ) * ei;
        s_y  = (s_y  - p_mass[i].m * pji.y ) * ei;
        s_z  = (s_z  - p_mass[i].m * pji.z ) * ei;
        particles[i].x  = pji.x  + s_x ;
        particles[i].y  = pji.y  + s_y ;
        particles[i].z  = pji.z  + s_z ;
        s_x  *= eta[i-1];
        s_y  *= eta[i-1];
        s_z  *= eta[i-1];
    }
    const double mi = 1./eta[0];
    particles[0].x  = s_x  * mi;
    particles[0].y  = s_y  * mi;
    particles[0].z  = s_z  * mi;
}

void reb_transformations_jacobi_to_inertial_acc(struct reb_particle* const particles, const struct reb_particle* const p_j, const double* const eta, const struct reb_particle* const p_mass, const int N){
    const double Mtotal  = eta[N-1];
    double s_ax  = p_j[0].ax  * Mtotal;
    double s_ay  = p_j[0].ay  * Mtotal;
    double s_az  = p_j[0].az  * Mtotal;
    for (unsigned int i=N-1;i>0;i--){
        const struct reb_particle pji = p_j[i];
        const double ei = 1./eta[i];
        s_ax  = (s_ax  - p_mass[i].m * pji.ax ) * ei;
        s_ay  = (s_ay  - p_mass[i].m * pji.ay ) * ei;
        s_az  = (s_az  - p_mass[i].m * pji.az ) * ei;
        particles[i].ax  = pji.ax  + s_ax ;
        particles[i].ay  = pji.ay  + s_ay ;
        particles[i].az  = pji.az  + s_az ;
        s_ax  *= eta[i-1];
        s_ay  *= eta[i-1];
        s_az  *= eta[i-1];
    }
    const double mi = 1./eta[0];
    particles[0].ax  = s_ax  * mi;
    particles[0].ay  = s_ay  * mi;
    particles[0].az  = s_az  * mi;
}

/******************************
 * Democratic heliocentric */

void reb_transformations_inertial_to_democratic_heliocentric_posvel(const struct reb_particle* const particles, struct reb_particle* const p_h, const int N){
    p_h[0].x  = 0.;
    p_h[0].y  = 0.;
    p_h[0].z  = 0.;
    p_h[0].vx = 0.;
    p_h[0].vy = 0.;
    p_h[0].vz = 0.;
    p_h[0].m  = 0.;
    for (unsigned int i=0;i<N;i++){
        p_h[0].x  += particles[i].x *particles[i].m;
        p_h[0].y  += particles[i].y *particles[i].m;
        p_h[0].z  += particles[i].z *particles[i].m;
        p_h[0].vx += particles[i].vx*particles[i].m;
        p_h[0].vy += particles[i].vy*particles[i].m;
        p_h[0].vz += particles[i].vz*particles[i].m;
        p_h[0].m  += particles[i].m;
    }
    p_h[0].x  /= p_h[0].m;
    p_h[0].y  /= p_h[0].m;
    p_h[0].z  /= p_h[0].m;
    p_h[0].vx /= p_h[0].m;
    p_h[0].vy /= p_h[0].m;
    p_h[0].vz /= p_h[0].m;
    
    const double m0 = particles[0].m;
    for (unsigned int i=1;i<N;i++){
        p_h[i].x  = particles[i].x  - particles[0].x ;
        p_h[i].y  = particles[i].y  - particles[0].y ;
        p_h[i].z  = particles[i].z  - particles[0].z ;
        double mf = m0 + particles[i].m;
        p_h[i].vx = mf*(particles[i].vx - p_h[0].vx)/m0;
        p_h[i].vy = mf*(particles[i].vy - p_h[0].vy)/m0;
        p_h[i].vz = mf*(particles[i].vz - p_h[0].vz)/m0;
        p_h[i].m  = particles[i].m;
    }
}

void reb_transformations_democratic_heliocentric_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_h, const int N){
    const double mtot = p_h[0].m;
    particles[0].x  = p_h[0].x;
    particles[0].y  = p_h[0].y;
    particles[0].z  = p_h[0].z;
    for (unsigned int i=1;i<N;i++){
        particles[0].x  -= p_h[i].x*particles[i].m/mtot;
        particles[0].y  -= p_h[i].y*particles[i].m/mtot;
        particles[0].z  -= p_h[i].z*particles[i].m/mtot;
    }
    for (unsigned int i=1;i<N;i++){
        particles[i].x = p_h[i].x+particles[0].x;
        particles[i].y = p_h[i].y+particles[0].y;
        particles[i].z = p_h[i].z+particles[0].z;
    }
}

void reb_transformations_democratic_heliocentric_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_h, const int N){
    reb_transformations_democratic_heliocentric_to_inertial_pos(particles,p_h,N);
    const double mtot = p_h[0].m;
    const double m0 = particles[0].m;
    for (unsigned int i=1;i<N;i++){
        double mf = particles[i].m + m0;
        particles[i].vx = m0*p_h[i].vx/mf+p_h[0].vx;
        particles[i].vy = m0*p_h[i].vy/mf+p_h[0].vy;
        particles[i].vz = m0*p_h[i].vz/mf+p_h[0].vz;
    }
    particles[0].vx = p_h[0].vx*(m0-mtot)/m0;
    particles[0].vy = p_h[0].vy*(m0-mtot)/m0;
    particles[0].vz = p_h[0].vz*(m0-mtot)/m0;
    for (unsigned int i=1;i<N;i++){
        double mf = particles[i].m/(particles[i].m + m0);
        particles[0].vx -= p_h[i].vx*mf;
        particles[0].vy -= p_h[i].vy*mf;
        particles[0].vz -= p_h[i].vz*mf;
    }
}