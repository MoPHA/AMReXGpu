#ifndef particleutils
#define particleutils
#include<AMReX.H>
#include"const_defs.hpp"

double uniform_density(const amrex::Geometry geom,int i ,int j ,int k);
void add_particle_density(const amrex::Geometry geom , CParticleContainer&P, double (*dist_func)(const amrex::Geometry,int,int,int),int ppc_max ,double v);

#endif
