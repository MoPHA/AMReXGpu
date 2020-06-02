#include "AMReX_Array.H"
#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_DistributionMapping.H"
#include "AMReX_Geometry.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParGDB.H"
#include "AMReX_ParallelDescriptor.H"
#include <AMReX_PlotFileUtil.H>
#include "AMReX_ParallelReduce.H"
#include "AMReX_REAL.H"
#include <iostream>
#include <math.h>
#include <utility>
#include <vector>
#include "AMReX_VisMF.H"
#include "const_defs.hpp"

double uniform_density(const amrex::Geometry geom,int i ,int j ,int k){
    return 1;
}


 void add_single_particle( CParticleTile&particlet ,amrex::RealArray pos , amrex::RealArray vel){ 
    CParticle p;
    p.id()   = CParticle::NextID();
    p.cpu()  = amrex::ParallelDescriptor::MyProc();
    p.pos(X) = pos[VX];
    p.pos(Y) = pos[VY];
    p.pos(Z) = pos[VZ];
    p.rdata(VX)=vel[VX];
    p.rdata(VY)=vel[VY];
    p.rdata(VZ)=vel[VZ];
    particlet.push_back(p);
}

void add_single_particle(CParticleContainer&P,amrex::RealArray pos , amrex::RealArray vel){
for(amrex::MFIter mfi= P.MakeMFIter(0) ;mfi.isValid();++mfi){
    
    // Each grid,tile has a their own local particle container
    auto& particles = P.GetParticles(0)[std::make_pair(mfi.index(),mfi.tileIndex())];
        if(mfi.index()==0){
            add_single_particle(particles,pos,vel);
        }
    }
    P.Redistribute();

}


void add_particle_density(const amrex::Geometry geom , CParticleContainer&P, double (*dist_func)(const amrex::Geometry,int,int,int),int ppc_max , double v){

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0,1);
    std::normal_distribution<double> vel(0,v); 

// For simplicity, particles are initialized with Maxwell–Boltzmann distribution


for(amrex::MFIter mfi= P.MakeMFIter(0) ;mfi.isValid();++mfi){
    
    // Each grid,tile has a their own local particle container
    auto& particles = P.GetParticles(0)[std::make_pair(mfi.index(),
                                        mfi.LocalTileIndex())];
    auto box=mfi.validbox();
   const auto lo = amrex::lbound(box);
   const auto hi = amrex::ubound(box);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
            
           double x = geom.ProbLo(X) + i*geom.CellSize(X);
           double y = geom.ProbLo(Y) + j*geom.CellSize(Y);
           double z = geom.ProbLo(Z) + k*geom.CellSize(Z);
           int num_particles = dist_func(geom,i,j,k)*ppc_max;  
            for(int p =0; p < num_particles ; p++){ 
                add_single_particle(particles,{x+dist(mt)*geom.CellSize(X),y+dist(mt)*geom.CellSize(Y),z+dist(mt)*geom.CellSize(Z)},{vel(mt),vel(mt),vel(mt)});
            }
       }
     }
   }
 }
P.Redistribute();

}
