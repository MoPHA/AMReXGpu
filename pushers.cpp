#include "const_defs.hpp"

void push_particle_position(amrex::Geometry geom,CParticleContainer &P,double dt){
    
    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        CParticle *  AMREX_RESTRICT particles= &(pti.GetArrayOfStructs()[0]);
        const int np = pti.numParticles();
        auto lb =geom.ProbLo();
        auto ub =geom.ProbHi();
        AMREX_PARALLEL_FOR_1D(np,i,{
                for(int comp=0;comp<3;comp++){
                    particles[i].pos(comp)+= particles[i].rdata(comp)*dt;
                    double pos=particles[i].pos(comp);
                if(pos > ub[comp]){
                    particles[i].pos(comp) = lb[comp]+(pos-ub[comp]); 
                }
                else if(pos < lb[comp]){
                    particles[i].pos(comp) = ub[comp]-(lb[comp]-pos); 
                } 
                    }
                }
                );
    }
    P.Redistribute();

}
