#include "const_defs.hpp"

void push_particle_position(amrex::Geometry geom,CParticleContainer &P,double dt){ 
    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        CParticle *  AMREX_RESTRICT particles= &(pti.GetArrayOfStructs()[0]);
        const int np = pti.numParticles();

// Normal arrays can not be used on the GPU!
        amrex::GpuArray<float,3> lb;
        amrex::GpuArray<float,3> ub;
        ub[X] =geom.ProbHi(X);
        ub[Y] =geom.ProbHi(Y);
        ub[Z] =geom.ProbHi(Z);
        lb[X] =geom.ProbLo(X);
        lb[Y] =geom.ProbLo(Y);
        lb[Z] =geom.ProbLo(Z);

  


        AMREX_PARALLEL_FOR_1D(np,i,{                    
                    particles[i].pos(X)+= particles[i].rdata(X)*dt; 
                    particles[i].pos(Y)+= particles[i].rdata(Y)*dt; 
                    particles[i].pos(Z)+= particles[i].rdata(Z)*dt; 
                        
                    for(int comp: {X,Y,Z}){
                        auto pos =particles[i].pos(comp);
                        if(particles[i].pos(comp) > ub[comp] ){
                            particles[i].pos(comp) = lb[comp]+(pos-ub[comp]);
                        }
                        else if(particles[i].pos(comp) < lb[comp]){
                            particles[i].pos(comp) = ub[comp]-(lb[comp]-pos);
                        }
                    }

                    }
                );
    }
    P.Redistribute();

}
