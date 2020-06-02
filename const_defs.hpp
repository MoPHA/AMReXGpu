// Some compiletime constants for convinience

#ifndef particledefs
#define particledefs
#include "AMReX_Particles.H"
#include<AMReX.H>


// CParticle is a general charged particle with the following info + position.

// mass 
// charge
// vx
// vy
// vz


// Se amrex documentation for more information
// https://amrex-codes.github.io/amrex/docs_html/Particle.html
#define C_NUM_REALS 3 
#define C_NUM_INTS 0 
#define C_NUM_SOA_REALS 0 
#define C_NUM_SOA_INTS 0 


#define X 0
#define Y 1
#define Z 2

// 
#define VX 0
#define VY 1
#define VZ 2


typedef amrex::ParIter< C_NUM_REALS  ,C_NUM_INTS,C_NUM_SOA_REALS,C_NUM_SOA_INTS> CParIter;
typedef amrex::ParticleContainer<C_NUM_REALS,C_NUM_INTS,C_NUM_SOA_INTS,C_NUM_SOA_REALS> CParticleContainer;
typedef amrex::Particle<C_NUM_REALS,C_NUM_INTS> CParticle;
typedef amrex::ArrayOfStructs<C_NUM_REALS,C_NUM_INTS> CParticles;
typedef amrex::ParticleTile<C_NUM_REALS  ,C_NUM_INTS,C_NUM_SOA_REALS,C_NUM_SOA_INTS> CParticleTile;

#endif
