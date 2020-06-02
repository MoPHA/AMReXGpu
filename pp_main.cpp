// MWE for particle to particle propagation
// on GPU using amrex

// Amrex
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>
#include <AMReX_CoordSys.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_Particle.H>
#include <AMReX_Utility.H>
#include <AMReX_ParticleUtil.H>
#include <AMReX_NeighborParticles.H>
#include "AMReX_Array4.H"


#include "const_defs.hpp"
#include "particle_utils.hpp"

void main_main();

int main(int argc, char* argv[])
{

    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();

}
void main_main()
{
    amrex::ParmParse pp; 
    // Number of cells in each coordinate direction (x,y,z)
    std::array<int,3> n_cell;
    // Max size for subgrids
    std::array<int,3> max_grid_size;
    
    // Particle per cell
    int ppc;
    // Particle mass and charge
    // Each Particle container only contains indentical particles
    // This test case only has one type of particles
    amrex::Real m;
    amrex::Real q;
    amrex::Real v;

    pp.get("n_cell",n_cell);
    pp.get("max_grid_size",max_grid_size);
    pp.get("ppc",ppc);
    pp.get(" m",m);
    pp.get(" q",q);
    pp.get(" v",v);

    // Periodicity 
    amrex::Vector<int> is_periodic({1,1,1});
    // Cell centered indexing 
    amrex::IndexType typ({AMREX_D_DECL(0,0,0)});


    amrex::IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    amrex::IntVect dom_hi(AMREX_D_DECL(n_cell[X]-1, n_cell[Y]-1, n_cell[Z]-1));
    amrex::Box domain(dom_lo, dom_hi,typ);
    amrex::BoxArray ba(domain);
    ba.maxSize({max_grid_size[X],max_grid_size[Y],max_grid_size[Z]}); 
     amrex::RealBox real_box({AMREX_D_DECL(0,0,0)},
                     {AMREX_D_DECL((double)n_cell[X] , (double)n_cell[Y],(double)n_cell[Z])});

    // This defines a Geometry object
    amrex::Geometry geom(domain,&real_box,amrex::CoordSys::cartesian,is_periodic.data());
    // How Boxes are distrubuted among MPI processes
    amrex::DistributionMapping dm(ba);
    CParticleContainer P(geom,dm,ba);
    add_particle_density(geom,P,uniform_density,ppc,v);    

    amrex::Print() << P.TotalNumberOfParticles();

}
