# Examples of very basic AMReX GPU functionality 

**Note** features are demonstrated in the GPU implementation of the [strugepic library](https://github.com/MoPHA/strugepic)

Below are very small example code snippets for the most common operation, see strugepic for more context. 


## Code snippets

Following typdefs are used 
```
typedef amrex::ParIter< N_REALS  ,N_INTS,N_SOA_REALS,N_SOA_INTS> CParIter;                 
typedef amrex::ParticleContainer<N_REALS,N_INTS,N_SOA_INTS,N_SOA_REALS> CParticleContainer;
typedef amrex::Particle<N_REALS,N_INTS> CParticle;                                         
typedef amrex::ArrayOfStructs<N_REALS,N_INTS> CParticles;                                  

```
Fields and particles 
```
CParticleContainer P(geom,dm,ba);
amrex::MultiFab E(ba,dm,Ncomp,Nghost);
amrex::MultiFab B(ba,dm,Ncomp,Nghost);
```


### Particle to Particle
```c++
    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        CParticle *  AMREX_RESTRICT particles= &(pti.GetArrayOfStructs()[0]);
        const int np = pti.numParticles();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (long i)
            {
               particles[i].pos(X)= particles[i].rdata(VX)*dt; 
               ...
              
            });

        }
        
```
### Grid to Grid
```c++
 for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        auto const& E_loc = E.array(mfi);
        auto const& B_loc = B.const_array(mfi);
         
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE  (int i,int j,int k )
        {
         E_loc(i,j,k,X) = B_loc(i+1,j,k,Y);  
         ...
        });
```

### Grid to Particle
```c++
  for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto box=E.box(pti.index());
        CParticle *  AMREX_RESTRICT particles= &(pti.GetArrayOfStructs()[0]);
        const int np = pti.numParticles();
        amrex::Array4<amrex::Real> const& E_loc = E[pti].array();
        
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (long i)
        {
        int x =floor(particles[i].pos(X));
        int y =floor(particles[i].pos(Y));
        int z =floor(particles[i].pos(Z));
        particles[i].rdata(VX)+=E_loc(x,y,z,X)*dt*particles[i].rdata(Q);
        ...
        });
        
```
### Particle to grid
```c++
for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto box=E.box(pti.index());
        CParticle *  AMREX_RESTRICT particles= &(pti.GetArrayOfStructs()[0]);
        const int np = pti.numParticles();
        amrex::Array4<amrex::Real> const& E_loc = E[pti].array();
        
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (long i)
        {
        int x =floor(particles[i].pos(X));
        int y =floor(particles[i].pos(Y));
        int z =floor(particles[i].pos(Z));
        amrex::Gpu::Atomic::Add(&E_loc(x,y,z,X),particles[i].rdata(VX))
        ...
        });
```
### Reduction
```c++
std::pair<amrex::Real,amrex::Real> get_total_energy(const amrex::Geometry geom,CParticleContainer&P,
                                                    amrex::MultiFab &E, amrex::MultiFab &B ){
    
    amrex::Gpu::DeviceVector<amrex::Real> EkinVect(20, 0.); 
    amrex::Real* EkinPtr = EkinVect.dataPtr();

    // Kinetic part 
    amrex::Real E_field=0;
    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        const int np  = pti.numParticles();
        auto particles = &(pti.GetArrayOfStructs()[0]);
        amrex::ParallelFor( np,[=] AMREX_GPU_DEVICE (long i)
        {
            amrex::Real tmp =particles[i].rdata(M)*0.5*
                    ( particles[i].rdata(VX)*particles[i].rdata(VX)+
                      particles[i].rdata(VY)*particles[i].rdata(VY)+
                      particles[i].rdata(VZ)*particles[i].rdata(VZ) 
                     );
            amrex::Gpu::Atomic::Add(EkinPtr + (i%20), tmp);
        });
    }

    // reduce the ekin vector to one final scalar
    amrex::Real E_kin = amrex::Reduce::Sum(20, EkinPtr);
    amrex::Gpu::synchronize();

    auto E_L2_norm=E.norm2({X,Y,Z});
    auto B_L2_norm=B.norm2({X,Y,Z});
    E_field+=E_L2_norm[X]*E_L2_norm[X]+E_L2_norm[Y]*E_L2_norm[Y]+E_L2_norm[Z]*E_L2_norm[Z];
    E_field+=B_L2_norm[X]*B_L2_norm[X]+B_L2_norm[Y]*B_L2_norm[Y]+B_L2_norm[Z]*B_L2_norm[Z];
    E_field*=0.5;

    amrex::ParallelAllReduce::Sum(E_kin, amrex::ParallelDescriptor::Communicator());
    return std::make_pair(E_field*geom.CellSize(X)*geom.CellSize(Y)*geom.CellSize(Z), E_kin);
}
```
