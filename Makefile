AMREX_LIBRARY_HOME ?= ${AMREX_INSTALL_DIR}

LIBDIR := $(AMREX_LIBRARY_HOME)/lib
INCDIR := $(AMREX_LIBRARY_HOME)/include

COMPILE_CPP_FLAGS ?= $(shell awk '/Cflags:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)
COMPILE_LIB_FLAGS ?= $(shell awk '/Libs:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)

CFLAGS := -I$(INCDIR) $(COMPILE_CPP_FLAGS)
LFLAGS := -L$(LIBDIR) $(COMPILE_LIB_FLAGS)


# Compile using gcc

all:
	nvcc  pp_main.cpp -c  $(CFLAGS) 
	nvcc  particle_utils.cpp -c $(CFLAGS) 
link:
	nvcc pp_main.o particle_utils.o -o main $(LFLAGS) 
