# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2012 CSIRO
#
# Description : Makefile for Input/Output
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

FACC = ./FlowAccumulation

include $(CONFFILE)

TARGET = SPMInit
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES =  FillDepressionDEM.f90 ExtForces_init.f90 \
	Strata_init.f90  Porosity_init.f90 TIN_init.f90 Sort_algo.f90  \
	TIN_gen.f90 Rain_init.f90  River_init.f90 \
	Param_init.f90 License.f90 Sim_init.f90 UpdateStrata.f90
	
CPPSOURCES= $(FACC)/Flow_accumulation.cpp

FOBJS=$(FSOURCES:.f90=.o)
CPPOBJS=$(CPPSOURCES:.cpp=.o)
 
OBJS=./*.o $(FACC)/*.o
  
.PHONY: dist module dust clobber

dist : module 

module : 	
	@echo
	@echo "*************************************************"
	@echo "Build modules and library for Initialisation"
	@echo	
	@echo "*************************************************"
	make $(TARGETLIB)
	
$(TARGETLIB) : $(CPPOBJS) $(FOBJS)
	#$(AR) $(TARGETLIB) $(OBJS) 
	cp *.o $(OBJDIR)
	cp $(FACC)/*.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done
				
%.o : %.f90
	$(F90) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) ${KDTREEFLAGS} -c $< -o $@ 

%.o : %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@ 
	
dust :
	rm -f *.o *.mod *.original $(FACC)/*.o

clobber : dust