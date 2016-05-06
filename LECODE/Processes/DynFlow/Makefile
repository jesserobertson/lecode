# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2012 CSIRO
#
# Description : Makefile for computing flow dyanmics
# ============================================================================
UP=../..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

TARGET = DynFlow
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES =  ODE_solver.f90 SlopeFlow.f90 DirectionCpt.f90 \
	RecordFW.f90 CleanFW.f90 

OBJS=$(FSOURCES:.f90=.o)
  
.PHONY: dist module dust clobber

dist : module 

module : 	
	@echo
	@echo "*************************************************"
	@echo "Build modules and library for flow dynamics"
	@echo	
	@echo "*************************************************"
	make $(TARGETLIB)
	
$(TARGETLIB) : $(OBJS) 
    #$(AR) $(TARGETLIB) $(OBJS)
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done
				
%.o : %.f90
	$(F90) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) ${KDTREEFLAGS} -c $< -o $@ 
		
dust :
	rm -f *.o *.mod *.original
	rm -f *~ *_genmod* 

clobber : dust