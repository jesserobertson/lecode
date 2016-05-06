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

include $(CONFFILE)

TARGET = SPMIO
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES =  StrataI.f90 ExtForcesI.f90 TINI.f90 TINO.f90 \
	FluxesI.f90 ParamI.f90 FlowWalkerI.f90 CheckPoint.f90 \
	FlowWalkerO.f90 SeaO.f90 StrataOut.f90 OceanIO.f90 \
	CarbOrgsI.f90

OBJS=$(FSOURCES:.f90=.o)
  
.PHONY: dist module dust clobber

dist : module 

module : 	
	@echo
	@echo "*************************************************"
	@echo "Build modules and library for Input/Output"
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

clobber : dust