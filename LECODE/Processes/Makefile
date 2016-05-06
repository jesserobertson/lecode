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

FDYN = DynFlow
FSTRA = LayerStatus
RSTRA = ReefGrowth

include $(CONFFILE)

TARGET = SPMProcess
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

DIRMOD= $(FDYN) $(FSTRA) $(RSTRA)

SOURCES = FunctionIdle.f90

OBJS=$(SOURCES:.f90=.o)

.PHONY: dist module dust clobber

dist : module 

module : 	
	@echo
	@echo "*************************************************"
	@echo "Build modules and library for Processes"
	@echo	
	@echo "*************************************************"
	for i in $(DIRMOD) ; do   \
    	  ( cd $$i ; make dist) ;       \
	done
	make $(TARGETLIB)
	
$(TARGETLIB) : $(OBJS)
	#$(AR) $(TARGETLIB) $(OBJS) $(FDYN)/*.o $(FSTRA)/*.o
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done
				
%.o : %.f90
	$(F90) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) ${KDTREEFLAGS} -c $< -o $@ 
	
dust :
	for i in $(DIRMOD) ; do   \
    	( cd $$i ; make dust) ;       \
	done
	rm -f *.o *.mod *.original

clobber : dust
	for i in $(DIRMOD) ; do   \
    	( cd $$i ; make clobber) ;       \
	done