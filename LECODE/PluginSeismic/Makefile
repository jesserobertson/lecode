# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2014 CSIRO
#
# Description : Makefile for Seismic plugin
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

TARGET = SPMSEIS
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES =  Seismic.f90

OBJS=$(FSOURCES:.f90=.o)
  
.PHONY: dist module dust clobber

dist : module 

module :
	@echo
	@echo "*************************************************"
	@echo "Build modules and library for Seismic plugin"
	@echo	
	@echo "*************************************************"
	make $(TARGETLIB)
	
$(TARGETLIB) : $(OBJS)
	#$(AR) $(TARGETLIB) $(OBJS) 
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR); \
		fi; \
	done
	
%.o : %.f90
	$(F90) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) ${KDTREEFLAGS} -c $< -o $@ 
	
dust :
	rm -f *.o *.mod *.original

clobber : dust