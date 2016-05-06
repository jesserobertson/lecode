# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2012 CSIRO
#
# Description : Makefile for sedimentary layer evolution
# ============================================================================
UP=../..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

TARGET = DynFlow
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES = DiffuseFct.f90 DiffuseSed.f90 \
	FunctionsFW.f90 LayersFW.f90 EDTransfer.f90 AdvanceFW.f90 

OBJS=$(FSOURCES:.f90=.o)
  
.PHONY: dist module dust clobber

dist : module 

module : 	
	@echo
	@echo "*************************************************"
	@echo "Build modules and library for sedimentary layer evolution"
	@echo	
	@echo "*************************************************"
	make $(TARGETLIB)
	
$(TARGETLIB) : $(OBJS) 
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