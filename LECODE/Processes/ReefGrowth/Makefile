# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2014 CSIRO
#
# Description : Makefile for computing reefs evolution
# ============================================================================
UP=../..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

TARGET = ReefDyn
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES =  FunctionCarbs.f90 BuildCarbs.f90 

OBJS=$(FSOURCES:.f90=.o)
  
.PHONY: dist module dust clobber

dist : module 

module : 	
	@echo
	@echo "*************************************************"
	@echo "Build modules and library for carbonate/organic evolution"
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