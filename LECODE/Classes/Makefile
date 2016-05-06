# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2012 CSIRO
#
# Description : Makefile for Strata Classes
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

FTRIA = ./Triangle

include $(CONFFILE)

TARGET = Classes
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

SOURCES= ClassGeneral.f90  ClassFile.f90 ClassStratal.f90 \
					ClassForces.f90 ClassError.f90 ClassTIN.f90 \
					ClassFluxes.f90  
							
CSOURCES= $(FTRIA)/TriangleCall.c

FOBJS=$(SOURCES:.f90=.o)
COBJS=$(CSOURCES:.c=.o)
OBJS=./*.o $(FTRIA)/*.o 

.PHONY: dist dust clobber 

dist : module

module :	
	@echo "*************************************************"
	@echo	
	@echo "Build modules and library for Classes."
	@echo	
	@echo "*************************************************"
	make $(TARGETLIB)
	
$(TARGETLIB): $(COBJS) $(FOBJS)
	#$(AR) $(TARGETLIB) $(OBJS) 
	cp *.o $(OBJDIR)
	cp $(FTRIA)/*.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done

%.o : %.f90
	$(F90) $(FFLAGS) ${KDTREEFLAGS} -c $< -o $@ 
			
%.o : %.c
	$(CC) $(CFLAGS) -I$(TRIANGLELIB) -c $< -o $@   
	
dust:
	rm -fv *.bak *~ *.o *.mod $(FTRIA)/*.o *.original 
	
clobber: dust
