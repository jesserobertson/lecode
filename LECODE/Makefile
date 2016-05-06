# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Created on: Aug 14, 2012
# Copyright (C) 2012 CSIRO
# Description : Makefile for LECODE
# ============================================================================
TOP=$(shell pwd)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

DIRMODS= Classes IOdeclaration PluginUDW Initialiser PluginSeismic PluginRMS PluginOcean Processes Finaliser 

SOURCES = SPMDriver.f90
OBJS=$(SOURCES:.f90=.o)

.PHONY : all dist plugin dust clobber

all: dist

dist: 
	@echo
	@echo "*************************************************"
	@echo "LECODE Author Tristan Salles "
	@echo "*************************************************"
	@echo
	@mkdir -p $(BUILDDIR)	
	@mkdir -p $(OBJDIR)
	@mkdir -p $(MODDIR)
	@mkdir -p $(LIBDIR)
	@mkdir -p bin
	for i in $(DIRMODS) ; do   \
    	  ( cd $$i ; make dist) ;       \
	done
	@echo "*************************************************"
	@echo	
	@echo "Build LECODE binary."
	@echo	
	@echo "*************************************************"
	@$(if $(wildcard SPMDriver.o),rm -f SPMDriver.o,)	
	make $(EXEC)

plugin : 
	cd $(DIRPLUG); make plugin; 
	@echo "*************************************************"
	@echo	
	@echo "LECODE shared library created."
	@echo	
	@echo "*************************************************"

$(EXEC) :	$(OBJS)
	$(F90) $(FFLAGS)  $(FOXFLAGS) $(H5FLAGS) ${KDTREEFLAGS} -I$(TRIANGLELIB) -o $@ $^ \
	$(TRIANGLELIB)/*.o  $(LDFLAGS) -lLECODE  $(H5LDFLAGS) $(H5LIBS) $(LDFOXFLAGS) ${KDTREELDFLAGS} $(KDTREELIBS) 
	#rm *.o
	@echo "*************************************************"
	@echo	
	@echo "LECODE updated in ./bin/."
	@echo	
	@echo "*************************************************"

%.o : %.f90
	$(F90) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) ${KDTREEFLAGS} -c $< -o $@ 
	$(AR) $(LIBDIR)/libLECODE.a $(OBJDIR)/*.o
	
dust :
	for i in $(DIRMODS) ; do   \
    	( cd $$i ; make dust) ;       \
	done
	$(foreach module,$(MODULES),cd $(module); cd - ; )
	rm -fv *~ *.bak *.o *.mod *.original

clobber : dust
	for i in $(DIRMODS) ; do   \
    	( cd $$i ; make clobber) ;   \
	done	
	rm -rfv $(BUILDDIR)
