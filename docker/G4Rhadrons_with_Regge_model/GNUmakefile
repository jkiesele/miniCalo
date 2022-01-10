# $Id: GNUmakefile,v 1.6 2006/07/11 08:26:05 mackepr Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

PROJECTID=G4Code
BACKUPLOCATIONS=heppc18:scrlocal/paperbackup lxplus:paperbackup/

name := exampleRhad
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = /usr/local/share/Geant4-10.6.2/geant4make
endif

##uncomment this and see exampleRhad.cc
##to enable hadronic interactions
#include hadronic_lists.gmk


.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

#LDFLAGS += -L$(XERCESROOT)/lib

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

purge:
	rm -f include/*~ src/*~ *~

backup:
	@backupdate=`date +'%y%m%d_%Hh%Mm'`;\
	thispwd=`pwd`;\
	thisdir=`basename $$thispwd`;\
	rm -rf tmp bin;\
	cd ../; tar  cvzf ${PROJECTID}_backup_$$backupdate.tgz --exclude=logfiles/* $$thisdir;\
	for i in ${BACKUPLOCATIONS}; do\
	echo Placing copy at $$i;\
	scp ${PROJECTID}_backup_$$backupdate.tgz $$i/;\
	done

