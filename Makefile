#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2011-2012 Clemson University.
# 
#  This file was originally written by Tianhong Yu. 
# 
#  This is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#  USA
# 
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file Makefile
#//! \brief A makefile to generate carbon dust examples.
#//!
#///////////////////////////////////////////////////////////////////////////////

SVNURL = http://svn.code.sf.net/p/nucnet-tools/code/trunk

NUCNET_TARGET = ../nucnet-tools-code
VENDORDIR = $(NUCNET_TARGET)/vendor
OBJDIR = $(NUCNET_TARGET)/obj

NNT_DIR = $(NUCNET_TARGET)/nnt
USER_DIR = $(NUCNET_TARGET)/user
MY_USER_DIR = $(NUCNET_TARGET)/my_user
BUILD_DIR = $(NUCNET_TARGET)/build

USE_SPARSE_SOLVER = no

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

#===============================================================================
# Svn.
#===============================================================================

ifndef SVN_BIN
  SVN_BIN = svn
endif

SVN_CHECKOUT := $(shell if [ ! -d $(NUCNET_TARGET) ]; then $(SVN_BIN) co $(SVNURL) $(NUCNET_TARGET); else $(SVN_BIN) update $(NUCNET_TARGET); fi )

#===============================================================================
# Includes.
#===============================================================================

include $(BUILD_DIR)/Makefile
include $(BUILD_DIR)/Makefile.sparse
include $(USER_DIR)/Makefile.inc

VPATH = $(BUILD_DIR):$(NNT_DIR):$(USER_DIR):$(MY_USER_DIR)

#===============================================================================
# Set the hydro.  CARBON_HYDRO_CODE is an environment variable.  In a bash 
# shell, set this by typing at the command line, for example,
# 'export CARBON_HYDRO_CODE=traj'.  This approach replaces the deprecated
# approach in which you set, for example, CARBON_HYDRO_CODE=traj within this 
# Makefile.
#===============================================================================

ifdef CARBON_HYDRO_CODE
  ifeq ($(CARBON_HYDRO_CODE), default)
        CFLAGS += -DHYDRO_DEFAULT
  else ifeq ($(CARBON_HYDRO_CODE), traj)
        CFLAGS += -DHYDRO_TRAJ
  endif
else
        CFLAGS += -DHYDRO_DEFAULT
endif

#===============================================================================
# Set the user hydro. NNT_HYDRO_CODE is an environment variable.  In a bash 
# shell, set this by typing at the command line, for example,
# 'export NNT_HYDRO_CODE=traj'.  This approach replaces the deprecated
# approach in which you set, for example, NNT_HYDRO_CODE=traj within this 
# Makefile.
#===============================================================================

#ifndef NNT_HYDRO_CODE
  #NNT_HYDRO_CODE=traj
#endif

#===============================================================================
# Carbon object. 
#===============================================================================

CARBON_OBJ = $(OBJDIR)/carbon_hydro.o              \
             $(OBJDIR)/carbon_rate_functions.o     \
             $(OBJDIR)/carbon_evolve.o             \
             $(OBJDIR)/carbon_flow_utilities.o     \
             $(OBJDIR)/my_bin_utilities.o          \
             $(OBJDIR)/carbon_molecule_utilities.o \
             $(OBJDIR)/carbon_reaction_utilities.o 

$(CARBON_OBJ): $(OBJDIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

#===============================================================================
# Objects and dependencies.
#===============================================================================

CARBON_OBJS = $(WN_OBJ)		\
               $(NNT_OBJ)	\
               $(SOLVE_OBJ)	\
               $(CARBON_OBJ)	\
               $(USER_OBJ)

CARBON_DEP = clean_carbon_hydro

ifeq ($(USE_SPARSE_SOLVER), yes)
	CFLAGS += -DSPARSE_SOLVER
	CARBON_OBJS += $(SP_OBJ)
	CARBON_DEP += sparse
	FLIBS= -L$(SPARSKITDIR) -lskit -lstdc++
	MC = $(FF)
else
	MC = $(CC)
endif

CARBON_DEP += $(CARBON_OBJS)

#===============================================================================
# Executables.
#===============================================================================

CARBON_EXEC = run_carbon               \
              run_generic              \
              compute_carbon_flows     \
              create_carbon_network    \
              print_carbon_flows

$(CARBON_EXEC): $(CARBON_DEP)
	$(CC) -c -o $(OBJDIR)/$@.o $@.cpp
	$(MC) $(CARBON_OBJS) $(OBJDIR)/$@.o $(CLIBS) $(FLIBS) -o $(BINDIR)/$@

.PHONY all_carbon : $(CARBON_EXEC)

#===============================================================================
# Make data.
#===============================================================================

carbon_data: 
	create_carbon_network

#===============================================================================
# Clean up.
#===============================================================================

.PHONY: clean_carbon cleanall_carbon clean_carbon_hydro

clean_carbon: 
	rm -f $(CARBON_OBJS)

cleanall_carbon: clean_carbon
	rm -f $(BINDIR)/$(CARBON_EXEC) $(BINDIR)/$(CARBON_EXEC).exe

clean_carbon_hydro:
	rm -f $(OBJDIR)/carbon_hydro.o

#===============================================================================
# Define.
#===============================================================================

CARBON_DEF = yes
