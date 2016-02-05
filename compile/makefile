#-------------------------------------------------------------------
# Makefile for the TIGGE CMORPH subsetting program
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Paths to the source directories and libraries
#--------------------------------------------------------------------

 ROOT_DIR    = /glade/u/home/tcram/cmorph_v0.x
 SOURCE      = $(ROOT_DIR)/source
 SOURCE_EXEC = $(ROOT_DIR)/bin
 NETCDF_BASE = /glade/apps/opt/netcdf/4.2/intel/default
 
#--------------------------------------------------------------------
# Name the executable
#--------------------------------------------------------------------

 EXECUTABLE = cmorph_v0.x_subset

#--------------------------------------------------------------------
# Compilers and options
#--------------------------------------------------------------------

# evans and mirage nodes configuration
 CPP      = /usr/bin/cpp
 CPP_OPTS = -P -traditional

 COMP     = ifort
 LOAD     = ifort

 FIX      = -fixed -extend_source 72
 FREE     = -free

 OPTS     = -c -real-size 32 -assume byterecl -check bounds

 LIBS        = -L$(NETCDF_BASE)/lib -lnetcdf
 INCL        = -I$(NETCDF_BASE)/include -I$(SOURCE)

 SUF      = f

 COMP_C   = /usr/bin/g++
 OPTS_C   = -c
 LIBS_C   =

#--------------------------------------------------------------------
# Object lists
#--------------------------------------------------------------------

 OBJS = \
      kinds.o \
      cmorph_v0.x_output_nc.o \
      cmorph_v0.x_subset.o
      
#--------------------------------------------------------------------
# compile-link
#--------------------------------------------------------------------

$(EXECUTABLE): $(OBJS)
	$(LOAD) $(OBJS) -o $(SOURCE_EXEC)/$(EXECUTABLE) $(LIBS) 

#--------------------------------------------------------------------
# rules for F90 files
#--------------------------------------------------------------------

kinds.o: $(SOURCE)/kinds.f90
	$(CPP) $(CPP_OPTS) $(SOURCE)/kinds.f90 > $(SOURCE_EXEC)/kinds.$(SUF)
	$(COMP) $(OPTS) $(FREE) $(SOURCE_EXEC)/kinds.$(SUF)

cmorph_v0.x_output_nc.o: $(SOURCE)/cmorph_v0.x_output_nc.f90
	$(CPP) $(CPP_OPTS) $(INCL) $(SOURCE)/cmorph_v0.x_output_nc.f90 > $(SOURCE_EXEC)/cmorph_v0.x_output_nc.$(SUF)
	$(COMP) $(OPTS) $(FREE) $(INCL) $(SOURCE_EXEC)/cmorph_v0.x_output_nc.$(SUF)	
	
cmorph_v0.x_subset.o: $(SOURCE)/cmorph_v0.x_subset.f90
	$(CPP) $(CPP_OPTS) $(INCL) $(SOURCE)/cmorph_v0.x_subset.f90 > $(SOURCE_EXEC)/cmorph_v0.x_subset.$(SUF)
	$(COMP) $(OPTS) $(FREE) $(INCL) $(SOURCE_EXEC)/cmorph_v0.x_subset.$(SUF)	
	
#--------------------------------------------------------------------
# clean
#--------------------------------------------------------------------
clean: 
	clear
	rm -f *.o $(SOURCE_EXEC)/*.f *.mod fort.* $(EXECUTABLE)
