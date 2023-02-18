.SUFFIXES:
.SUFFIXES: .o .F

include user_build_options


OBJS=	module_sf_noahlsm_glacial_only.o \
	module_sf_noahlsm.o \
	module_model_constants.o \
	kwm_date_utilities.o \
	module_io.o \
	module_ascii_io.o \
	module_sfcdif_wrf.o \
	module_Noahlsm_utility.o\
	module_optimal.o

CMD = driver.exe urban_driver.exe

FFLAGS 	=	$(FREESOURCE) $(F90FLAGS)

all:	$(CMD)

driver.exe:	$(OBJS)
	$(COMPILERF90) -o $(@) $(OBJS) simple_driver.o $(NETCDFLIB)


urban_driver.exe:	$(OBJS)
	$(COMPILERF90) -o $(@) $(OBJS) simple_driver_urban.o $(NETCDFLIB)

.F.o:
	$(RM) $(*).f90
	$(CPP) $(CPPMACROS) $(NETCDFINC) $(*).F > $(*).f90
	$(COMPILERF90) -c $(FFLAGS) $(NETCDFINC) $(*).f90
	#$(RM) $(*).f90
	echo ""

clean:
	$(RM) $(OBJS) simple_driver.o simple_driver_urban.o $(CMD) *.mod *~ *.f90

urban_driver.exe:			simple_driver_urban.o
driver.exe:				simple_driver.o
module_sf_noahlsm.o:			module_model_constants.o
module_sf_noahlsm.o:			module_optimal.o
module_sf_noahlsm_glacial_only.o:	module_model_constants.o
module_sf_noahlsm_glacial_only.o:	module_sf_noahlsm.o
simple_driver.o:			module_io.o
simple_driver.o:			kwm_date_utilities.o
simple_driver.o:			module_Noahlsm_utility.o
module_io.o:				module_ascii_io.o
