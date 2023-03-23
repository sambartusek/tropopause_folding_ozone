### name of the executable that will be produced
PROG       = $(INSTALLDIR)/3d_labelling_and_fold_id.exe

DATE=$(shell date +"%Y%m%d")

# complete list of all f90 source files
SRCS  = mo_f2kcli.f90 netcdf_tools.f90 libipo.f 3d_labelling_and_fold_id.f90

# --------------------------------------------------------------------
OBJS = mo_f2kcli.o netcdf_tools.o libipo.o 3d_labelling_and_fold_id.o

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(F90FLAGS) $(OBJS) $(LIBS) -o $@

list:
	@echo "------------------------------------------------"
	@echo "SRCS = $(SRCS)"
	@echo "------------------------------------------------"
	@echo
	@echo "------------------------------------------------"
	@echo "OBJS = $(OBJS)"
	@echo "------------------------------------------------"

zip:
	zip -r 3d_labelling_and_fold_id_$(DATE).zip *
clean:
	rm -f *.mod
	rm -f *.o


distclean: clean
	rm -f $(PROG)

help:
	@echo ''
	@echo '  posible make targets:'
	@echo '  -------------------------------------------------------'
	@echo '   make                 : build all'
	@echo '   make all             : (= make)'
	@echo ''
	@echo '   make clean           : delete libs/mod/obj'
	@echo '   gmake distclean      : clean up distribution'
	@echo ''
	@echo '   gmake zip             : zip source-code'
	@echo '  -------------------------------------------------------'
	@echo ''


%.o: %.f90
	$(F90) $(F90FLAGS) $(LIBS) $(INCLUDES) -c $<
%.o: %.f
	$(FC) $(FFLAGS) $(LIBS) $(INCLUDES) -c $<

## ------------------------------------------------------------------
