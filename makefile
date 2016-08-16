DLSODAR = opkda1.o opkda2.o opkdmain.o
SLSODAR = opksa1.o opksa2.o opksmain.o
EQUATIONS = eqs_cowell.o eqs_ks.o eqs_edromo.o eqs_gdromo.o
EVENTS = evts_cowell.o evts_ks.o evts_edromo.o evts_gdromo.o sundman.o

OBJECTS = kinds.o constants.o io.o settings.o \
auxiliaries.o third_bodies.o state_init.o transform.o processing.o \
propagate_cowell.o reference_trajectory.o integration.o step_sizes.o \
dp_trajectory.o xeverhart.o $(DLSODAR) $(SLSODAR) $(EQUATIONS) $(EVENTS) NAPLES.o

# Main target - binary
NAPLES.x: NAPLES.f90 $(OBJECTS)
	gfortran -o NAPLES.x -g -fcheck=bounds  $(OBJECTS)

# OBJECT FILES
kinds.o: kinds.f90
	gfortran -c kinds.f90 -g -fcheck=bounds 

constants.o: constants.f90 kinds.o
	gfortran -c constants.f90 -g -fcheck=bounds 

settings.o: settings.f90 kinds.o
	gfortran -c settings.f90 -g -fcheck=bounds 

type_definitions.o: type_definitions.f90 kinds.o
	gfortran -c type_definitions.f90 -g -fcheck=bounds 

auxiliaries.o: auxiliaries.f90 kinds.o
	gfortran -c auxiliaries.f90 -g -fcheck=bounds 

third_bodies.o: third_bodies.f90 kinds.o
	gfortran -c third_bodies.f90 -g -fcheck=bounds 

state_init.o: state_init.f90 kinds.o
	gfortran -c state_init.f90 -g -fcheck=bounds 

transform.o: transform.f90 kinds.o
	gfortran -c transform.f90 -g -fcheck=bounds 

processing.o: processing.f90 kinds.o
	gfortran -c processing.f90 -g -fcheck=bounds

sundman.o: sundman.f90 kinds.o settings.o transform.o
	gfortran -c sundman.f90 -g -fcheck=bounds

# EQUATIONS OF MOTION
eqs_cowell.o: eqs_cowell.f90 kinds.o constants.o third_bodies.o auxiliaries.o
	gfortran -c eqs_cowell.f90 -g -fcheck=bounds 

eqs_ks.o: eqs_ks.f90 kinds.o transform.o third_bodies.o constants.o auxiliaries.o
	gfortran -c eqs_ks.f90 -g -fcheck=bounds 

eqs_edromo.o: eqs_edromo.f90 kinds.o settings.o transform.o third_bodies.o \
constants.o auxiliaries.o
	gfortran -c eqs_edromo.f90 -g -fcheck=bounds 

eqs_gdromo.o: eqs_gdromo.f90 kinds.o settings.o transform.o third_bodies.o \
constants.o auxiliaries.o
	gfortran -c eqs_gdromo.f90 -g -fcheck=bounds 

# EVENTS
evts_cowell.o: evts_cowell.f90 kinds.o constants.o
	gfortran -c evts_cowell.f90 -g -fcheck=bounds 

evts_ks.o: evts_ks.f90 kinds.o auxiliaries.o constants.o
	gfortran -c evts_ks.f90 -g -fcheck=bounds 

evts_edromo.o: evts_edromo.f90 kinds.o auxiliaries.o constants.o settings.o
	gfortran -c evts_edromo.f90 -g -fcheck=bounds 

evts_gdromo.o: evts_gdromo.f90 kinds.o auxiliaries.o constants.o settings.o
	gfortran -c evts_gdromo.f90 -g -fcheck=bounds 

# PROPAGATION, MISC
step_sizes.o: step_sizes.f90 kinds.o constants.o
	gfortran -c step_sizes.f90 -g -fcheck=bounds 

propagate_cowell.o: propagate_cowell.f90 kinds.o constants.o\
eqs_cowell.o evts_cowell.o settings.o auxiliaries.o third_bodies.o $(DLSODAR)
	gfortran -c propagate_cowell.f90 -g -fcheck=bounds 

reference_trajectory.o: reference_trajectory.f90 kinds.o constants.o \
settings.o propagate_cowell.o auxiliaries.o
	gfortran -c reference_trajectory.f90 -g -fcheck=bounds 

dp_trajectory.o: dp_trajectory.f90 kinds.o constants.o \
state_init.o settings.o integration.o auxiliaries.o step_sizes.o $(EQUATIONS) \
$(EVENTS)
	gfortran -c dp_trajectory.f90 -g -fcheck=bounds 

integration.o: integration.f90 kinds.o auxiliaries.o settings.o constants.o \
transform.o $(SLSODAR) xeverhart.o
	gfortran -c integration.f90 -g -fcheck=bounds 

io.o: io.f90 kinds.o settings.o
	gfortran -c io.f90 -g -fcheck=bounds 

NAPLES.o: NAPLES.f90 kinds.o constants.o io.o settings.o auxiliaries.o io.o \
reference_trajectory.o dp_trajectory.o processing.o third_bodies.o integration.o
	gfortran -c NAPLES.f90 -g -fcheck=bounds 


# LSODAR - DOUBLE PRECISION:
opksa1.o: LSODAR/opksa1.f
	gfortran -c LSODAR/opksa1.f -g -fdefault-real-8 -std=legacy

opksa2.o: LSODAR/opksa2.f
	gfortran -c LSODAR/opksa2.f -g -fdefault-real-8 -std=legacy

opksmain.o: LSODAR/opksmain.f opksa1.o opksa2.o
	gfortran -c LSODAR/opksmain.f -g -fdefault-real-8 -std=legacy
	

# LSODAR - QUADRUPLE PRECISION (compatible with double):
opkda1.o: LSODAR/opkda1.f
	gfortran -c LSODAR/opkda1.f -g -fdefault-real-8 -std=legacy

opkda2.o: LSODAR/opkda2.f
	gfortran -c LSODAR/opkda2.f -g -fdefault-real-8 -std=legacy

opkdmain.o: LSODAR/opkdmain.f opkda1.o opkda2.o
	gfortran -c LSODAR/opkdmain.f -g -fdefault-real-8 -std=legacy

# RADAU
#everhart.o: everhart.f90 kinds.o
#	gfortran -c everhart.f90 -g -fcheck=bounds

xeverhart.o: xeverhart.f90 kinds.o sundman.o
	gfortran -c xeverhart.f90 -g -fcheck=bounds 

.PHONY: clean
clean:
	-rm *.x *.o *.mod
