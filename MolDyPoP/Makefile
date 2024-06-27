ROOT=/usr/local/
ROOT=$(HOME)/bin/
#LIB=libspp
#LIBS= $(LIB)2d.a $(LIB)3d.a
EXEC=ParticleDynamics
#SRCS=	computations.cpp	dumbdebug.cpp	parameters.cpp	topology.cpp	xyparticle.cpp	xygroup.cpp	sampler.cpp	integrator.cpp	inputoutput.cpp	routines.cpp	main.cpp
 #SRCS=	computations.cpp	dumbdebug.cpp	parameters.cpp	topology.cpp	partition.cpp	group.cpp	neighbor_list.cpp	sampler.cpp	integrator.cpp	inputoutput.cpp	routines.cpp	main.cpp
 SRCS=	computations.cpp	parameters.cpp	topology.cpp	partition.cpp	group.cpp	neighbor_list.cpp	sampler.cpp	integrator.cpp	inputoutput.cpp	routines.cpp	main.cpp


#SRCS=	computations.cpp	topology.cpp	particle.cpp	group.cpp	alternativemain.cpp
#SRCS=	particle.cpp main.cpp	computations.cpp	topology.cpp
OBJS=$(SRCS:.cpp=.o)
#OBJS2D=$(SRCS:.cpp=_2d.o)
#OBJS3D=$(SRCS:.cpp=_3d.o)
HDRS=$(SRCS:.cpp=.h)
COMP= g++
# CFLAGS= -c -Wall -O3 -ffast-math -fopenmp
# LFLAGS= -Wall -O3 -ffast-math -fopenmp
# CFLAGS= -c -Wall -O3 -ffast-math
# LFLAGS= -Wall -O3 -ffast-math
#CFLAGS= -std=c++14 -c -Wall -O0 -g
#LFLAGS= -Wall -O0 -g

# Standard Flags
CFLAGS= -std=c++14 -c -Wall -O2 -g
LFLAGS= -Wall -O2 -g
#For debugging
#CFLAGS= -c -Wall -g
#LFLAGS= -Wall -g

all: $(EXEC)

$(EXEC): $(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(EXEC)

$(OBJS): $(SRCS)
	$(COMP) $(CFLAGS) $(SRCS)
#particle.o: particle.cpp
#	$(COMP) $(CFLAGS) particle.cpp
#
#main.o: main.cpp
#	$(COMP) $(CFLAGS) main.cpp

# g++ -c particle.cpp -o particlemd
#.c.o:
#	$(COMP) -c $(LFLAGS) $(OBJS)

clean:
	rm $(OBJS)
	rm ${EXEC}

allclean:
	rm $(OBJS) $(EXEC) *.xyz *.out *.dat

doxy:
	doxygen Doxyfile
	cd ../doc/latex/
	make
	cd -
