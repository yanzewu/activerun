CXXFLAGS = -std=c++11 -lpthread -O3
CXX = g++
VPATH = activerun

objects = main.o BrownianForce.o Dumper.o Input.o LangevinIntegrator.o mathutil.o MorseForce.o SwimForce.o SwimForce3d.o

activerun : $(objects)
	$(CXX) $(CXXFLAGS) -o activerun-dev $(objects)

main.o : activerun.h Input.h Force.h Dumper.h System.h Context.h Integrator.h thermostat.h \
	threadpool.h dict.h mathutil.h includes.h neighlist.h pbc.h vec.h dict.h dimension.h

BrownianForce.o : Force.h threadpool.h dict.h mathutil.h System.h Context.h includes.h Input.h \
	neighlist.h pbc.h vec.h dict.h dimension.h

Dumper.o : System.h Context.h includes.h Input.h neighlist.h pbc.h vec.h dict.h dimension.h

Input.o : Input.h includes.h vec.h dict.h dimension.h

LangevinIntegrator.o : Integrator.h System.h Context.h includes.h Input.h neighlist.h pbc.h \
	vec.h dict.h dimension.h

mathutil.o : mathutil.h

MorseForce.o : Force.h threadpool.h dict.h mathutil.h System.h Context.h includes.h Input.h neighlist.h pbc.h \
	vec.h dimension.h

SwimForce.o : Force.h threadpool.h dict.h mathutil.h System.h Context.h includes.h Input.h \
	neighlist.h pbc.h vec.h dict.h dimension.h 

SwimForce3d.o : Force.h threadpool.h dict.h mathutil.h System.h Context.h includes.h Input.h \
	neighlist.h pbc.h vec.h dict.h dimension.h 

.PHONY : clean
clean : 
	rm activerun-dev $(objects)