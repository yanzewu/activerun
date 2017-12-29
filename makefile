CXXFLAGS = -std=c++11 -lpthread -O3
CXX = g++
VPATH = activerun

objects = main.o BrownianForce.o Compute.o Dumper.o Input.o LangevinIntegrator.o mathutil.o MorseForce.o SwimForce.o SwimForce3d.o

activerun : $(objects)
	$(CXX) $(CXXFLAGS) -o activerun-dev $(objects)

main.o : Input.h Force.h Dumper.h System.h Compute.h Context.h Integrator.h Thermo.h \
	threadpool.h dict.h mathutil.h includes.h neighlist.h pbc.h vec.h dict.h dimension.h \
	resources.h output.h

BrownianForce.o : Force.h threadpool.h dict.h mathutil.h System.h Context.h includes.h Input.h \
	neighlist.h pbc.h vec.h dict.h dimension.h resources.h output.h

Compute.o : Compute.h System.h Context.h includes.h Input.h neighlist.h pbc.h threadpool.h \
	resources.h output.h 

Dumper.o : System.h Context.h includes.h Input.h neighlist.h pbc.h vec.h dict.h dimension.h \
	output.h

Input.o : Input.h includes.h vec.h dict.h dimension.h

LangevinIntegrator.o : Integrator.h System.h Context.h includes.h Input.h neighlist.h pbc.h \
	vec.h dict.h dimension.h resources.h output.h

mathutil.o : mathutil.h

MorseForce.o : Force.h threadpool.h dict.h mathutil.h System.h Context.h includes.h Input.h neighlist.h pbc.h \
	vec.h dimension.h resources.h output.h

SwimForce.o : Force.h threadpool.h dict.h mathutil.h System.h Context.h includes.h Input.h \
	neighlist.h pbc.h vec.h dict.h dimension.h resources.h output.h

SwimForce3d.o : Force.h threadpool.h dict.h mathutil.h System.h Context.h includes.h Input.h \
	neighlist.h pbc.h vec.h dict.h dimension.h resources.h output.h 

.PHONY : clean
clean : 
	rm activerun-dev $(objects)
