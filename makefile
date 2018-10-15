output: Source2.o nrerror2.o dmatrix2.o free_dmatrix2.o dvector.o free_dvector.o mgfas5.o mgfas4.o nrfunc4.o nrfunc3.o nrfunc2.o trapzd.o integrat.o Stat3.o Stat2.o
	g++ Source2.o nrerror2.o dmatrix2.o free_dmatrix2.o dvector.o free_dvector.o mgfas5.o mgfas4.o nrfunc4.o nrfunc3.o nrfunc2.o trapzd.o integrat.o Stat3.o Stat2.o -o output
Source2.o: Source2.cpp
	g++ -c Source2.cpp
nrerror2.o: nrerror2.cpp Header1.h
	g++ -c nrerror2.cpp
dmatrix2.o: dmatrix2.cpp Header2.h
	g++ -c dmatrix2.cpp
dvector.o: dvector.cpp Header2.h
	g++ -c dvector.cpp
free_dmatrix2.o: free_dmatrix2.cpp Header2.h
	g++ -c free_dmatrix2.cpp
free_dvector.o: free_dvector.cpp Header2.h
	g++ -c free_dvector.cpp
mgfas5.o: mgfas5.cpp
	g++ -c mgfas5.cpp
mgfas4.o: mgfas4.cpp
	g++ -c mgfas4.cpp
nrfunc4.o: nrfunc4.cpp
	g++ -c nrfunc4.cpp
nrfunc3.o: nrfunc3.cpp
	g++ -c nrfunc3.cpp
nrfunc2.o: nrfunc2.cpp
	g++ -c nrfunc2.cpp
trapzd.o: trapzd.cpp Header5.h
	g++ -c trapzd.cpp
integrat.o: integrat.cpp Header6.h
	g++ -c integrat.cpp
Stat3.o: Stat3.cpp
	g++ -c Stat3.cpp
Stat2.o: Stat2.cpp
	g++ -c Stat2.cpp
clean:
	rm *.o output
