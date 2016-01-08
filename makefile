all:
	g++ cfd2d.cpp riemann.h  
clean:
	rm -f -r a.out *.gch *~ *output.csv
