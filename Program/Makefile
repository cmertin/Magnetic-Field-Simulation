CC = g++
CFLAGS = -g -I. -c
LFLAGS = -g

Simulation:   		main.o SolenoidBField.o
	   		$(CC) $(LFLAGS) -o Simulation main.o SolenoidBField.o
main.o:    		main.cpp SolenoidBField.o
	    		$(CC) $(CFLAGS) main.cpp
SolenoidBField.o: 	SolenoidBField.h SolenoidBField.cpp
	    		$(CC) $(CFLAGS) SolenoidBField.h
	    		$(CC) $(CFLAGS) SolenoidBField.cpp
			$(CC) $(CFLAGS) Variables.h

clean:
	rm *.o Simulation *~ *.gch \#*
