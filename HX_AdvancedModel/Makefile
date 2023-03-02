################################ OPTIONS #######################################
CC           = g++           #Compilador
FLAG_11      = -std=c++11    #Por si quiero compilar con el c++11
FLAG_OPT     = -O3           #-O3 : Compilacion mas lenta, ejecucion mas rapida
                             #      tambien se pueden usar flags -O1, O2, Os ...
                             #      segun si queremos compil. o ejec. rapida
################################################################################

CODE = Base

all:
	$(CC) $(FLAG_OPT) $(FLAG_11) $(CODE).cpp -o $(CODE).o

clean:
	rm *.o $(CODE)
