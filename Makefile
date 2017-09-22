CC=g++
CFLAGS=`root-config --cflags --libs` -Wall
#CFLAGS=`root-config --cflags --libs`
INC_PATH=-Isrc 

SRC = src
KIN = Kinematics

#Debug with OPT=-g
#Reduce compilation time and be able to debug with OPT=-O0
ifdef OPT
 OPTM=$(OPT)
else
 OPTM=-O3
endif

OBJS_SRC = $(SRC)/align_wrt_beam_photon.o $(SRC)/boost_CS.o \
	$(SRC)/BinData.o

main: main.C $(OBJS_SRC) $(SRC)/functions.h
	$(CC) -o $@ $(OPTM) $@.C $(OBJS_SRC) $(INC_PATH) $(CFLAGS)

#SRC files
$(SRC)/align_wrt_beam_photon.o: $(SRC)/align_wrt_beam_photon.cpp $(SRC)/functions.h
	$(CC) -o $@ -c $< $(OPTM) -I $(SRC) $(CFLAGS) 

$(SRC)/boost_CS.o: $(SRC)/boost_CS.cpp $(SRC)/functions.h
	$(CC) -o $@ -c $< $(OPTM) -I $(SRC) $(CFLAGS) 

$(SRC)/BinData.o: $(SRC)/BinData.cpp $(SRC)/functions.h
	$(CC) -o $@ -c $< $(OPTM) -I $(SRC) $(CFLAGS) 


clean:
	$(RM) main $(OBJS) $(OBJS_SRC)
