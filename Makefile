CC=g++
CFLAGS=`root-config --cflags --libs` -Wall -std=c++11 -Wfatal-errors -Wpedantic 
INC_PATH=-Isrc 

SRC = src

#Debug with OPT=-g
#Reduce compilation time and be able to debug with OPT=-O0
ifdef OPT
 OPTM=$(OPT)
else
 OPTM=-O3
endif

OBJS_SRC = $(SRC)/align_wrt_beam_photon.o $(SRC)/boost_CS.o \
	$(SRC)/BinData.o $(SRC)/setup.o
HEADER = $(SRC)/common.hxx $(SRC)/setup.h $(SRC)/function.h

all: main RealData
.PHONY : all

main: main.C $(OBJS_SRC) 
	$(CC) -o $@ $(OPTM) $@.C $(OBJS_SRC) $(INC_PATH) $(CFLAGS)

RealData: RealData.C $(OBJS_SRC) 
	$(CC) -o $@ $(OPTM) $@.C $(OBJS_SRC) $(INC_PATH) $(CFLAGS)

#SRC files
$(SRC)/align_wrt_beam_photon.o: $(SRC)/align_wrt_beam_photon.cpp $(SRC)/functions.h
	$(CC) -o $@ -c $< $(OPTM) -I $(SRC) $(CFLAGS) 

$(SRC)/boost_CS.o: $(SRC)/boost_CS.cpp $(SRC)/functions.h
	$(CC) -o $@ -c $< $(OPTM) -I $(SRC) $(CFLAGS) 

$(SRC)/BinData.o: $(SRC)/BinData.cpp $(SRC)/functions.h
	$(CC) -o $@ -c $< $(OPTM) -I $(SRC) $(CFLAGS)

$(SRC)/setup.o: $(SRC)/setup.cpp $(SRC)/setup.h
	$(CC) -o $@ -c $< $(OPTM) -I $(SRC) $(CFLAGS)

clean:
	$(RM) main RealData $(OBJS_SRC)
