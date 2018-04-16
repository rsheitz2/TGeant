CC=g++
#CFLAGS=`root-config --cflags --libs` -Wall
CFLAGS=`root-config --cflags --libs` -std=c++11 -Wfatal-errors
INC_PATH=-ISrc

SRC = Src

#Debug with OPT=-g
#Reduce compilation time and be able to debug with OPT=-O0
ifdef OPT
 OPTM=$(OPT)
else
 OPTM=-O3
endif

OBJS_SRC = $(SRC)/functions.o
HEADER = $(SRC)/common.hxx $(SRC)/function.h
EXE = avgBinBounds genAvgBinBounds

all: $(EXE)
.PHONY : all


avgBinBounds: avgBinBounds.C $(OBJS_SRC) 
	$(CC) -o $@ $(OPTM) $@.C $(OBJS_SRC) $(INC_PATH) $(CFLAGS)

genAvgBinBounds: genAvgBinBounds.C $(OBJS_SRC) 
	$(CC) -o $@ $(OPTM) $@.C $(OBJS_SRC) $(INC_PATH) $(CFLAGS)

#SRC files
$(SRC)/functions.o: $(SRC)/functions.cpp
	$(CC) -o $@ -c $< $(OPTM) -I $(SRC) $(CFLAGS) 

clean:
	$(RM) $(EXE) $(OBJS_SRC)
