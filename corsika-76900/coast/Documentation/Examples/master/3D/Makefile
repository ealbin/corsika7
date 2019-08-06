CCC = g++       # compiler
DICT = rootcint # root dictionary-generator
AR = ar         # archiver
LD = g++        # linker

ifndef COAST_DIR
  COAST_DIR=../../../../..
endif

ARFLAGS = rcs

LDFLAGS += -fPIC  -ggdb3
LDFLAGS += $(shell root-config --libs)

LDFLAGS += -L$(COAST_DIR)/lib -L$(COAST_DIR)/lib/unknown 
LDFLAGS += -lCorsikaIntern
LDFLAGS += -lCorsikaROOT
LDFLAGS += -lCorsikaFileIO
LDFLAGS += -lCorsikaToROOT

CPPFLAGS += -c -fPIC -ggdb3 -Wall
CPPFLAGS += -I. -I$(COAST_DIR)/include
CPPFLAGS += $(shell root-config --cflags)

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, ${SRCFILES})

OBJSO = libCOAST.so

EXE = 

all: ${EXE} ${OBJSO}

.cc.o:
	${CCC} ${CPPFLAGS} $^

libCOAST.so: Interface.o TPlotter.o 
	${LD} ${LDFLAGS} --shared $^ ${LIBS} -o $@ 

clean: 
	@rm -f *.o *.a *.so *~ core *.root MakeAnim Makefile.dep
	@rm -rfv $(OBJSO) $(EXE)
