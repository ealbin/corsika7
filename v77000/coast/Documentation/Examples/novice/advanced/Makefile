CCC = g++       # compiler
DICT = rootcint # root dictionary-generator
AR = ar         # archiver
LD = g++        # linker

ifndef COAST_DIR
  COAST_DIR=../../../../..
endif

ARFLAGS = rcs

LDFLAGS += -fPIC -ggdb3 -Wl,--no-as-needed

LDFLAGS += -L$(COAST_DIR)/lib -L$(COAST_DIR)/lib/unknown 
LDFLAGS += -lCorsikaToROOT
LDFLAGS += -lCorsikaROOT
LDFLAGS += -lCorsikaFileIO
LDFLAGS += -lCorsikaIntern
LDFLAGS += $(shell root-config --libs)

LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaIntern
LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaROOT
LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaFileIO
LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaToROOT

CPPFLAGS += -c -fPIC -ggdb3
CPPFLAGS += -I. -I$(COAST_DIR)/include
CPPFLAGS += $(shell root-config --cflags)

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, ${SRCFILES})

OBJSO = libCOAST.so #libCorsikaROOTPlotterInterface.so
EXE = CorsikaROOTPlotter

all: ${OBJSO} ${EXE}

.cc.o:
	${CCC} ${CPPFLAGS} $^

${EXE}: CorsikaROOTPlotter.o TCorsikaPlotter.o TSimpleReader.o
	${CCC}$^ -o $@  ${LDFLAGS} 

lib%.so: CorsikaROOTPlotterInterface.o TCorsikaPlotter.o
	${LD} ${LDFLAGS} --shared $^ ${LIBS} -o $@ 

clean: 
	@rm -f *.o *.a *.so *~ core *.root $(EXE)
