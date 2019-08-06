CCC = g++       # compiler
DICT = rootcint # root dictionary-generator
AR = ar         # archiver
LD = g++        # linker

ifndef COAST_DIR
  COAST_DIR=../../..
endif

ARFLAGS = rcs

LDFLAGS += -fPIC -ggdb3 -Wl,--no-as-needed

LDFLAGS += -L$(COAST_DIR)/lib -L$(COAST_DIR)/lib/unknown
LDFLAGS += -lCorsikaFileIO
LDFLAGS += -lCorsikaIntern

CPPFLAGS += -c -fPIC  -ggdb3
CPPFLAGS += -I$(COAST_DIR)/include

ifdef ROOTSYS
LDFLAGS += $(shell root-config --libs)
CPPFLAGS += $(shell root-config --cflags)
endif

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, ${SRCFILES})

PLOTTER = CorsikaPlotter
READER = CorsikaReader


ifdef ROOTSYS
all: ${PLOTTER} ${READER}

plotter: ${PLOTTER}
else
all: ${READER}
endif
reader: ${PLOTTER}

.cc.o:
	${CCC} ${CPPFLAGS} $^

ifdef ROOTSYS
${PLOTTER}: CorsikaPlotter.o
	${CCC} $^ -o $@ ${LDFLAGS}
endif

${READER}: CorsikaReader.o
	${CCC} $^ -o $@ ${LDFLAGS}


clean: 
	@rm -f *.o *.a *.so *~ core *.root $(EXE)
