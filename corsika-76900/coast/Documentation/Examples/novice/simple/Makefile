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
LDFLAGS += -lCorsikaFileIO
LDFLAGS += -lCorsikaIntern
LDFLAGS += $(shell root-config --libs)

CPPFLAGS += -c -fPIC  -ggdb3
CPPFLAGS += -I$(COAST_DIR)/include
CPPFLAGS += $(shell root-config --cflags)

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, ${SRCFILES})

EXE = CorsikaPlotter

all: ${EXE}

.cc.o:
	${CCC} ${CPPFLAGS} $^

${EXE}: CorsikaPlotter.o
	${CCC} $^ -o $@ ${LDFLAGS} 

clean: 
	@rm -f *.o *.a *.so *~ core *.root $(EXE)
