CCC = g++       # compiler
DICT = rootcint # root dictionary-generator
AR = ar         # archiver
LD = g++        # linker

ARFLAGS = rcs

LDFLAGS += -Wl,--no-as-needed -fPIC -ggdb3
LDFLAGS += $(shell root-config --libs)

LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaIntern
LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaFileIO

CXXFLAGS += -c -fPIC  -ggdb3
CXXFLAGS += -I$(COAST_DIR)/include
CXXFLAGS += $(shell root-config --cflags)

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, ${SRCFILES})

EXE = CorsikaPlotter # CorsikaDumper

all: ${EXE}

CorsikaPlotter: CorsikaPlotter.o
# CorsikaDumper: CorsikaDumper.o

clean: 
	@rm -f *.o *.a *.so *~ core *.root $(EXE)