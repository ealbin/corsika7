CCC = g++       # compiler
DICT = rootcint # root dictionary-generator
AR = ar         # archiver
LD = g++        # linker

ARFLAGS = rcs

LDFLAGS += -fPIC
LDFLAGS += $(shell root-config --libs)

LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaIntern
LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaFileIO

CPPFLAGS += -c -fPIC 
CPPFLAGS += -I$(COAST_DIR)/include
CPPFLAGS += $(shell root-config --cflags)

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, ${SRCFILES})

EXE = CorsikaPlotter

all: ${EXE}

.cc.o:
	${CCC} ${CPPFLAGS} $^

${EXE}: CorsikaPlotter.o
	${CCC} ${LDFLAGS} $^ -o $@

clean: 
	@rm -f *.o *.a *.so *~ core *.root $(EXE)