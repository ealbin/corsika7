CCC = g++       # compiler
DICT = rootcint # root dictionary-generator
AR = ar         # archiver
LD = g++        # linker

ifdef COAST_DIR
  ROOTDIR=$(COAST_DIR)
else
  exit
endif

ARFLAGS = rcs

LDFLAGS += -Wl,--no-as-needed -fPIC  -ggdb3
LDFLAGS += $(shell root-config --libs)

LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaIntern
LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaROOT
LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaFileIO
LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaToROOT

CPPFLAGS += -c -fPIC -ggdb3 -Wall
CPPFLAGS += -I. -I$(ROOTDIR)/include
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
