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

LDFLAGS += -Wl,--no-as-needed -fPIC -ggdb3
LDFLAGS += $(shell root-config --libs)

LDFLAGS += -L$(ROOTDIR)/lib -lCorsikaIntern
LDFLAGS += -L$(ROOTDIR)/lib -lCorsikaROOT
LDFLAGS += -L$(ROOTDIR)/lib -lCorsikaFileIO
LDFLAGS += -L$(ROOTDIR)/lib -lCorsikaToROOT

CPPFLAGS += -fPIC  -ggdb3 -Wall
CPPFLAGS += -I. -I$(ROOTDIR)/include
CPPFLAGS += $(shell root-config --cflags)

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, $(SRCFILES))

OBJSO = libCOAST.so
EXE = 

all: $(OBJSO) $(EXE)


lib%.so: InclinedPlane.o TPlotter.o 
	$(LD) $(LDFLAGS) --shared $^ $(LIBS) -o $@ 


clean: 
	@rm -fv *.o *.a *.so *~ core *.root MakeAnim Makefile.dep
	@rm -rfv $(OBJSO) $(EXE)

