CCC = g++       # compiler
DICT = rootcint # root dictionary-generator
AR = ar         # archiver
LD = g++        # linker

ROOTDIR=$(COAST_DIR)

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
EXE = MakeAnim rootplot

all: $(OBJSO) $(EXE)

lib%.so: HistogramInterface.o TPlotter.o TUser.o TCorsika.o
	$(LD) $(LDFLAGS) --shared $^ $(LIBS) -o $@ 

MakeAnim: MakeAnim.o
	$(LD) $(CPPFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^ 

includes:


clean: 
	@rm -fv *.o *.a *.so *~ core *.root MakeAnim
	@rm -rfv $(OBJSO) $(EXE)
