CCC = g++       # compiler
DICT = rootcint # root dictionary-generator
AR = ar         # archiver
LD = g++        # linker

ifdef COAST_DIR
  ROOTDIR=$(COAST_DIR)
else
  ROOTDIR=..
endif

ARFLAGS = rcs

LDFLAGS += -Wl,--no-as-needed -fPIC 
LDFLAGS += $(shell root-config --libs)

LDFLAGS += -L$(ROOTDIR)/lib -lCorsikaIntern
LDFLAGS += -L$(ROOTDIR)/lib -lCorsikaROOT
LDFLAGS += -L$(ROOTDIR)/lib -lCorsikaFileIO
LDFLAGS += -L$(ROOTDIR)/lib -lCorsikaToROOT

CPPFLAGS += -fPIC -ggdb3
CPPFLAGS += -I. -I$(ROOTDIR)/include
CPPFLAGS += $(shell root-config --cflags)

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, $(SRCFILES))

OBJSO = libCOAST.so

all: $(OBJSO) Plot2DTest Plot2DGround

.cc.o:
	$(CCC) -c $(CPPFLAGS) $^

lib%.so: Plot2DInterface.o TPlotter.o
	$(LD) $(LDFLAGS) -Wl,--rpath=$(shell root-config --libdir) -lASImage --shared $^ $(LIBS) -o $@ 

Plot2DTest: Plot2DTest.cc TPlotter.o
	@echo "========> You need ROOT's libASImage for this program"
	$(CCC)  $(CPPFLAGS) $(LDFLAGS)  -lASImage $^ -o $@

Plot2DGround: Plot2DGround.cc TPlotter.o
	@echo "========> You need ROOT's libASImage for this program"
	$(CCC)  $(CPPFLAGS) $(LDFLAGS)  -lASImage $^ -o $@

includes:


clean: 
	@rm -f *.o *.a *.so *~ core *.root Plot2DTest Plot2DGround
	@rm -rf $(OBJSO)
