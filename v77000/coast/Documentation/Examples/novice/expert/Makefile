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

CPPFLAGS += -c -fPIC -ggdb3
CPPFLAGS += -I. -I$(COAST_DIR)/include
CPPFLAGS += $(shell root-config --cflags)

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, ${SRCFILES})

OBJSO = libCOAST.so #libExpertInterface.so

EXE = MakeAnim MakeAnim2D CompareAnim

all: ${EXE} ${OBJSO}

CompareAnim: CompareAnim.o
	${CCC} $^ ${LDFLAGS} -o $@

MakeAnim: MakeAnim.o
	${CCC} $^ ${LDFLAGS} -o $@

MakeAnim2D: MakeAnim2D.o
	${CCC} $^ -o $@ ${LDFLAGS} 

.cc.o:
	${CCC} ${CPPFLAGS} $^

lib%.so: ExpertInterface.o TPlotter.o
	${LD} ${LDFLAGS}  --shared $^ ${LIBS} -o $@ 

clean: 
	@rm -f *.o *.a *.so *~ core *.root ${EXE}
