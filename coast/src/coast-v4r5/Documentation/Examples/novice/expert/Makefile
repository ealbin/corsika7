CCC = g++       # compiler
DICT = rootcint # root dictionary-generator
AR = ar         # archiver
LD = g++        # linker

ARFLAGS = rcs

LDFLAGS += -fPIC 
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

OBJSO = libCOAST.so #libExpertInterface.so

EXE = MakeAnim MakeAnim2D CompareAnim

all: ${EXE} ${OBJSO}

CompareAnim: CompareAnim.o
	${CCC} $^ ${LDFLAGS} -o $@

MakeAnim: MakeAnim.o
	${CCC} $^ ${LDFLAGS} -o $@

MakeAnim2D: MakeAnim2D.o
	${CCC} $^ ${LDFLAGS} -o $@

.cc.o:
	${CCC} ${CPPFLAGS} $^

lib%.so: ExpertInterface.o TPlotter.o
	${LD} ${LDFLAGS} --shared $^ ${LIBS} -o $@ 

clean: 
	@rm -f *.o *.a *.so *~ core *.root ${EXE}