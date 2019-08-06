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
LDFLAGS += -lboost_program_options

CPPFLAGS += -c -fPIC  -ggdb3 -D_FILE_OFFSET_BITS=64
CPPFLAGS += -I$(COAST_DIR)/include

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, ${SRCFILES})

MERGE = merge_corsika
CHECK = check_corsika


all: ${MERGE} ${CHECK}

.cc.o:
	${CCC} ${CPPFLAGS} $^


${MERGE}: merge_corsika.o
	${CCC} $^ -o $@ ${LDFLAGS}

${CHECK}: check_corsika.o
	${CCC} $^ -o $@ ${LDFLAGS}


clean: 
	@rm -f *.o *.a *.so *~ core *.root $(EXE)
