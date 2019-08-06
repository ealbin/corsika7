CXXFLAGS+=$(shell root-config --cflags)
LDFLAGS+=$(shell root-config --libs)

ifndef COAST_DIR
  COAST_DIR=../../..
endif

CXXFLAGS+=-I$(COAST_DIR)/include

all: libCOAST.so

libCOAST.so: coast.cc
	$(CXX) --shared -fPIC $(CXXFLAGS)  -o $@ $^

clean:
	rm -f *.so *~
