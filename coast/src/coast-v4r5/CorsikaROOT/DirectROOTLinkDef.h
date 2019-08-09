#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;           // added for namespace
#pragma link C++ nestedtypedef;         // added for namespace

#pragma link C++ namespace crsIO;       // added for namespace

#pragma link C++ class crsIO::TParticle+;
#pragma link C++ class crsIO::TLongitudinal+;
#pragma link C++ class crsIO::TCherenkov+;
#pragma link C++ class crsIO::TShower+;
#pragma link C++ class crsIO::TRun+;

#endif

