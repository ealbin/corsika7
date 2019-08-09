#ifndef _INCLUDE_TCORSIKAPLOTTER__
#define _INCLUDE_TCORSIKAPLOTTER__

#include <string>

namespace crs {
  class TSubBlock;
  class MParticleBlock;
};

class TFile;
class TH1D;


class TCorsikaPlotter {
    
public:
  TCorsikaPlotter ();
  TCorsikaPlotter (const std::string& filename);
  ~TCorsikaPlotter ();
  
  void PlotSubBlock (const crs::TSubBlock &Block);
  
private:
  void PlotParticles (const crs::MParticleBlock &Data);
  
private:
  std::string fFileName;
  TFile *fFile;
  TH1D *fEnergy;
};



#endif
