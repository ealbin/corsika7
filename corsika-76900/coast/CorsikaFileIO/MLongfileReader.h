#ifndef _M_LONGFILEREADER_H
#define _M_LONGFILEREADER_H

#include <string>
#include <istream>
#include <vector>



namespace crsRead {
  
  class TSubBlock;
  
  /** 
      \class MLongfileReader
      \brief Implementation of a CORSIKA longitudinal file reader
      
      An easy to use class for CORSIKA long file reading
      
      \author Ralf Ulrich
      \date Fr 20. Jul 15:05:26 CEST 2007
      \version $Id: MLongfileReader.h 5116 2016-01-04 19:09:04Z darko $ 
  */
  
  class MLongfileReader {
      
  public:
      struct ParticleEntry {
	  ParticleEntry(double x, 
			double ng, 
			double np, 
			double ne,
			double nmb,
			double nm,
			double nh,
			double nc,
			double nn,
			double nckov) 
	      : X(x), 
		Ngamma(ng),
		Npositron(np),
		Nelectron(ne),
		Nmu_bar(nmb),
		Nmu(nm),
		Nhadrons(nh),
		Ncharged(nc),
		Nnuclei(nn),
		Ncherenkov(nckov) {}
	  double X;
	  double Ngamma;
	  double Npositron;
	  double Nelectron;
	  double Nmu_bar;
	  double Nmu;
	  double Nhadrons;
	  double Ncharged;
	  double Nnuclei;
	  double Ncherenkov;
      };

      struct EdepEntry {
	  EdepEntry(double x, 
		    double g, 
		    double em, 
		    double emCut,
		    double mu,
		    double muCut,
		    double had,
		    double hadCut,
		    double neut,
		    double sum) 
	      : X(x), 
		dEdXgamma(g),
		dEdXemIoniz(em),
		dEdXemCut(emCut),
		dEdXmuIoniz(mu),
		dEdXmuCut(muCut),
		dEdXhadIoniz(had),
		dEdXhadCut(hadCut),
		dEdXneutrino(neut),
		dEdXsum(sum) {}
	  double X;
	  double dEdXgamma;
	  double dEdXemIoniz;
	  double dEdXemCut;
	  double dEdXmuIoniz;
	  double dEdXmuCut;
	  double dEdXhadIoniz;
	  double dEdXhadCut;
	  double dEdXneutrino;
	  double dEdXsum;
      };
      
  public:
    MLongfileReader (const std::string &filename, int verbose=0);
    
    bool ReadShower(int n);
    
    // gaisser hillas parameters
    double GetXmax() const {return fXmax;}
    double GetX0() const {return fX0;}
    double GetNmax() const {return fNmax;}
    double GetA() const {return fA;}
    double GetB() const {return fB;}
    double GetC() const {return fC;}
    double GetChi2() const {return fChi2;}
    
    // particle profiles
    bool IsSlant() const {return fIsSlant;}
    bool IsVertical() const {return !fIsSlant;}
    double GetNBins() const {return fNBins;}
    double GetDeltaX() const {return fDeltaX;}
    const ParticleEntry& GetEntry(int i) const {return fParticleEntry[i];}
    
    // energy deposit profiles
    bool IsSlantEdep() const {return fIsSlantEdep;}
    bool IsVerticalEdep() const {return !fIsSlantEdep;}
    double GetNBinsEdep() const {return fNBinsEdep;}
    double GetDeltaXEdep() const {return fDeltaXEdep;}
    const EdepEntry& GetEntryEdep(int i) const {return fEdepEntry[i];}

  private:
    bool FindShower(int n);
    bool ParseParticleTable(int n);
    bool ParseEdepTable(int n);
    
  private:
    std::string fFileName;
    std::ifstream* fFile;
    int fVerbose;
    
    double fXmax;
    double fX0;
    double fNmax;
    double fA;
    double fB;
    double fC;
    double fChi2;
    
    bool fIsSlant;
    double fNBins;
    double fDeltaX;
    std::vector<ParticleEntry> fParticleEntry;
    
    bool fIsSlantEdep;
    double fNBinsEdep;
    double fDeltaXEdep;
    std::vector<EdepEntry> fEdepEntry;
    
  };
  
};

#endif
