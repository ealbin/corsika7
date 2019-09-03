#include <crs/TSubBlock.h>
#include <crs/MRunHeader.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MParticle.h>

#include <crsRead/MCorsikaReader.h> // for reading
#include <crsRead/TCorsikaReader.h> // for reading
#include <crsRead/TSubBlockIO.h> // for writing

#include <boost/program_options.hpp>

#include <fstream>
#include <string>
#include <iostream>

using namespace std;
namespace po = boost::program_options;


class CMerger : public crsRead::TCorsikaReader {
public:
  CMerger(const std::string &filename, 
          std::ofstream &output, crsRead::TSubBlockIO& outSB,
          const bool writeHeaderBlocks, const bool writeEndBlocks,
          const int verbose) :
    TCorsikaReader(filename, verbose), fOutput(output), fSB(outSB), 
    fWriteHeaderBlocks(writeHeaderBlocks), fWriteEndBlocks(writeEndBlocks) {
  }
  
  virtual void HandleSubBlock(const crs::TSubBlock& sb) { 
    switch(sb.GetBlockType()) {
    case crs::TSubBlock::eRUNH:
    case crs::TSubBlock::eEVTH:
      if (fWriteHeaderBlocks) {
        fSB.SetBuffer(sb);
        fOutput << fSB;
      }
      break;
    case crs::TSubBlock::eRUNE:
    case crs::TSubBlock::eEVTE:
      if (fWriteEndBlocks) {
        fSB.SetBuffer(sb);
        fOutput << fSB;
      }
      break;
    case crs::TSubBlock::ePARTDATA:
      fSB.SetBuffer(sb);
      fOutput << fSB;
    default:
      break;
    }
  }
  
private:
  std::ofstream& fOutput;
  crsRead::TSubBlockIO& fSB;
  bool fWriteHeaderBlocks;
  bool fWriteEndBlocks;
};


void help();

int 
main(int argc, char **argv) 
{ 
  po::options_description desc("Options of the merge_corsika program");
  desc.add_options()
    ("help,h", "description of options")
    ("verbosity,v", po::value<int>(), "verbosity")
    ("input,i", po::value<vector<string> >(), "input file name(s)")
    ("output,o", po::value<string>(), "output file name")
    ;
  
  po::positional_options_description pd;
  pd.add("input", -1); // unlimited

  po::variables_map opt;
  po::store(po::command_line_parser(argc, argv).
            options(desc).positional(pd).run(), opt);
  po::notify(opt);
  
  // check if user needs help-printout
  if (opt.count("help") || 
      (!opt.count("input") && !opt.count("output"))) {
    help();
    exit(1);
  }  

  int verbosity = 0;
  if (opt.count("verbosity") != 0)
    verbosity = opt["verbosity"].as<int>();
  
  vector<string> fileNames;
  fileNames = opt["input"].as<vector<string> >();  
  const string outputName = opt["output"].as<string>();
  
  
  // prepare output file
  std::ofstream binaryFile(outputName.c_str(),
                           std::ios::binary | std::ios::out | std::ios::trunc);
  crsRead::TSubBlockIO subBlock;
  
  
  
  // loop input files
  for(unsigned int iFile = 0; iFile<fileNames.size(); ++iFile) {
    
    const string fname = fileNames[iFile];
    const bool writeHeaderBlocks = iFile == 0;
    const bool writeEndBlocks = iFile == fileNames.size() - 1;
        
    CMerger merger(fname, binaryFile, subBlock, writeHeaderBlocks, writeEndBlocks, verbosity);
    
    // configure binary output sub-block based on first input file
    if (iFile == 0) {
      subBlock.InitSubBlockSize(merger.GetThinning(), merger.GetNBlockSizeInfo());
    }
    
    // read one file and merge it into output
    merger.Read();
    
  } // loop input files
  

  subBlock.FinishBlock(binaryFile); // fill up last block
  binaryFile.flush ();
  binaryFile.close ();
  
  return 0;
}


void 
help() 
{    
  std::cout << " ******************************************************** \n"
	    << " CORSIKA merging tool by Ralf Ulrich, March 3, 2011 \n"
	    << " Karlsruhe Institute of Technologie, ralf.ulrich@kit.edu \n\n"
	    << "    use: ./merge_corsika -i <input-files> -o <output-name> \n"
	    << "    output: one merged file for all input files.\n"
            << "            the first input file will provide the RUN+EVT information\n"
	    << " ******************************************************** \n"
	    << std::endl;
}
