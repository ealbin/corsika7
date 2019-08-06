/***************************************************************************
                          sumparallelresults.cpp  -  description
                             -------------------
    begin                : Apr 11 2014
    copyright            : (C) 2014 by Tim Huege
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *  Copyright and any other appropriate legal protection of these          *
 *  computer programs and associated documentation are reserved in         *
 *  all countries of the world.                                            *
 *                                                                         *
 *  These programs or documentation may not be reproduced by any method    *
 *  without prior written consent of Karlsruhe Institute of Technology or its    *
 *  delegate. Commercial and military use are explicitly forbidden.        *
 *                                                                         *
 *  Karlsruhe Institute of Technology welcomes comments concerning the REAS      *
 *  code but undertakes no obligation for maintenance of the programs,     *
 *  nor responsibility for their correctness, and accepts no liability     *
 *  whatsoever resulting from the use of its programs.                     *
 *                                                                         *
 ***************************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <string>
#include <iomanip>
#include "threevector.h"

#include <sys/stat.h> // for mkdir
#include <sys/types.h> // for mkdir
#include <errno.h>  // for mkdir


using namespace std;

int main(int argc, char *argv[])
{
  map<string, map<double,ThreeVector> > results;

  // no arguments on command line
  if (argc==1) {
	cout << "\nSyntax: sumparallelresults *.bins\n" << endl;
    return 10;
  }

  const string dirname = string(argv[1]).substr(0,9)+"_coreas";

  // create results directory
  errno = 0;
  int result = mkdir(dirname.c_str(), 0777);  // mkdir resultdir
  if ( (result != 0) && (errno != EEXIST) )
  {
    cout << "Error creating or accessing directory '" << argv[1] << "'." << endl;
    throw 2;
  }

  for (int i=1; i<argc; ++i)
  {
    cout << "Processing " << argv[i] << " ...\n";
    ifstream binfile(argv[i]);
    string filename;
    string binname;
    double dummy;
  
    while (binfile) {
      filename.clear();
      binfile >> filename >> dummy >> dummy >> dummy >> dummy >> dummy;

      unsigned int pos = filename.find("_coreas/",0);
      if (pos >= filename.size())
        continue;

//    cout << pos << "\t" << filename << "\t" << filename.size() << endl;

      binname = filename.substr(pos+8,filename.size()-pos-8);
//    cout << binname << endl;

      double time;
      ThreeVector efield;

      ifstream datfile(filename.c_str());
      while (datfile) {
        datfile >> time >> efield;
        if (datfile) {
//        cout << time << "\t" << efield << endl;
          results[binname][time]+=efield;
	    }
      }
      datfile.close();
    }
    binfile.close();
  }


  
  for (map<string, map<double,ThreeVector> >::const_iterator it=results.begin(); it!=results.end(); ++it) {
    ofstream resfile((dirname+"/"+it->first).c_str());
    for (map<double,ThreeVector>::const_iterator sp=it->second.begin(); sp!=it->second.end(); ++sp) {
      resfile << sp->first << "\t" << setprecision(14) << sp->second << "\n";
    }
    resfile.close();
  }

  std::cout << "Done. Created summed up simulation files in directory " << dirname << ". Use it with any of the .bins files of the parallel simulation.\n";
  
}
