#include <TSimpleReader.h>
#include <crsRead/TSubBlockIO.h>

#include <iostream>
using namespace std;


int main (int argc, char **argv) {

  if (argc<2) {
    cout << "  please specify Corsika file " << endl;
    return 0;
  }
  
  string fname (argv [1]);
  
  TSimpleReader Reader (fname, 50);
  Reader.Read ();
  
  return 0;
}
