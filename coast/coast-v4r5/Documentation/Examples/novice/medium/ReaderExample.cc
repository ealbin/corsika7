#include <TReaderExample.h>

#include <iostream>
using namespace std;


int main (int argc, char **argv) {

    if (argc<2) {
	cout << "  please specify Corsika file " << endl;
	return 0;
    }

    string fname (argv [1]);

    TReaderExample Reader (fname);
    Reader.Read ();

    return 0;
}
