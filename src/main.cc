#include <fstream>
#include <iostream>
#include "version.cc"
#include "grid.cc"

using namespace std;
#define NUMBER_OF_DASH 64
#define DEBUG true

string version = "0.1";
int main(){
    cout << string(NUMBER_OF_DASH, '-') << endl;
    cout << "SFEAC Version: " << version << endl;

    if(DEBUG){printCppVersion(NUMBER_OF_DASH);}

    Operator A(100);
    A.delPhiiDelPhij();
    A.dirichlet(0,0);
    A.neumann(99,1);
    // A.print();

    Operator f(100);
    f.Phii(1.0/99);
    f.dirichlet(0,0);
    f.neumann(99,1);
    // f.print();

    auto invA = A.inverse();
    // invA.print();

    auto res = invA.multi(f);
    // res.print();
    
    ofstream resFile("./output/res.csv");
    res.opToCSV(resFile);
    resFile.close();

    cout << string(NUMBER_OF_DASH, '-') << endl;
    cout << "No Errors. Exiting" << endl;
    cout << string(NUMBER_OF_DASH, '-') << endl;
    return 0;
}
