#include <fstream>
#include <iostream>
#include <ostream>
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

    ///////////////////////////////////////////////////
    ///                 Example code                 //
    ///////////////////////////////////////////////////

    // Matrix2D mat(3);
    // mat.values = {
    //     { 2, 1,-1},
    //     {-3,-1, 2},
    //     {-2, 1, 2}
    // };
    // mat.inverse();
    // mat.output();

    // ofstream outputfile("./output/mat.csv");
    // mat.output(outputfile);
    // outputfile.close();

    // Grid1D grid(0,1,11);
    // grid.printPoints();
    // grid.printSimplex();
    // grid.printBoundarySimplices();

    // Grid2D grid2d(0,1,0,1,3,3);
    // grid2d.printPoints();
    // grid2d.printSimplex();
    // grid2d.printBoundarySimplices();

    ///////////////////////////////////////////////////
    ///       Solving Poisson equation in 1-D        //
    ///////////////////////////////////////////////////

    Grid1D base(0,1,10);

    Matrix1D f(base);
    f.dirichlet(0, -1);
    f.dirichlet(f.values.size()-1, 1.5);

    Matrix2D A(base);
    A.dirichlet(0, -1);
    A.dirichlet(f.values.size()-1, 1.5);

    //Doing inv(A)
    A.inverse();

    //Doing inv(A) * f
    Matrix1D res = A.mul(f);

    ofstream file("./output/res.csv");
    res.output(file);
    file.close();

    ///////////////////////////////////////////////////
    ///                  Next Steps                  //
    ///////////////////////////////////////////////////

    //Make adjustments for 2D scenarios

    cout << string(NUMBER_OF_DASH, '-') << endl;
    cout << "No Errors. Exiting" << endl;
    cout << string(NUMBER_OF_DASH, '-') << endl;
    return 0;
}
