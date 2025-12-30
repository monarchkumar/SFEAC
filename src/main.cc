#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>
#include "version.cc"
// #include "op.cc"
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

    // matrix mat(2);
    // mat.values = {
    //     { 2, 1,-1},
    //     {-3,-1, 2},
    //     {-2, 1, 2}
    // };
    // mat.inverse();
    // mat.output(std::cout);

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

    Grid1D base(0,1,11);
    base.printPoints();
    // base.printSimplex();
    // base.printBoundarySimplices();

    Matrix1D f(11,1);
    cout<<endl;
    f.output(std::cout);


    //Next Step to Program->
    //AssmDelUDelV, will take u matrix and give A matrix. Apply boundary conditions.
    //AssmLin, will take u matrix and give f matrix. Apply boundary conditions.
    //Solution will take A matrix, u matrix, f matrix and calculate u matrix values. Maybe accept RHS and LHS.

    cout << string(NUMBER_OF_DASH, '-') << endl;
    cout << "No Errors. Exiting" << endl;
    cout << string(NUMBER_OF_DASH, '-') << endl;
    return 0;
}


// Old non-generic way to solve equation.

// void possionEquation(){
//     Operator A(10);
//     A.delPhiiDelPhij();
//     A.dirichlet(0,0);
//     A.neumann(9,1);
//     // A.print();

//     Operator f(10);
//     f.Phii(1.0/9);
//     f.dirichlet(0,0);
//     f.neumann(9,1);
//     // f.print();

//     auto invA = A.inverse();
//     // invA.print();

//     auto res = invA.multi(f);
//     res.print();
    
//     ofstream resFile("./output/res.csv");
//     res.opToCSV(resFile);
//     resFile.close();
// }