#include <fstream>
#include <iostream>
#include <ostream>
#include "version.cc"
#include "op.cc"
#include "grid.cc"

using namespace std;
#define NUMBER_OF_DASH 64
#define DEBUG true

string version = "0.1";

void possionEquation(){
    Operator A(10);
    A.delPhiiDelPhij();
    A.dirichlet(0,0);
    A.neumann(9,1);
    // A.print();

    Operator f(10);
    f.Phii(1.0/9);
    f.dirichlet(0,0);
    f.neumann(9,1);
    // f.print();

    auto invA = A.inverse();
    // invA.print();

    auto res = invA.multi(f);
    res.print();
    
    ofstream resFile("./output/res.csv");
    res.opToCSV(resFile);
    resFile.close();
}

int main(){
    cout << string(NUMBER_OF_DASH, '-') << endl;
    cout << "SFEAC Version: " << version << endl;

    if(DEBUG){printCppVersion(NUMBER_OF_DASH);}

    // possionEquation();


    // Grid2D grid2d(0,1,0,1,3,3);
    // grid2d.printPoints();
    // cout<<grid2d.points.size()<<endl;
    // cout<<grid2d.simplices.size()<<endl;
    // cout<<grid2d.boundarySimplices.size()<<endl;
    // grid2d.printSimplex();
    // grid2d.printBoundarySimplices();

    // kernel for delta testfun_i() *  delta testfun_j()
    // std::array<std::array<double, 2>, 2> ker = {{
    //     {1,-1},
    //     {-1,1}
    // }};
    // Assembly1D a(ker,5,1);
    // a.print();

    matrix mat(2);
    mat.values = {
        { 2, 1,-1},
        {-3,-1, 2},
        {-2, 1, 2}
    };
    mat.inverse();
    mat.output(std::cout);

    ofstream outputfile("./output/mat.csv");
    mat.output(outputfile);
    outputfile.close();




    cout << string(NUMBER_OF_DASH, '-') << endl;
    cout << "No Errors. Exiting" << endl;
    cout << string(NUMBER_OF_DASH, '-') << endl;
    return 0;
}
