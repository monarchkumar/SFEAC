#include <fstream>
#include <iostream>
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

    /*
    mesh.tri.points #[x,y] pair of points in the grid.
    mesh.tri.nsimplex #Number of simplex (triangles)
    mesh.tri.npoints #Number of points
    mesh.tri.simplices # List of "index" of points which make a simplex. 3points x numOfSimplex
    mesh.tri.convex_hull # the list of outer most edges (edges so vertex pair)
    */
    Grid2D grid2d(0,1,0,1,3,3);
    // grid2d.printPoints();
    // cout<<grid2d.simplices.size()<<endl;
    // cout<<grid2d.points.size()<<endl;
    // grid2d.printSimplex();



    cout << string(NUMBER_OF_DASH, '-') << endl;
    cout << "No Errors. Exiting" << endl;
    cout << string(NUMBER_OF_DASH, '-') << endl;
    return 0;
}
