#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>
#include "version.cc"
#include "grid.cc"

using namespace std;
#define NUMBER_OF_DASH 32
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

    // Grid1D base(0,1,8);

    // Matrix1D f1(base);
    // f1.dirichlet(0, 0);
    // f1.dirichlet(f1.values.size()-1, 0);

    // Matrix2D A1(base);
    // A1.dirichlet(0, 0);
    // A1.dirichlet(f1.values.size()-1, 0);

    // //Doing inv(A)
    // A1.inverse();

    // //Doing inv(A) * f
    // Matrix1D res1 = A1.mul(f1);

    // ofstream file1D("./output/res.csv");
    // base.output(file1D);
    // res1.output(file1D);
    // file1D.close();

    ///////////////////////////////////////////////////
    ///       Solving Poisson equation in 2-D        //
    ///////////////////////////////////////////////////
    Grid2D base2D(0,1,0,1,3,3);
    // base2D.printPoints();
    // base2D.printSimplex();
    // base2D.printBoundarySimplices();

    Matrix1D f2(base2D);
    Matrix2D A2(base2D);

    //Applying Dirichlet Conditions at (0,0)
    //Using u(x,y) = (x^2 + y^2) / 4 = 0
    f2.values[base2D.points[0].index] = 0;
    std::vector<double> temp(A2.values.size(),0);
    A2.values[base2D.points[0].index] = temp;
    A2.values[base2D.points[0].index][base2D.points[0].index] = 1;

    //Applying Neumann Conditions at Boundaries
    //Using u(x,y) = (x^2 + y^2) / 4 
    for (Grid2D::point2D p: base2D.points){
        if(p.boundary){
            f2.values[p.index]={(p.x*p.x + p.y*p.y)/4.0};
        }
    }

    A2.output(cout);
    cout<<endl;
    // f2.output(cout);
    // cout<<endl;

    //Solving
    A2.inverse();
    Matrix1D res2 = A2.mul(f2);

    res2.output();
    cout<<" NOTE: f matrix not defined correctly yet."<<endl;

    ///////////////////////////////////////////////////
    ///                  Next Steps                  //
    ///////////////////////////////////////////////////

    // Fixing 2D Poisson.
    // Adding more flexibility to equations (accessible to more than poisson equations.)
    // Adding better documentation and README


    cout << string(NUMBER_OF_DASH, '-') << endl;
    cout << "No Errors. Exiting" << endl;
    cout << string(NUMBER_OF_DASH, '-') << endl;
    return 0;
}
