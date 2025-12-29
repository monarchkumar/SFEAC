#include <array>
#include <format>
#include <iostream>
#include <istream>
#include <map>
#include <ostream>
#include <sys/types.h>
#include <utility>
#include <vector>


#include <fstream>

class Grid2D{
    public:
        //An indexed x-y point2D
        struct point2D{
            uint index;
            double x;
            double y;
            bool boundary = false;
        };
        std::vector<point2D> points; //The list of points in 2D grid.

        //An indexed triangle element
        struct simplexTriangle{
            uint index;
            std::array<uint, 3> points;
        };
        struct simplexLineSegment{
            uint index;
            std::array<uint, 2> points;
        };
        std::vector<simplexTriangle> simplices; //The list of simplex triangles in 2D grid.
        std::vector<simplexLineSegment> boundarySimplices; //The list of boundary simplex lines in 2D grid.

        //Make grid when initiated.
        Grid2D(double xMin, double xMax, double yMin, double yMax, uint xNum, uint yNum){
            double xLen = xMax-xMin;
            double yLen = yMax-yMin;
            
            for(uint i=0; i<xNum; i++){
                for(uint j=0; j<yNum; j++){
                    point2D p; 
                    p.x=((xLen*i)/(xNum-1)) + xMin;
                    p.y=((yLen*j)/(yNum-1)) + yMin;
                    p.index = j + xNum * i;
                    if(p.x == xMax || p.x == xMin || p.y == yMax || p.y == yMin) {p.boundary=true;}
                    points.push_back(p);
                }
            }
            calcSimplex(xMax,yMax,xNum,yNum);
            calcSimplexBoundary(xMax,yMax,xNum,yNum);
        }

        //Calculate and list simplices
        void calcSimplex(double xMax, double yMax, uint xNum, uint yNum){
            //If point2D's x or y at max do nothing; else, make two triangles point2D's x->x+1 and y-y+1 region.
            uint index = 0;
            for(point2D p: points){
                if(p.x!=xMax && p.y!=yMax){
                    simplices.push_back({index,p.index, p.index+yNum, p.index+1});
                    index++;
                    simplices.push_back({index,p.index+1+yNum, p.index+1, p.index+yNum});
                    index++;
                }
            }
        }
        //Calculate and list boundary simplices
        void calcSimplexBoundary(double xMax, double yMax, uint xNum, uint yNum){
            for(point2D p: points){
                if(p.boundary){ //If point2D is boundary, look x+1 and y+1 items; if both boundary add pair to boundarySimplices
                    if(p.y != yMax) {
                        if(p.boundary == points.at(p.index+1).boundary){
                            boundarySimplices.push_back({static_cast<uint>(boundarySimplices.size()),p.index, p.index+1});
                        }
                    }
                    if(p.x != xMax) {
                        if(p.boundary == points.at(p.index+yNum).boundary){
                            boundarySimplices.push_back({static_cast<uint>(boundarySimplices.size()),p.index, p.index+yNum});
                        }
                    }
                }             
            }
        }
        void printPoints(){
            printf("Points:[index-boundary-> x , y]\n");
            for (point2D p: points) {std::cout<<p.index<<"-"<<p.boundary<<"-> x:"<<p.x<<", y:"<<p.y<<std::endl;}
        }
        void printSimplex(){
            printf("Triangular Simplices:[index-> 3x Point Indicies]\n");
            for (simplexTriangle s: simplices){std::cout<<s.index<<"-> "<<s.points[0]<<" "<<s.points[1]<<" "<<s.points[2]<<std::endl;}    
        }
        void printBoundarySimplices(){
            printf("Boundary Line Segment Simplices:[index-> 2x Point Indicies]\n");
            for (simplexLineSegment s: boundarySimplices){std::cout<<s.index<<"-> "<<s.points[0]<<" "<<s.points[1]<<std::endl;}    
        }
};

class Assembly1D{

    public:   
        std::map<std::pair<int, int>,double> matrix;

        Assembly1D(auto kernel, uint size, uint offset){
            double h=1.0/(size-1);

            for (uint o = 0; o<size-1; o+=offset){
                for (uint i = 0; o<kernel.size(); i+=1 ){
                    for (uint j = 0; j<kernel[0].size(); j+=1 ){
                        matrix[{i+o,j+o}]= kernel[i][j];
                    }
                }
            }
        }
        void print(){
            for(auto x:matrix){
                std::cout<<x.first.first<<" "<<x.first.second<<" : "<<x.second<<std::endl;
            }
        }
};
// void delPhiiDelPhij(){
//     op.resize(size);
//     for (int i = 0; i<size; i++){
//         op[i].resize(size);
//         op[i][i]=2;
//         if (i==0){op[i][i+1]=-1;}
//         if (i==op.size()-1){op[i][i-1]=-1;}
//         if (i>0 && i<op.size()-1) {
//             op[i][i-1]=-1;
//             op[i][i+1]=-1;
//         }
//         for(double& i : op[i]) { i *= (1/h);}
//     }
//     op[0][0]/=2;
//     op[size-1][size-1]/=2;
//     defined = true;
// }

class matrix{
    public:
    bool inverted = false;
    std::vector<std::vector<double>> values;
    matrix(uint size){
        //Creates a blank matrix with size x size dimension.
        std::vector<std::vector<double>> values(size, std::vector<double>(size,0));
        matrix::values = std::move(values);
    }
    void output(std::ostream& output = std::cout){
        //Prints the 'value' variable to the ostream item. Defaults to std::cout.
        for (int i = 0; i<matrix::values.size(); i++){
            std::string row = "";
            for (double val: matrix::values.at(i)){
                std::string cell =  std::format("{:.3e},",val);
                if(cell[0]!='-') cell = ' ' + cell;
                row+=cell;
            }
            row.pop_back();
            output<<row<<std::endl;
        }
    } 
    // https://www.math-cs.gordon.edu/courses/ma342/handouts/gauss.pdf
    void inverse(){
        uint n = matrix::values.size();

        //Augment identity matrix
        for(uint i = 0; i < n; i++){
            std::vector<double> irow(n,0);irow[i]=1;
            matrix::values.at(i).insert(matrix::values.at(i).end(), irow.begin(), irow.end());
        }

        //Row wise elimination
        for (int icol = 0; icol < n; icol++) {

            //Get the row with abs(max) value
            int maxRow = icol;
            for (int irow = icol + 1; irow < n; irow++) {
                if (std::abs(matrix::values[irow][icol]) > std::abs(matrix::values[maxRow][icol])) {
                    maxRow = irow;
                }
            }
            //Put it at the 'icol'
            swap(matrix::values[icol], matrix::values[maxRow]);

            //Say if not invertible.
            if (matrix::values[icol][icol] == 0) {
                std::cout << "Matrix is singular." << std::endl;return;  
            }

            //Scale to make matrix::values[icol][icol] = 1
            double scaleFactor = matrix::values[icol][icol];
            for (int jcol = 0; jcol < 2 * n; jcol++) {
                matrix::values[icol][jcol] /= scaleFactor;
            }
            //Scale and substract the row from other rows.
            for (int irow = 0; irow < n; irow++) {
                if (irow != icol) {
                    double factor = matrix::values[irow][icol];
                    for (int jcol = 0; jcol < 2 * n; jcol++) {
                        matrix::values[irow][jcol] -= factor * matrix::values[icol][jcol];
                    }
                }
            }
        }
        //Done, remove the identity matrix
        for(uint i = 0; i < n; i++){
            matrix::values.at(i).erase(matrix::values.at(i).begin(),matrix::values.at(i).begin()+n);
        }
        inverted= !inverted;
    }
};


