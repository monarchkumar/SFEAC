#include <array>
#include <format>
#include <iostream>
#include <map>
#include <ostream>
#include <sys/types.h>
#include <utility>
#include <vector>

class Grid1D{
    public:
    //An indexed x-y point1D
    struct point1D{
        uint index;
        double x;
        bool boundary = false;
    };
    std::vector<point1D> points; //The list of points in 1D grid.

    //An indexed line element
    struct simplexLineSegment{
        uint index;
        std::array<uint, 2> points;
    };
    std::vector<simplexLineSegment> simplices; //The list of simplex line segments in 1D grid.
    std::vector<point1D> boundarySimplices; //The list of boundary points in 1D grid.

    //Make grid when initiated.
    Grid1D(double xMin, double xMax, uint xNum){
        double xLen = xMax-xMin;
        
        for(uint i=0; i<xNum; i++){
            point1D p; 
            p.x=((xLen*i)/(xNum-1)) + xMin;
            p.index = i;
            if(p.x == xMax || p.x == xMin) {p.boundary=true;}
            points.push_back(p);
        }
        calcSimplex(xMax,xNum);
        calcSimplexBoundary(xMax,xNum);
    }

    //Calculate and list simplices
    void calcSimplex(uint xMax, uint xNum){
        //Select point pairs and add to list (barring last point)
        uint index = 0;
        for(point1D p: points){
            if(p.x!=xMax){
                simplices.push_back({index,p.index, p.index+1});
                index++;
            }
        }
    }
    // //Calculate and list boundary simplices
    void calcSimplexBoundary(double xMax, uint xNum){
        for(point1D p: points){
            if(p.boundary){ 
                boundarySimplices.push_back({p.index,p.x});
            }             
        }
    }
    void printPoints(){
        printf("Points:[index-boundary-> x]\n");
        for (point1D p: points) {std::cout<<p.index<<"-"<<p.boundary<<"-> x:"<<p.x<<std::endl;}
    }
    void printSimplex(){
        printf("Line Segment Simplices:[index-> 2x Point Indicies]\n");
        for (simplexLineSegment s: simplices){std::cout<<s.index<<"-> "<<s.points[0]<<" "<<s.points[1]<<std::endl;}    
    }
    void printBoundarySimplices(){
        printf("Boundary Point Simplices:[index-> Point Indicies]\n");
        for (point1D s: boundarySimplices){std::cout<<s.index<<"-> "<<s.x<<std::endl;}    
    }
};

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

class Matrix1D{
    public:
    std::vector<double> values;

    Matrix1D(uint size,double value = 0){
        //Creates a blank matrix with size x size dimension.
        std::vector<double> values(size, value);
        Matrix1D::values = std::move(values);
    }
    void output(std::ostream& output = std::cout){
        //Prints the 'value' variable to the ostream item. Defaults to std::cout.
        std::string row = "";
        for (double val: Matrix1D::values){
            std::string cell =  std::format("{:.3e},",val);
            if(cell[0]!='-') cell = ' ' + cell;
            row+=cell;
        }
        row.pop_back();
        output<<row<<std::endl;
    } 
};

class Matrix2D{
    public:
    bool inverted = false;
    std::vector<std::vector<double>> values;

    Matrix2D(uint size,double value = 0){
        //Creates a blank matrix with size x size dimension.
        std::vector<std::vector<double>> values(size, std::vector<double>(size,value));
        Matrix2D::values = std::move(values);
    }
    void output(std::ostream& output = std::cout){
        //Prints the 'value' variable to the ostream item. Defaults to std::cout.
        for (int i = 0; i<Matrix2D::values.size(); i++){
            std::string row = "";
            for (double val: Matrix2D::values.at(i)){
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
        uint n = Matrix2D::values.size();

        //Augment identity matrix
        for(uint i = 0; i < n; i++){
            std::vector<double> irow(n,0);irow[i]=1;
            Matrix2D::values.at(i).insert(Matrix2D::values.at(i).end(), irow.begin(), irow.end());
        }

        //Row wise elimination
        for (int icol = 0; icol < n; icol++) {

            //Get the row with abs(max) value
            int maxRow = icol;
            for (int irow = icol + 1; irow < n; irow++) {
                if (std::abs(Matrix2D::values[irow][icol]) > std::abs(Matrix2D::values[maxRow][icol])) {
                    maxRow = irow;
                }
            }
            //Put it at the 'icol'
            swap(Matrix2D::values[icol], Matrix2D::values[maxRow]);

            //Say if not invertible.
            if (Matrix2D::values[icol][icol] == 0) {
                std::cout << "Matrix2D is singular." << std::endl;return;  
            }

            //Scale to make matrix::values[icol][icol] = 1
            double scaleFactor = Matrix2D::values[icol][icol];
            for (int jcol = 0; jcol < 2 * n; jcol++) {
                Matrix2D::values[icol][jcol] /= scaleFactor;
            }
            //Scale and substract the row from other rows.
            for (int irow = 0; irow < n; irow++) {
                if (irow != icol) {
                    double factor = Matrix2D::values[irow][icol];
                    for (int jcol = 0; jcol < 2 * n; jcol++) {
                        Matrix2D::values[irow][jcol] -= factor * Matrix2D::values[icol][jcol];
                    }
                }
            }
        }
        //Done, remove the identity matrix
        for(uint i = 0; i < n; i++){
            Matrix2D::values.at(i).erase(Matrix2D::values.at(i).begin(),Matrix2D::values.at(i).begin()+n);
        }
        Matrix2D::inverted= !Matrix2D::inverted;
    }
};


class AssmDelUDelV{
    //Create assembled matrix for Del U . Del V operation.
    public:   
    enum dim {ONEDIM, TWODIM};

        AssmDelUDelV(dim dimension){
            

        };
};