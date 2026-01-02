#include <algorithm>
#include <array>
#include <cstddef>
#include <format>
#include <iostream>
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
        double length;
    };
    std::vector<simplexLineSegment> simplices; //The list of simplex line segments in 1D grid.
    std::vector<point1D> boundarySimplices; //The list of boundary points in 1D grid.

    //Make grid when initiated. Uniform distribution.
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
                simplices.push_back({index,p.index,p.index+1,points[p.index+1].x-p.x});
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
        printf("Line Segment Simplices:[index-> 2x Point Indicies, length]\n");
        for (simplexLineSegment s: simplices){std::cout<<s.index<<"-> "<<s.points[0]<<" "<<s.points[1]<<" len:"<<s.length<<std::endl;}    
    }
    void printBoundarySimplices(){
        printf("Boundary Point Simplices:[index-> x]\n");
        for (point1D s: boundarySimplices){std::cout<<s.index<<"-> x:"<<s.x<<std::endl;}    
    }
    void output(std::ostream& output = std::cout){
        //Prints the 'value' variable to the ostream item. Defaults to std::cout.
        std::string row = "";
        for (point1D p: Grid1D::points){
            std::string cell =  std::format("{:.3e},",p.x);
            if(cell[0]!='-') cell = ' ' + cell;
            row+=cell;
        }
        row.pop_back();
        output<<row<<std::endl;
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
        double area;
    };
    struct simplexLineSegment{
        uint index;
        std::array<uint, 2> points;
        double length;
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
        std::array<uint,3> set;
        for(point2D p: points){
            if(p.x!=xMax && p.y!=yMax){
                simplices.push_back({index,p.index, p.index+yNum, p.index+1});
                set = simplices.back().points;
                simplices.back().area = 0.5 * std::abs(
                 points[set[0]].x * (points[set[1]].y - points[set[2]].y) +
                    points[set[1]].x * (points[set[2]].y - points[set[0]].y) +
                    points[set[2]].x * (points[set[0]].y - points[set[1]].y)
                ); 
                index++;

                simplices.push_back({index,p.index+yNum+1, p.index+1, p.index+yNum});
                set = simplices.back().points;
                simplices.back().area = 0.5 * std::abs(
                 points[set[0]].x * (points[set[1]].y - points[set[2]].y) +
                    points[set[1]].x * (points[set[2]].y - points[set[0]].y) +
                    points[set[2]].x * (points[set[0]].y - points[set[1]].y)
                ); 
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
                        boundarySimplices.push_back({static_cast<uint>(boundarySimplices.size()),p.index, p.index+1,std::abs(points[p.index+1].y-p.y)});
                    }
                }
                if(p.x != xMax) {
                    if(p.boundary == points.at(p.index+yNum).boundary){
                        boundarySimplices.push_back({static_cast<uint>(boundarySimplices.size()),p.index, p.index+yNum,std::abs(points[p.index+yNum].x-p.x)});
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
        printf("Triangular Simplices:[index-> 3x Point Indicies, Area]\n");
        for (simplexTriangle s: simplices){std::cout<<s.index<<"-> "<<s.points[0]<<" "<<s.points[1]<<" "<<s.points[2]<<", "<<s.area<<std::endl;}    
    }
    void printBoundarySimplices(){
        printf("Boundary Line Segment Simplices:[index-> 2x Point Indicies, length]\n");
        for (simplexLineSegment s: boundarySimplices){std::cout<<s.index<<"-> "<<s.points[0]<<" "<<s.points[1]<<" len:"<<s.length<<std::endl;}    
    }
    void output(std::ostream& output = std::cout){
        //Prints the 'points' variable to the ostream item. Defaults to std::cout.
        for (int i = 0; i<Grid2D::points.size(); i++){
            std::vector<std::string> rows(3);
            for (auto val: Grid2D::points){
                std::string index =  std::format("{},",val.index);
                rows[0] = std::move(index);
                std::string x =  std::format("{:.3e},",val.x);
                if(x[0]!='-') x = ' ' + x;
                rows[1] = std::move(x);
                std::string y =  std::format("{:.3e},",val.y);
                if(y[0]!='-') y = ' ' + y;
                rows[2] = std::move(y);
            }
            rows[1].pop_back();
            rows[2].pop_back();
            for(auto row:rows){output<<row<<std::endl;}
        }
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
    Matrix1D(Grid1D grid1d){
        //Creates a matrix and calculates values based on slopes in 'grid' given. Indexes values according to 'grid' indices.  
        std::vector<double> values(grid1d.points.size(), 0);
        for (Grid1D::simplexLineSegment i: grid1d.simplices){
            values[i.points.front()]+=0.5*i.length;
            values[i.points.back()]+=0.5*i.length;
        }
        Matrix1D::values = std::move(values);
    }
    //Finite elements : theory, fast solvers, and applications in solid mechanics by Braess, Dietrich, 2001, pg 54
    //NOTE: Needs work
    Matrix1D(Grid2D grid2d){
        for(Grid2D::point2D p:grid2d.points){
            Matrix1D::values.push_back(0);
        }
    }
    void dirichlet(uint index, double value){values[index]=value;}
    void neumann(uint index, double value){values[index]+=value;}
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
    Matrix2D(Grid2D grid2d){
        // Creates a matrix and calculates values based on slopes in 'grid' given. Indexes values according to 'grid' indices. 
        // Procedure:
        // Loop through simplexes.
        // Calculate gradients (ui/dx , ui/dy) for each point (i) in the simplex.
        // Then we should have all variables for required equation LHS calculation.

        std::vector<std::vector<double>> values(grid2d.points.size(), std::vector<double>(grid2d.points.size(),0));
        for (Grid2D::simplexTriangle s: grid2d.simplices){

            Grid2D::point2D p0 = grid2d.points[s.points[0]];
            Grid2D::point2D p1 = grid2d.points[s.points[1]];
            Grid2D::point2D p2 = grid2d.points[s.points[2]];
            // slope in x = a & slope in y = b when we calculate
            // u1   x1 y1 1   a
            // u2 = x2 y2 1 X b
            // u3   x3 y3 1   c

            //ui values will change each point in simplex, will define later
            Matrix2D pij(3);
            pij.values={
                {p0.x, p0.y, 1},
                {p1.x, p1.y, 1},
                {p2.x, p2.y, 1},
            };
            pij.inverse();

            //Example
            //ui.values = {0,1,0};
            //std::vector<double> slope = pij.mul(ui).values;
            //We get {del ui / del x, del ui / del y, dontcare} of simplex.

            std::vector<std::vector<double>> slope;
            for(uint idx =0; idx < s.points.size(); idx++){
                // Get slope when each point is basis
                Matrix1D ui(3);
                ui.values[idx] = 1;
                std::vector<double> slopeIdx = pij.mul(ui).values;
                slope.push_back(slopeIdx);
            }
            for(uint i = 0;i <s.points.size();i++){
                for(uint j = 0;j <s.points.size();j++){
                    // OPERATION FOR POISSON 2D//
                    // For Poisson we need, Int Int {(del ui / del x) + (del uj / del x)} * {(del ui / del j) + (del uj / del y)} dx dy
                    values[s.points[i]][s.points[j]] += ((slope[i][0] * slope[j][0]) + (slope[i][1] * slope[j][1])) * (s.area);
                }
            }
        }
        Matrix2D::values = std::move(values);
    }
    Matrix2D(Grid1D grid){
        //Creates a matrix and calculates values based on slopes in 'grid' given. Indexes values according to 'grid' indices.  
        //[ 1,-1] * 1/h
        //[-1, 1] * 1/h
        std::vector<std::vector<double>> values(grid.points.size(), std::vector<double>(grid.points.size(),0));
        for (Grid1D::simplexLineSegment i: grid.simplices){
            double h = grid.points[i.points.back()].x-grid.points[i.points.front()].x;
            values[i.points.front()][i.points.front()] += 1.0/h;
            values[i.points.front()][i.points.back()] -= 1.0/h;
            values[i.points.back()][i.points.front()] -= 1.0/h;
            values[i.points.back()][i.points.back()] += 1.0/h;
        }
        Matrix2D::values = std::move(values);
    }
    Matrix1D mul(Matrix1D mat){
        uint size = mat.values.size();
        Matrix1D res(size,0);
        for(uint i = 0;i<size;i++){
            for(uint j = 0;j<size; j++){
                res.values[i]+=Matrix2D::values[i][j]*mat.values[j]; 
            }
        }
        return res;
    }
    void boundary(void (*boundaryAssignment)(Grid2D grid)){

    }
    void dirichlet(uint index, double value=0){
        std::vector<double> row(Matrix2D::values.size(),0);
        row[index] = 1;
        Matrix2D::values[index]=std::move(row);
    }
    void neumann(uint index=0, double value=0){
        NULL;
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
    void inverse(){
    // https://www.math-cs.gordon.edu/courses/ma342/handouts/gauss.pdf
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


