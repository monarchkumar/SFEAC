#include <array>
#include <iostream>
#include <vector>

class Grid2D{
    public:
        //An indexed x-y point
        struct point{
            uint index;
            double x;
            double y;
            bool boundary = false;
        };
        std::vector<point> points; //The list of points in 2D grid.

        //An indexed triangle element
        struct simplex{
            uint index;
            std::array<uint, 3> points;
        };
        std::vector<simplex> simplices; //The list of simplex triangles in 2D grid.

        //Make grid when initiated.
        Grid2D(double xMin, double xMax, double yMin, double yMax, uint xNum, uint yNum){
            double xLen = xMax-xMin;
            double yLen = yMax-yMin;
            
            for(uint i=0; i<xNum; i++){
                for(uint j=0; j<yNum; j++){
                    point p; 
                    p.x=((xLen*i)/(xNum-1)) + xMin;
                    p.y=((yLen*j)/(yNum-1)) + yMin;
                    p.index = j + xNum * i;
                    if(p.x == xMax || p.x == xMin || p.y == yMax || p.y == yMin) {p.boundary=true;}
                    points.push_back(p);
                }
            }
            calcSimplex(xMax,yMax,xNum,yNum);
        }
        void calcSimplex(double xMax, double yMax, uint xNum, uint yNum){
            //If point's x or y at max do nothing; else, make two triangles point's x->x+1 and y-y+1 region.
            uint index = 0;
            for(point p: points){
                if(p.x!=xMax && p.y!=yMax){
                    simplex s;
                    s = {index,p.index, p.index+yNum, p.index+1};
                    simplices.push_back(s);
                    index++;
                    s = {index,p.index+1+yNum, p.index+1, p.index+yNum};
                    simplices.push_back(s);
                    index++;
                }
            }
        }
        void printPoints(){
            printf("Points:[index-boundary-> x , y]\n");
            for (point p: points) {std::cout<<p.index<<"-"<<p.boundary<<"-> x:"<<p.x<<", y:"<<p.y<<std::endl;}
        }
        void printSimplex(){
            printf("Simplices:[index-> 3x Point Indicies]\n");
            for (simplex s: simplices){std::cout<<s.index<<"-> "<<s.points[0]<<" "<<s.points[1]<<" "<<s.points[2]<<std::endl;}    
        }
        
};