#include <iostream>
#include <ostream>
#include <vector>
#include <optional>
#include <cmath> 

class Operator{
    public:
        std::vector<std::vector<double>> op; //Operator 
        bool defined = false;
        int size; double h; 

        //Dirichlet both-boundary
        Operator(int arraySize){
            size = arraySize; h=1.0/(arraySize-1);
        }
        void delPhiiDelPhij(){
            op.resize(size);
            for (int i = 0; i<size; i++){
                op[i].resize(size);
                op[i][i]=2;
                if (i==0){op[i][i+1]=-1;}
                if (i==op.size()-1){op[i][i-1]=-1;}
                if (i>0 && i<op.size()-1) {
                    op[i][i-1]=-1;
                    op[i][i+1]=-1;
                }
                for(double& i : op[i]) { i *= (1/h);}
                
            }
            defined = true;
        }
        void Phii(double f){
            op.resize(1);
            op[0].resize(size);
            for(double& i : op[0]) {i=f;}     
            op[0].at(0)/=2;           
            op[0].at(size-1)/=2;           
            defined = true;
        }

        void dirichlet(int location, double value){
            if(!defined){std::cout<<"Operator not defined"<<std::endl; return;}
            if(op.size()==1){op[0][location]=value;}
            else if (op.size()>1){
                for (double& i: op[location]){i=0;}
                op.at(location).at(location)=1.0;
            }
        }
        //C.M. da Fonseca / Journal of Computational and Applied Mathematics 200 (2007) 283 – 286
        // void inverse(){
        //     double Theta = 1;
        //     double inv[op.size()][op.size()];
        //     for(int i=0;i<op.size();i++){
        //         for (int j = 0; j<op[i].size(); j++) {
        //             Theta = inv[i][j] * Theta - inv[i-1][i] - inv[i-1][i]
        //             if (i<=j){
        //                 inv[i][j] = std::pow(-1.0,i+j);
        //                 for(int k=i;i<j-1;k++){inv[i][j] *= inv[k][k+1];} 
        //                 inv[i][j] *= inv[i][i]
        //             }
        //             else {}
        //         }
        //     }
        // }
        
        void print(){
            std::cout<<"operator: "<<std::endl;
            for (int i = 0; i<op.size(); i++){
                for (auto val: op.at(i)){
                    if (val<0) {std::cout<<" "<<val;}
                    else  {std::cout<<"  "<<val;}
                }
                std::cout<<std::endl;
            }
            std::cout<<"size: "<<size<<std::endl;
            std::cout<<"h: "<<h<<std::endl;
        }
};


class Functions: protected Operator{
    public:
        std::vector<double> positions;
        std::vector<std::optional<double>> values;
        double h; 

        enum dimensional {D1,D2,D3};

        //Dirichlet both-boundary
        Functions(dimensional dimension, int length, double xMin, double xMax, double xMinVal, double xMaxVal):Operator(length){
            if (dimension == D1){
                h=(xMax-xMin)/(length-1);
                for (int i=0;i<length;i++){positions.push_back(i*h);};
                values.resize(length);
                values[0]=xMinVal;
                values[length-1]=xMaxVal;
            }
        }

        void printValues(){
            std::cout<<"Printing Function Parameters"<<std::endl;
            std::cout<<std::endl;

            std::cout<<"h:"<<h<<std::endl;
            std::cout<<std::endl;

            std::cout<<"position:"<<std::endl;
            for (auto &val: positions){std::cout<<val<<" ";}
            std::cout<<std::endl;

            std::cout<<std::endl<<"value:"<<std::endl;
            for (auto &val: values){
                if(val.has_value()){std::cout<<val.value()<<" ";}
                else {std::cout<<"NoValue ";}
            }
            std::cout<<std::endl;
        }
};