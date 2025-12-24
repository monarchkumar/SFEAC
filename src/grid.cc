#include <iostream>
#include <ostream>
#include <vector>
#include <optional>
#include <cmath> 
#include <map> 

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
        std::vector<std::vector<double>> inverse(){
            // double inv[size][size];

            auto a = [this](int i){return op[i-1][i-1];};
            auto b = [this](int i){return op[i-1][i];};
            auto c = [this](int i){return op[i][i-1];};
            // for (int i=1;i<6;i++){std::cout<<a(op,i)<<std::endl;}
            // for (int i=1;i<5;i++){std::cout<<b(op,i)<<std::endl;}
            // for (int i=1;i<5;i++){std::cout<<c(op,i)<<std::endl;}

            std::map<int,double> Theta; Theta[0] = 1; Theta[1] = a(1);
            for (int i=2; i<size; i++) {Theta[i] = a(i)*Theta[i-1] - b(i-1)* c(i-1)*Theta[i-2];}
            // for (auto val: Theta){std::cout<<val.first<<" "<<val.second<<std::endl;}
        
            double det = Theta[size-1];
            if (det < 10e-10){std::cout<<"Zero Determinant\n";}

            std::map<int,double> Phi; Phi[size+1] = 1; Phi[size] = a(size);
            for (int i=size-1; i>0; i--) {Phi[i] = a(i)*Phi[i+1] - b(i)* c(i)*Phi[i+2];}
            // for (auto val: Phi){std::cout<<val.first<<" "<<val.second<<std::endl;}
            
            std::vector<std::vector<double>> inv(size, std::vector<double>(size));
            for(int i=1;i<=size;i++){
                for (int j=1;j<=size;j++) {
                    inv.at(i-1).at(j-1) = pow(-1.0, i+j)/det;
                    // std::cout<<"-> "<<i<<" "<<j<<std::endl;

                    if(i<=j){
                        inv.at(i-1).at(j-1) = pow(-1.0, i+j)/det;

                        for(int k = i; k <= j-1; k++) {inv.at(i-1).at(j-1) *= b(k);}
                        inv.at(i-1).at(j-1) *= Theta[i-1];
                        inv.at(i-1).at(j-1) *= Phi[j+1];
                    } 
                    else {
                        inv.at(i-1).at(j-1) = pow(-1.0, i+j)/det; 
                        for(int k = j; k <= i-1; k++) {inv.at(i-1).at(j-1) *= c(k);}
                        inv.at(i-1).at(j-1) *= Theta[j-1];
                        inv.at(i-1).at(j-1) *= Phi[i+1];
                    } 
                }
            }
            return inv;
            // for(int i=1;i<=size;i++){
            //     for (int j=1;j<=size;j++) {
            //         if(inv.at(i-1).at(j-1)<=0) printf("%.3e ", inv.at(i-1).at(j-1));
            //         else printf(" %.3e ", inv.at(i-1).at(j-1));
            //     }
            //     std::cout<<std::endl;
            // }

        }   
        
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