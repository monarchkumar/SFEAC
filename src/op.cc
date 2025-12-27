#include <fstream>
#include <iostream>
#include <vector>
#include <cmath> 
#include <map> 
#include <format> 

class Operator{
    public:
        std::vector<std::vector<double>> op; //Operator 
        // std::vector<std::vector<double>> inv; //Inverse
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
            op[0][0]/=2;
            op[size-1][size-1]/=2;
            defined = true;
        }
        void Phii(double f){
            op.resize(size);
            for(auto& row : op) {row={f};}     
              
            op.at(0)[0]/=2;           
            op.at(size-1)[0]/=2;  
            
            defined = true;
        }

        void dirichlet(int location, double value){
            if(!defined){std::cout<<"Operator not defined"<<std::endl; return;}
            if(op.front().size()==1){op[location].front()=value;}
            else {
                for (double& i: op[location]){i=0;}
                op.at(location).at(location)=1.0;
            }
        }

        void neumann(int location, double value){
            if(!defined){std::cout<<"Operator not defined"<<std::endl; return;}
            if(op.front().size()==1){op[location].front()+=value;}
        }


        //C.M. da Fonseca / Journal of Computational and Applied Mathematics 200 (2007) 283 – 286
        Operator inverse(){
            Operator resOp(this->size);

            auto a = [this](int i){return op[i-1][i-1];};
            auto b = [this](int i){return op[i-1][i];};
            auto c = [this](int i){return op[i][i-1];};

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
            resOp.op = inv;
            return resOp;
           

        }   

        Operator multi(Operator withOp){

            Operator resOp(this->size);
            if(this->op.front().size() != withOp.op.size()) {std::cout<<"Not compatable."<<std::endl;return resOp;} 


            std::vector<std::vector<double>> res(op.size(), std::vector<double>(withOp.op.front().size()));

            for (int i =0; i<this->op.size();i++){
                for(int j=0; j<withOp.op.front().size(); j++){
                    res[i][j] = 0;
                    for (int k = 0; k < op.front().size(); k++) {
                        res[i][j] += op[i][k] * withOp.op[k][j];
                    }
                }
            }
            resOp.op = res;
            return resOp;
        }
        
        void print(){
            std::cout<<"operator: "<<std::endl;
            for (int i = 0; i<op.size(); i++){
                for (double val: op.at(i)){
                    std::string s =  std::format("{:.3e} ",val);
                    if(s[0]!='-') s = ' ' + s;
                    std::cout<<s;
                }
                std::cout<<std::endl;
            }

        }
        void opToCSV(std::ofstream& outputFile){
            for (int i = 0; i<op.size(); i++){
                std::string row = "";
                for (double val: op.at(i)){
                    std::string cell =  std::format("{:.3e},",val);
                    if(cell[0]!='-') cell = ' ' + cell;
                    row+=cell;
                }
                row.pop_back();
                outputFile<<row<<std::endl;
            }

        } 
};