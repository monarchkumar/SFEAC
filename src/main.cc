#include <iostream>
#include "version.cc"
#include "grid.cc"

using namespace std;
#define NUMBER_OF_DASH 64
#define DEBUG true

string version = "0.1";
int main(){
    cout << string(NUMBER_OF_DASH, '-') << endl;
    cout << "SFEAC Version: " << version << endl;

    if(DEBUG){printCppVersion(NUMBER_OF_DASH);}

    //Setup
    // Functions testFunction(Functions::D1,5,0,1,0,0);
    // testFunction.gradU_gradV();
    // testFunction.printValues();

    Operator A(5);
    A.delPhiiDelPhij();
    A.dirichlet(0,0);
    A.dirichlet(4,1);
    // A.print();
    // cout << string(NUMBER_OF_DASH, '-') << endl;


    Operator f(5);
    f.Phii(1.0/4);
    f.dirichlet(0,0);
    f.dirichlet(4,1);
    // f.print();

       
    auto invA= A.inverse();

    // for(int i=1;i<=A.size;i++){
    //     for (int j=1;j<=A.size;j++) {
    //         if(invA[i-1][j-1]<=0) printf("%.3e ", invA[i-1][j-1]);
    //         else printf(" %.3e ", invA[i-1][j-1]);
    //     }
    //     std::cout<<std::endl;
    // }




    cout << string(NUMBER_OF_DASH, '-') << endl;
    cout << "Exiting" << endl;
    cout << string(NUMBER_OF_DASH, '-') << endl;
    return 0;
}
/*

   

    if(DEBUG){printCppVersion(NUMBER_OF_DASH);}

    // I/O initiated;
    string inputFilePath;
    string outputFilePath;
    ofstream outputFile;
    
    // Input files are YAML.
    // cout << "File directory: ";
    // cin >>inputFilePath;
    // TINY_YAML::Yaml inputFile(inputFilePath + "input.yaml");
    // TINY_YAML::Yaml inputFile("/home/mkumar/SFEAC/input/""input.yaml");
    cout << "[HARDCODED] File directory: /home/mkumar/SFEAC/input/" << endl;
    cout << string(NUMBER_OF_DASH, '-') << endl;
    
    cout << "Input file found." << endl;
    // cout << "Title: " << inputFile["metadata"]["label"].getData<string>() << endl;
    // outputFilePath = inputFile["metadata"]["output"].getData<string>();
    cout << "Output: " << outputFilePath << endl;
    outputFile.open(outputFilePath+"/output.txt");
    outputFile << "Undefined Version: " << version << endl;
    cout << string(NUMBER_OF_DASH, '-') << endl;

    // Define solver
    cout<<"here"<<endl;
    // if(inputFile["solver"]["label"].getData<string>() == "ES1D"){
        // G1D solver;
        // solver.m = 1;
        // solver.xMin = 0;
        // solver.xMax = 10;
        // solver.xGridPointsNum = 5;
        // solver.noPartNum = 10;
        // solver.initDist = solver.UniformDistribution;
        // solver.initG1D();
        // solver.solveG1D();
        //---------------------------------------------------------------------------
        // string mass = inputFile["particle"]["mass"].getData<string>();
        // std::size_t* pos = 0 ;
        // cout << std::stod(mass,pos)<< endl;
        //---------------------------------------------------------------------------



    // }
    // else {
    //     cout << "Error: No solver specified." << endl;
    // }

    outputFile.close();
    cout << string(NUMBER_OF_DASH, '-') << endl;
    cout << "Exiting" << endl << endl;
    return 0;
}*/