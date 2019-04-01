#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>

void readCsv(std::string filename, std::vector<double>* a){
    a->resize(0);
   std::ifstream infile(filename);
   if(! infile)
       throw std::runtime_error("Could not open file " + filename);
   std::string line;
   int row = 0;
   while (std::getline(infile, line))
   {
     if (line.back() == ';') {
       a->push_back(std::stod(line));
       //std::cout << "number: " << (*a)[row] << std::endl;
       row++;
     }
   }
}

void writeCsv(std::string filename, std::vector<double> a){
    std::ofstream myfile;
    myfile.open (filename);
    for(int i = 0; i < a.size(); i++){
        myfile << a[i] << ";" << std::endl;
    }
    myfile << std::endl;
    myfile.close();
}


void writeCsv(std::string filename, std::vector<std::vector<double>> a){
    std::ofstream myfile;
    myfile.open (filename);
    for(int i = 0; i < a.size(); i++){
        for(int j = 0; j < a[i].size(); j++){
            myfile << a[i][j] << ";";
        }
        myfile <<  std::endl;
    }
    myfile << std::endl;
    myfile.close();
}
