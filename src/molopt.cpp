#include "mol.hpp"
#include "molopt.hpp"
#include "constants.hpp"

#include <iostream>
#include <cstring>
#include <fstream>
#include <iomanip>

void Molecule::read_geom(const std::string& filename) {

    string linecontent;
    int reading;
    string x;
    double r[3];

    std::cout << "======= Reading Molecule =======" << std::endl;
    std::ifstream infile(filename);

    reading = 0;
    natoms = 0;
    if (infile.is_open()) {
        while (getline(infile, linecontent) && reading>=0) {
            std::stringstream ss(linecontent);
            // Controls
            //ss >> x;
            ss >> x >> r[0] >> r[1] >> r[2];
            if (x=="[MolEnd]")
                reading = -1;
            if (reading == 1) {
                //std::cout << x << " " << r[0] << " " << r[1] << " " << r[2] << std::endl;
                Structure.push_back(Atoms(x, r));
                Nele+=Atoms(x, r).id;
                natoms+=1;
                }
            if (x=="[Mol]")
                reading = 1;
        }
    }
    //std::cout << " ====== End of Molecule Structure ====== " << std::endl;
    std::cout << "Mol structure read. " << std::endl;
    infile.close();
}
