#include "controls.hpp"

// Read Controls
void readControls(const std::string& filename) {

    string linecontent;
    int reading = 1;
    string x;
    string buf_int;

    std::cout << "Input file detected: " << filename << std::endl;
    std::cout << "======= Reading Control =======" << std::endl;
    std::ifstream infile(filename);

    // General Control
    if (infile.is_open()) {
        while (getline(infile, linecontent)) {
            std::stringstream ss(linecontent);
            // Controls
            ss >> x >> buf_int;
            if (reading == 1) {
                // General
                if (x=="PrintLevel")
                    prtlevel = stoi(buf_int);
                if (x=="Guess")
                    Guess_Method = buf_int;
                if (x=="Maxiter")
                    scfMaxiter = stoi(buf_int);
                if (x=="Scfcon")
                    scfcon = stoi(buf_int);
                if (x=="Dencon")
                    dencon = stoi(buf_int);
                }
            if (x=="[Controls]")
                std::cout << "Control read: " << std::endl;
            if (x=="[SCF]")
                std::cout << "SCF Control read: " << std::endl;
        }
    }
    std::cout << "Print Level: " << prtlevel << std::endl;
    infile.close();

}