#include "controls.hpp"

// Read Controls
void readControls(const std::string& filename) {

    string linecontent;
    int reading;
    string x;
    int buf_int;

    std::cout << "Input file detected: " << filename << std::endl;
    std::cout << "======= Reading Control =======" << std::endl;
    std::ifstream infile(filename);

    reading = 0;
    if (infile.is_open()) {
        while (getline(infile, linecontent) && reading>=0) {
            std::stringstream ss(linecontent);
            // Controls
            //ss >> x;
            ss >> x >> buf_int;
            if (x=="[End]")
                reading = -1;
            if (reading == 1) {
                if (x=="PrintLevel")
                    prtlevel = buf_int;
                }
            if (x=="[Controls]")
                reading = 1;
        }
    }
    std::cout << "Control read: " << std::endl;
    std::cout << "Print Level: " << prtlevel << std::endl;
    infile.close();

}