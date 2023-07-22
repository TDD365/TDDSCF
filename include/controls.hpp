#ifndef _CTRL_HPP
#define _CTRL_HPP

#include <iostream>
#include <cstring>
#include <fstream>
#include <iomanip>
using namespace std;

// Options
inline int prtlevel = 0;
inline int scfMaxiter = 64;
inline int scfcon = 5;
inline int dencon = 5;
inline std::string Guess_Method = "Dumb";

// Read
void readControls(const std::string& filename);

#endif // _CTRL_HPP