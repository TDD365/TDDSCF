#include "mol.hpp"
#include "constants.hpp"

#include <iostream>
#include <cstring>
#include <fstream>
#include <iomanip>

// Constructors for Atoms
Atoms::Atoms(string s, double r[]) {

    nid = s;
    if (s == "H") {
        id = 1;
        mass = 1.00782504;
    } else if (s == "C") {
        id = 6;
        mass = 12.;
    } else if (s == "N") {
        id = 7;
        mass = 14.00307401;
    } else if (s == "O") {
        id = 8;
        mass = 15.994914;
    } else if (s == "F") {
        id = 9;
        mass = 17.002095238;
    } else {
        std::cerr << "Element " << s <<" not implemented!";
    }

    x = r[0]*au2bohr;
    y = r[1]*au2bohr;
    z = r[2]*au2bohr;

}
Atoms::Atoms(int s, double r[]) {

    id = s;
    if (s == 1) {
        nid = "H";
        mass = 1.00782504;
    } else if (s == 6) {
        nid = "C";
        mass = 12.;
    } else if (s == 7) {
        nid = "N";
        mass = 14.00307401;
    } else if (s == 8) {
        nid = "O";
        mass = 15.994914;
    } else if (s == 9) {
        nid = "F";
        mass = 17.002095238;
    } else {
        std::cerr << "Element ID " << s <<" not implemented!";
    }

    x = r[0]*au2bohr;
    y = r[1]*au2bohr;
    z = r[2]*au2bohr;

}

// Functions in Molecule
// Calculate Nuclear Repulsion Energy
void Molecule::calc_enuc() {
    enuc = 0.;
    for (auto i=0; i<natoms; i++) {
        for (auto j=i+1; j<natoms; j++) {
            auto rij = sqrt(
                (Structure[i].x - Structure[j].x)*(Structure[i].x - Structure[j].x)
               +(Structure[i].y - Structure[j].y)*(Structure[i].y - Structure[j].y)
               +(Structure[i].z - Structure[j].z)*(Structure[i].z - Structure[j].z));
            enuc += Structure[i].id * Structure[j].id / rij;
        }
    }
    cout << "Nuclear Repulsion Energy (A.U.) :" 
    << setw(9) << setprecision(6) << enuc
    << std::endl;
}

// Print Structure
void Molecule::print_xyz() {
    for (auto i=0; i<natoms; i++) {
        cout << setw(4) << Structure[i].nid << " " 
        << setw(9) << setprecision(6) << Structure[i].x << " " 
        << setw(9) << setprecision(6) << Structure[i].y << " " 
        << setw(9) << setprecision(6) << Structure[i].z << " " 
        << std::endl;
    }
}
