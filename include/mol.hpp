#ifndef _MOLE_HPP
#define _MOLE_HPP

#include <vector>
#include <string>
#include <iostream>
#include "maths.hpp"

using namespace std;

// Atoms
class Atoms
{
    private:
    public:
        // Variables
        int id;            // Atomic Number
        double mass;       // Atomic mass
        double x, y, z;    // Position
        string nid;        // Name

        // Constructor
        Atoms() {
            id=0;
            mass=0.;
            x=0.;
            y=0.;
            z=0.;
            nid="X";
        }
        Atoms(string s, double r[]);
        Atoms(int s, double r[]);

        // Destructor
        ~Atoms(){};
};

// Molecule
class Molecule{

    private:
    public:
        // Property - Classical
        int natoms;   // Number of Atoms
        double enuc;  // Nuclear Repulsion Energy
        double ekin;  // Nuclear Kinetic Energy
        double eele;  // Electronic Energy
        double etol;  // Total Energy

        // Property - Quamtum Chemistry
        int NBF = 0;                  // Total size of Basis
        int Nele = 0;                 // Number of electrons
        MatrixRowMajor Ovlp;          // Overlap
        MatrixRowMajor HVel;          // Core Hamiltonian - Velocity Part
        MatrixRowMajor HNuc;          // Core Hamiltonian - Nuclear Part
        MatrixRowMajor HCore;         // Core Hamiltonian
        TensorRowMajor ERI;           // ERI

        // Object Variables
        std::vector<Atoms> Structure;

        // Constructor
        Molecule(){
            natoms=-1;
            enuc=-99999999.;
        };

        // Destructor
        ~Molecule(){
            cout << "Cleaning Molecule Object..." << std::endl;}

        // Functions Classical
        void read_geom(const std::string& filename);
        void calc_enuc();
        void print_xyz();

};


#endif // _MOLE_HPP