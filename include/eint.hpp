#ifndef _INT_HPP
#define _INT_HPP

using namespace std;
#include "mol.hpp"
#include "maths.hpp"
#include <iostream>
#include <libint2.hpp>

class Eint {
    private:
    public:
        // Molecule Pointer
        Molecule* Mol_pointer;

        // Properties
        int NBF;
        string Basis_method, Basis_name;

        // Using libint's default
        libint2::Atom atom_temp;
        vector<libint2::Atom> int_atoms;
        libint2::BasisSet BSs;

        // If reading from external
        string infile = "Infile.h5";

        // Constructor
        Eint(Molecule* QMole){
            // Assign Pointer
            Mol_pointer = QMole;
            // Build Atom Object
            for (auto i=0; i<QMole->Structure.size(); i++) {
                atom_temp.atomic_number = QMole->Structure[i].id;
                atom_temp.x = QMole->Structure[i].x;
                atom_temp.y = QMole->Structure[i].y;
                atom_temp.z = QMole->Structure[i].z;
                int_atoms.push_back(atom_temp);
            }
        }

        // Destructor
        ~Eint(){
            Mol_pointer = nullptr;
            std::cout << "Cleaning Integral Object..." << endl;}

        // Functions
        void read_int(const std::string& filename);
        void build_shell();
        void compute_ovlp();
        void compute_h1e_nuc();
        void compute_h1e_kin();
        void compute_eri();
        void compute();

};

#endif //_INT_HPP
