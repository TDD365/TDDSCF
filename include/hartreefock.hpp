#ifndef _HF_HPP
#define _HF_HPP

// Hartree-Fock
#include "mol.hpp"
#include "eint.hpp"
#include "maths.hpp"

// Single Slater Object
class HartreeFock {
    private:
    public:
        // Molecule Pointer
        Eint* IntMole;
        Molecule* Mol_pointer;

        // Properties
        int NBF;
        MatrixRowMajor SHalf;         // S^{-1/2}
        MatrixRowMajor coeff;         // MO Coefficient
        MatrixRowMajor Den;           // Density Matrix
        MatrixRowMajor Fock;          // Fock Matrix
        MatrixRowMajor oDen;          // Orthogonal Density Matrix
        MatrixRowMajor oFock;         // Orthogonal Fock Matrix

        // Controls
        double conv = 1e-6;           // Convergence
        int maxiter = 512;            // Iterations

        // Constructor
        HartreeFock(Eint* Buffer){
            // Assign Pointer
            IntMole = Buffer;
            Mol_pointer = Buffer->Mol_pointer;
            // Matrix size
            NBF = IntMole->NBF;
        }

        // Destructor
        ~HartreeFock(){
            IntMole = nullptr;
            Mol_pointer = nullptr;
            std::cout << "Cleaning HF Object..." << endl;
        }

        // Functions
        // Before SCF
        void get_HCore();
        void get_Guess();

        // During SCF
        void get_Fock();
        void get_SHalf();
        void get_oFock();
        void get_oDen();

        // SCF
        void do_SCF();

        // After SCF
        void get_charge();
        void print_result();

};

#endif // _HF_HPP