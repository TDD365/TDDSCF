#ifndef _HF_HPP
#define _HF_HPP

// Hartree-Fock
#include "mol.hpp"
#include "eint.hpp"
#include "maths.hpp"

// Single Slater Object
// Should set it as template afterwards, for RHF, UHF, X2C, ...
class HartreeFock {
    private:
    public:
        // Molecule Pointer
        Eint* IntMole;
        Molecule* Mol_pointer;

        // Properties
        int NBF;
        MatrixRowMajor  SHalf;        // S^{-1/2}
        MatrixRowMajor  Coeff;        // MO Coefficient
        MatrixRowMajor  Den;          // Density Matrix
        MatrixRowMajor  Fock;         // Fock Matrix
        MatrixRowMajor  oDen;         // Orthogonal Density Matrix
        MatrixRowMajor  oFock;        // Orthogonal Fock Matrix
        Eigen::VectorXd EOrbs;        // Eigevalues

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
        void get_Guess_Den();
        void get_SHalf();

        // During SCF
        void get_Den();
        void get_Fock();
        void solve_FCeSC();
        void get_Energy();

        // SCF
        void SCFStep();
        void do_SCF();

        // After SCF
        void get_charge();
        void print_result();

};

#endif // _HF_HPP