// All subroutines to build Fock Matrix
#include "hartreefock.hpp"
#include "controls.hpp"
#include "maths.hpp"
using namespace std;

// Pre-SCF
// 1. Get HCore
void HartreeFock::get_HCore() {
    Mol_pointer->HCore = Mol_pointer->HNuc + Mol_pointer->HVel;
    if (prtlevel>1) {
        std::cout << "=== HCore ===" << std::endl;
        std::cout << Mol_pointer->HCore <<endl;
        std::cout << "=========" << std::endl;
    }
    get_SHalf();  
}
// Get S{-1/2}
void HartreeFock::get_SHalf() {
    SHalf = Mol_pointer->Ovlp.sqrt().inverse();
    if (prtlevel>1) {
        std::cout << "=== S^-1/2 ===" << std::endl;
        std::cout << this->SHalf <<endl;
        std::cout << "=========" << std::endl;
    }    
}

// 2. Form Guess
void HartreeFock::get_Den() {
    if (Mol_pointer->Nele % 2)
        std:cerr << "Non-Singlet NYI!";
    // Get Den
    auto Coeff_occ = this->Coeff.leftCols(Mol_pointer->Nele/2);
    this->Den = 2.0 * Coeff_occ * Coeff_occ.transpose();
    if (prtlevel>1) {
        std::cout << "=== Current Density ===" << std::endl;
        std::cout << this->Den << std::endl;
        std::cout << "=========" << std::endl;
    }    
}
void HartreeFock::get_Guess_Den() {
    this->Coeff = MatrixRowMajor::Zero(NBF,NBF);
    this->Den = MatrixRowMajor::Zero(NBF,NBF);
    std::cout << "=== Getting Initial Guess ===" << std::endl;
    std::cout << "Method is " << Guess_Method << std::endl;
    if (Mol_pointer->Nele % 2)
        std:cerr << "Non-Singlet NYI!";
    if (Guess_Method=="Dumb") {
        for(auto i=0; i<Mol_pointer->Nele/2; i++)
            this->Coeff(i,i) = 1.0;
    } else if (Guess_Method=="Core") {
        this->Fock = Mol_pointer->HCore;
        this->solve_FCeSC();
    }
    this->get_Den();
    this->get_Energy();
}

// During-SCF
// 3. Get Fock
void HartreeFock::get_Fock() {
	// Do contraction
    // This should be simplified, but I'm just lazy... come bite me...
    this->Fock = Mol_pointer->HCore;
	for(auto i=0; i<NBF; i++) {
		for(auto j=0; j<=i; j++) {
			for(auto k=0; k<NBF; k++) {
				for(auto l=0; l<NBF; l++) {
					this->Fock(i,j) += 1.0 * this->Den(k,l) * Mol_pointer->ERI(i,j,k,l);
					this->Fock(i,j) -= 0.5 * this->Den(k,l) * Mol_pointer->ERI(i,k,l,j);
				}
			}
			if (j != i)
				this->Fock(j,i) = this->Fock(i,j);
		}
    }
    if (prtlevel>0) {
        std::cout << "=== Fock ===" << std::endl;
        std::cout << this->Fock << std::endl;
        std::cout << "=========" << std::endl;
    }
}

// 4. Diag F
void HartreeFock::solve_FCeSC() {
	// Solve Eigenvalue Problem
    Eigen::GeneralizedSelfAdjointEigenSolver<MatrixRowMajor> HFSolver(this->Fock, 
                                                                      Mol_pointer->Ovlp);
    this->EOrbs = HFSolver.eigenvalues();
    this->Coeff = HFSolver.eigenvectors();
    if (prtlevel>0) {
        std::cout << "=== after FCeSC ===" << std::endl;
        std::cout << this->Coeff <<endl;
        std::cout << "=========" << std::endl;
    }
}

// 5. Calculate Energy
void HartreeFock::get_Energy() {
	Mol_pointer->eele = 0.;
    // RHF
    // Can we add something like SCF.Type='RHF'/'UHF'?
	for(auto i=0; i<NBF; i++)
        for(auto j=0; j<NBF; j++)
            Mol_pointer->eele += 0.5 * this->Den(i,j) 
                                     * (Mol_pointer->HCore(j,i) + this->Fock(j,i));
    Mol_pointer->etol = Mol_pointer->enuc + Mol_pointer->eele;
}
