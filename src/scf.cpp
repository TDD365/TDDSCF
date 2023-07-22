// SCF Process
#include "hartreefock.hpp"
#include "controls.hpp"
#include "maths.hpp"
using namespace std;

void HartreeFock::SCFStep() {

    // We only have Density at last step! Get Fock first!
    this->get_Fock();

    // Solve Eigenvalue Problem
    Eigen::GeneralizedSelfAdjointEigenSolver<MatrixRowMajor> HFSolver(this->Fock, 
                                                                      Mol_pointer->Ovlp);
    this->EOrbs = HFSolver.eigenvalues();
    this->Coeff = HFSolver.eigenvectors();

    // Get New Density
    this->get_Den();

    // Get New Energy
    this->get_Energy();

}


void HartreeFock::do_SCF() {

    int iter = 0;
    double energy_diff;
    double den_diff;
    double conv = pow(10,-scfcon);
    double denconv = pow(10,-dencon);
    

    // Print SCF Control Info
    std::cout << "=== SCF ===" << std::endl;
    std::cout << "=== SCF Setting ===" << std::endl;

    // TODO: Add Density-based Gradients; add DIIS;
    std::cout << "Method is " << "Conventional" << std::endl;
    std::cout << "SCF Iter " << scfMaxiter << std::endl;
    std::cout << "SCF Convergence " << conv << std::endl;
    std::cout << "Density Convergence " << denconv << std::endl;

    // Get Guess -- should I put get_HCore outside?
    this->get_HCore();
    this->get_Guess_Den();

    // Start SCF
    std::cout << "=== Do SCF ===" << std::endl;
    do {
        // Keep in mind, energy For guess step makes no sense
        auto energy_last = Mol_pointer->eele;
        MatrixRowMajor den_last = this->Den;

        // Update MO
        SCFStep();

        // Diff
        energy_diff = Mol_pointer->eele - energy_last;
        den_diff = (this->Den - den_last).norm();

        // Print
        // Suppress a print level
        prtlevel--;
        if (iter==0)
            printf("%6s\t %20s\t %20s\t %20s\n", 
                   "Iter", "E(elec)", "Delta(E)", "Delta(D)");
        if (iter>0) {
            printf("%4d\t %20.12f\t %20.12f\t %20.12f\n", 
                   iter, Mol_pointer->eele, energy_diff, den_diff);
        } else {
            printf("%4d\t %20.12f\t %20s\t %20s\n", 
                   iter, Mol_pointer->eele, "N/A", "N/A");
        }

        iter++;
    } while (iter <= maxiter && 
             (fabs(energy_diff) >= conv ||
             fabs(den_diff) >= denconv));

    // Print Result
    if (iter <= maxiter) {
        std::cout << "=== SCF Done in " << iter << " Cycles ===" << std::endl;
    } else {
        std::cout << "=== SCF Process Still Not Converged ===" << std::endl;
    }
    std::cout << "=== Final Energy ===" << std::endl;
    printf("E-Elec:\t %20.12f\nE-Nuc:\t %20.12f\nE-Tot:\t %20.12f\n",
            Mol_pointer->eele, Mol_pointer->enuc, Mol_pointer->etol);
}