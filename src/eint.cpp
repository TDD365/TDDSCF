#include "mol.hpp"
#include "eint.hpp"
#include "H5Cpp.h"

// Read Basis & Int options
void Eint::read_int(const std::string& filename) {

    string linecontent;
    int reading;
    string x, opt_1, opt_2;

    std::cout << "======= Reading Basis Set =======" << std::endl;
    std::ifstream infile(filename);

    if (infile.is_open()) {
        while (getline(infile, linecontent) && reading>=0) {
            std::stringstream ss(linecontent);
            // Controls
            ss >> x >> opt_1 >> opt_2;
            //std::cout << x << " " << opt_1 << " " << opt_2 << std::endl;
            if (x=="Basis") {
                if (opt_1=="inLib") {
                    Basis_method = opt_1;
                    Basis_name = opt_2;
                } else {
                    std::cerr << "Custom Basis NYI!";
                }
            }
        }
    }
    std::cout << "Basis Set info read." << std::endl;
    infile.close();
}

void Eint::build_shell() {

    // Set up shell
    if (Basis_method=="inLib") {
        BSs = libint2::BasisSet(Basis_name, int_atoms);
        NBF = BSs.nbf();
    } else {
        std::cerr << "Custom Basis NYI!";
    }

    // Show shells
    std::cout << "=== Shell Info ===" << std::endl;
    std::cout << "Basis: " << Basis_method << " " << Basis_name << std::endl;
    std::cout << "Number of Basis Functions: " << BSs.nbf() << std::endl;
    std::cout << "Number of Primitive Gaussian Functions: " << BSs.max_nprim() << std::endl;
    std::cout << "Maximum Angular Momentum: " << BSs.max_l() << std::endl;

    // Print Basis Info in detail
    for(auto& Bs: BSs) {
    }

}

// Overlap Matrix
void Eint::compute_ovlp() {

    MatrixRowMajor Result(NBF, NBF);
    libint2::Engine h1e_engine(libint2::Operator::overlap,
                               BSs.max_nprim(),
                               BSs.max_l()
                               );
    const auto& buf = h1e_engine.results();
    // Do integral
    for(auto s1=0; s1!=BSs.size(); ++s1) {
        for(auto s2=0; s2<=s1; ++s2) {
            // Do integral
            h1e_engine.compute(BSs[s1], BSs[s2]);
            // Map it
            // Libint give a row-major order result
            // Eigen Map: Array -> (rows, cols)
            Eigen::Map<const MatrixRowMajor> buffermat(buf[0], 
                                                       BSs[s1].size(), 
                                                       BSs[s2].size());
            // Map the sub-matrix to full matrix
            // Eigen block: (i,j,p,q) -> Matrix（p,q) at position (i,j)
            // Here BSs.shell2bf() is starting place for all shells.
            Result.block(BSs.shell2bf()[s1], 
                         BSs.shell2bf()[s2], 
                         BSs[s1].size(), 
                         BSs[s2].size()) = buffermat;
            if (s1 != s2) {
                Result.block(BSs.shell2bf()[s2], 
                             BSs.shell2bf()[s1], 
                             BSs[s2].size(), 
                             BSs[s1].size()) = buffermat.transpose();
            }
        }
    }
    Mol_pointer->Ovlp = Result.replicate(1,1);
    std::cout << "=== S ===" << std::endl;
    std::cout << Mol_pointer->Ovlp <<endl;
    std::cout << "=========" << std::endl;

}

// H1E - Kinetic
void Eint::compute_h1e_kin() {

    MatrixRowMajor Result(NBF, NBF);
    libint2::Engine h1e_engine(libint2::Operator::kinetic,
                               BSs.max_nprim(),
                               BSs.max_l()
                               );
    const auto& buf = h1e_engine.results();
    // Do integral
    for(auto s1=0; s1!=BSs.size(); ++s1) {
        for(auto s2=0; s2<=s1; ++s2) {
            // Do integral
            h1e_engine.compute(BSs[s1], BSs[s2]);
            // Map it
            // Libint give a row-major order result
            // Eigen Map: Array -> (rows, cols)
            Eigen::Map<const MatrixRowMajor> buffermat(buf[0], 
                                                       BSs[s1].size(), 
                                                       BSs[s2].size());
            // Map the sub-matrix to full matrix
            // Eigen block: (i,j,p,q) -> Matrix（p,q) at position (i,j)
            // Here BSs.shell2bf() is starting place for all shells.
            Result.block(BSs.shell2bf()[s1], 
                         BSs.shell2bf()[s2], 
                         BSs[s1].size(), 
                         BSs[s2].size()) = buffermat;
            if (s1 != s2) {
                Result.block(BSs.shell2bf()[s2], 
                             BSs.shell2bf()[s1], 
                             BSs[s2].size(), 
                             BSs[s1].size()) = buffermat.transpose();
            }
        }
    }
    Mol_pointer->HVel = Result.replicate(1,1);
    std::cout << "=== T ===" << std::endl;
    std::cout << Mol_pointer->HVel <<endl;
    std::cout << "=========" << std::endl;

}

// H1E - Nuclear
void Eint::compute_h1e_nuc() {

    MatrixRowMajor Result(NBF, NBF);
    libint2::Engine h1e_engine(libint2::Operator::nuclear,
                               BSs.max_nprim(),
                               BSs.max_l()
                               );
    const auto& buf = h1e_engine.results();

    // Prepare atomic info
    std::vector<std::pair<double,std::array<double,3>>> pair_info;
    for(const auto& atom : int_atoms) {
        pair_info.push_back({static_cast<double>(atom.atomic_number), 
                            {{atom.x, atom.y, atom.z}}});
    }
    h1e_engine.set_params(pair_info);

    // Do integral
    for(auto s1=0; s1!=BSs.size(); ++s1) {
        for(auto s2=0; s2<=s1; ++s2) {
            // Do integral
            h1e_engine.compute(BSs[s1], BSs[s2]);
            // Map it
            // Libint give a row-major order result
            // Eigen Map: Array -> (rows, cols)
            Eigen::Map<const MatrixRowMajor> buffermat(buf[0], 
                                                       BSs[s1].size(), 
                                                       BSs[s2].size());
            // Map the sub-matrix to full matrix
            // Eigen block: (i,j,p,q) -> Matrix（p,q) at position (i,j)
            // Here BSs.shell2bf() is starting place for all shells.
            Result.block(BSs.shell2bf()[s1], 
                         BSs.shell2bf()[s2], 
                         BSs[s1].size(), 
                         BSs[s2].size()) = buffermat;
            if (s1 != s2) {
                Result.block(BSs.shell2bf()[s2], 
                             BSs.shell2bf()[s1], 
                             BSs[s2].size(), 
                             BSs[s1].size()) = buffermat.transpose();
            }
        }
    }
    Mol_pointer->HNuc = Result.replicate(1,1);
    std::cout << "=== V ===" << std::endl;
    std::cout << Mol_pointer->HNuc <<endl;
    std::cout << "=========" << std::endl;

}

// ERI
// In human words, Electronic repulsion integral
void Eint::compute_eri() {

    MatrixRowMajor Result(NBF*NBF, NBF*NBF);
    libint2::Engine h2e_engine(libint2::Operator::coulomb,
                               BSs.max_nprim(),
                               BSs.max_l()
                               );
    const auto& buf = h2e_engine.results();

    // Do integral
    // Here we will calculate (ij|kl)
    // (12|34) = (21|34) = (12|43) = (21|43) = (34|12) = (43|12) = (34|21) = (43|21)
    // Yes I know there are 4-fold symmetries, we will never encounter that in this code!

    // Go over all shells...
    for(auto i=0; i<BSs.size(); i++) {
        // No need to calculate (ji|kl)
        for(auto j=0; j<i; j++) {
            // No need to calculate (kl|ij)...
            for(auto k=0; k<i; k++) {
                // No need to calculate (kl|ji)...
                for(auto l=0; l<k; l++) {
                    //h2e_engine.compute(BSs[i], BSs[j], BSs[k], BSs[l]);
                    //std::cout << h2e_engine<<endl;
                }
            }
        }
    }

}

// Get everything in a single command
void Eint::compute() {

    // Initialize libint
    libint2::initialize();

    // H1E
    this->compute_ovlp();
    this->compute_h1e_kin();
    this->compute_h1e_nuc();
    this->Mol_pointer->HCore = this->Mol_pointer->HVel 
                             + this->Mol_pointer->HNuc;
    // H2E
    this->compute_eri();

    // End libint 
    libint2::finalize();

}
