#include "mol.hpp"
#include "eint.hpp"
#include "controls.hpp"
using namespace std;
  
int main()
{
    // Read Controls
    readControls("Input");

    // Construct Molecule Object
    Molecule* QMole = new Molecule();

    // Get Basic information
    QMole->read_geom("Input");
    QMole->print_xyz();
    QMole->calc_enuc();

    // Setup Basis set
    Eint* Integrals = new Eint(QMole);
    Integrals->read_int("Input");
    Integrals->build_shell();

    // Compute
    Integrals->compute();

    return 0;
}