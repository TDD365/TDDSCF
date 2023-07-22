#include "controls.hpp"
#include "hartreefock.hpp"
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

    // HF
    HartreeFock* SCF = new HartreeFock(Integrals);
    SCF->do_SCF();

    return 0;
}