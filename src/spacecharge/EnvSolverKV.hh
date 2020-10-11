/** KV envelope solver

 Reference: 
 S. M. Lund and B. Bukh, Phys. Rev. ST Accel. Beams 7, 024801 (2004).
*/

#ifndef ENVSOLVER_KV_H
#define ENVSOLVER_KV_H

#include "Bunch.hh"
#include "CppPyWrapper.hh"

using namespace std;

class EnvSolverKV: public OrbitUtils::CppPyWrapper
{
    public:
    
    double ex, ey, Q;

    /** Constructor */
    EnvSolverKV(double eps_x, double eps_y, double perveance);

    /** Apply space charge kick.
     
     The method tracks a single particle who's x and y coordinates represent
     the horizontal and vertical beam radii of the ellipse.
    */
    void trackBunch(Bunch* bunch, double length);
};
//end of ENVSOLVER_KV_H
#endif
