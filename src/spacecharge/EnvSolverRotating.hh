/** Envelope solver for the 2D rotating distribution.
 
 Reference: 
 V. Danilov, S. Cousineau, S. Henderson, and J. Holmes, Self-consistent time dependent
 two dimensional and three dimensional space charge distributions with linear force,
 Physical Review Special Topics - Accelerators and Beams 6, 74–85 (2003).
*/

#ifndef ENVSOLVER_ROTATING_H
#define ENVSOLVER_ROTATING_H

#include "Bunch.hh"
#include "CppPyWrapper.hh"

using namespace std;

class EnvSolverRotating: public OrbitUtils::CppPyWrapper
{
    public:
    
    double Q; // beam perveance

    /** Constructor */
    EnvSolverRotating(double perveance);

    /** Apply space charge kick.
     
     The boundary of a tilted ellipse in real space can be parameterized as:
        x = a*cos(psi) + b*sin(psi),
        y = e*cos(psi) + f*sin(psi),
     where 0 <= psi <= 2pi. The method tracks two particles. The particle coordinates
     map to the parameters {a, b, e, f} and their slopes in the following way:
        particle 1 --> {x:a, x':a', y:e, y':e'},
        particle 2 --> {x:e, x':e', y:f, y':f'}.
    */
    void trackBunch(Bunch* bunch, double length);
};
//end of ENVSOLVER_ROTATING_H
#endif






















/**---------------------------------------------------------------------------------*/




// #ifndef ENVSOLVER_H
// #define ENVSOLVER_H

// #include "Bunch.hh"

// namespace envelopesolver
// {
//     /** Apply space charge kick according to the KV envelope equations.
     
//      The method tracks a single particle who's x and y coordinates represent
//      the horizontal and vertical beam radii of a uniform density, upright ellipse.

//      Reference:
//      S. M. Lund and B. Bukh, Phys. Rev. ST Accel. Beams 7, 024801 (2004).

//     * @param bunch: Bunch object
//     * @param ex{y}: r.m.s. emittance in the x{y} dimension
//     * @param Q: beam perveance
//     * @param length: length of the applied kick
//     */
//     void KVSolver(Bunch* bunch, double ex, double ey, double Q, double length);

    
//     /** Apply space charge kick according to rotating envelope equations.
//     *
//     * The coordinates of the boundary of a tilted ellipse can be parameterized as
//     *     x = a * cos(psi) + b * sin(psi),
//     *     y = e * cos(psi) + f * sin(psi),
//     * with 0 <= psi <= 2pi. The method tracks two particles. The coordinates relate
//     * to the a, b, e, and f in the following way:
//     *     particle 1 -> {x:a, x':a', y:e, y':e'},
//     *     particle 2 -> {x:b, x':b', y:f, y':f'}.
//     * The kicks are applied assuming a uniform particle density within the ellipse.
//     * 
//     * Reference: 
//     *  V. Danilov, S. Cousineau, S. Henderson, and J. Holmes, Self-consistent time 
//     *  dependent two dimensional and three dimensional space charge distributions 
//     *  with linear force, Physical Review Special Topics - Accelerators and Beams 
//     *  6, 74–85 (2003).
//     *
//     * @param bunch: Bunch object
//     * @param Q: beam perveance
//     * @param length: length of the applied kick
//     */
//     void RotatingSolver(Bunch* bunch, double Q, double length);

// }

// #endif // ENVSOLVER_H

















