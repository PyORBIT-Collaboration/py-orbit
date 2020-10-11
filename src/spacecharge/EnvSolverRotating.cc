#include "EnvSolverRotating.hh"

EnvSolverRotating::EnvSolverRotating(double perveance): CppPyWrapper(NULL)
{
    Q = perveance;
}
    
void EnvSolverRotating::trackBunch(Bunch* bunch, double length)
{
    double a, b, e, f;
    double phi, cosP, sinP, cosP2, sinP2;
    double cx, cy, factor;

    a = bunch->x(0);
    b = bunch->x(1);
    e = bunch->y(0);
    f = bunch->y(1);

    phi = -0.5 * atan2(2*(a*e + b*f), (a*a + b*b - e*e - f*f));
    cosP = cos(phi);
    sinP = sin(phi);
    cosP2 = cosP*cosP;
    sinP2 = sinP*sinP;
    
    cx = sqrt(pow(a*f-b*e, 2) / ((e*e+f*f)*cosP2 + (a*a+b*b)*sinP2 +  2*(a*e+b*f)*cosP*sinP));
    cy = sqrt(pow(a*f-b*e, 2) / ((a*a+b*b)*cosP2 + (e*e+f*f)*sinP2 -  2*(a*e+b*f)*cosP*sinP));
    factor = 2 * Q / (cx + cy);

    bunch->xp(0) += (factor * ((a*cosP2 - e*sinP*cosP)/cx + (a*sinP2 + e*sinP*cosP)/cy)) * length;
    bunch->xp(1) += (factor * ((b*cosP2 - f*sinP*cosP)/cx + (b*sinP2 + f*sinP*cosP)/cy)) * length;
    bunch->yp(0) += (factor * ((e*cosP2 + a*sinP*cosP)/cy + (e*sinP2 - a*sinP*cosP)/cx)) * length;
    bunch->yp(1) += (factor * ((f*cosP2 + b*sinP*cosP)/cy + (f*sinP2 - b*sinP*cosP)/cx)) * length;
}
