#include "EnvSolverKV.hh"

EnvSolverKV::EnvSolverKV(double eps_x, double eps_y, double perveance): CppPyWrapper(NULL)
{ 
    ex = eps_x;
    ey = eps_y;
    Q = perveance;
}
    
void EnvSolverKV::trackBunch(Bunch* bunch, double length)
{
    double a = std::abs(bunch->x(0));
    double b = std::abs(bunch->y(0));

    bunch->xp(0) += ((16*ex*ex / (a*a*a)) + 2*Q/(a+b)) * length;
    bunch->yp(0) += ((16*ey*ey / (b*b*b)) + 2*Q/(a+b)) * length;
}
