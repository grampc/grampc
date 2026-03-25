/* This file is part of GRAMPC - (https://github.com/grampc/grampc)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright 2014-2025 by Knut Graichen, Andreas Voelz, Thore Wietzke,
 * Tobias Englert (<v2.3), Felix Mesmer (<v2.3), Soenke Rhein (<v2.3),
 * Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * All rights reserved.
 *
 * GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */
#include "mobile_robot_problem_description.hpp"

/** Utility macros for integer powers */
#define POW2(x) (x)*(x)
#define POW3(x) (x)*(x)*(x)
#define POW4(x) (x)*(x)*(x)*(x)


MobileRobotProblemDescription::MobileRobotProblemDescription(const std::vector<typeRNum>& pSys, const std::vector<typeRNum>& pCost, const std::vector<typeRNum>& pCon)
: ProblemDescription(3, 2, 0, 0, 4, 0, 0),
  pSys_(pSys),
  pCost_(pCost),
  pCon_(pCon)
{
}

void MobileRobotProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{    
    out[0] = u[0] * cos(x[2]);
    out[1] = u[0] * sin(x[2]);
    out[2] = u[1];
}

void MobileRobotProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{
    out[0] = 0;
    out[1] = 0;
    out[2] = u[0] * (vec[1] * cos(x[2]) - vec[0] * sin(x[2]));
}

void MobileRobotProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{
    out[0] = vec[0] * cos(x[2]) + vec[1] * sin(x[2]);
    out[1] = vec[2];
}

void MobileRobotProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{
    // desired trajectory for position and orientation
    typeRNum xd = cos(*param.t0 + t);
    typeRNum yd = sin(*param.t0 + t);
    typeRNum thetad = M_PI_2 + *param.t0 + t;
    // longitudinal and lateral errors
    typeRNum elon =  cos(thetad) * (x[0] - xd) + sin(thetad) * (x[1] - yd);
    typeRNum elat = -sin(thetad) * (x[0] - xd) + cos(thetad) * (x[1] - yd);
    // non-quadratic cost function
    out[0] = pCost_[0] * POW4(elon)
           + pCost_[1] * POW2(elat)
           + pCost_[2] * POW4(sin((x[2] - thetad) / 2))
           + pCost_[3] * POW4(u[0] - 1.0)
           + pCost_[4] * POW4(u[1] - 1.0);
}

void MobileRobotProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{
    // desired trajectory for position and orientation
    typeRNum xd = cos(*param.t0 + t);
    typeRNum yd = sin(*param.t0 + t);
    typeRNum thetad = M_PI_2 + *param.t0 + t;
    // longitudinal and lateral errors
    typeRNum elon =  cos(thetad) * (x[0] - xd) + sin(thetad) * (x[1] - yd);
    typeRNum elat = -sin(thetad) * (x[0] - xd) + cos(thetad) * (x[1] - yd);
    // gradient of non-quadratic cost function
    out[0] = pCost_[0] * 4 * POW3(elon) * cos(thetad) + pCost_[1] * 2 * elat * -sin(thetad);
    out[1] = pCost_[0] * 4 * POW3(elon) * sin(thetad) + pCost_[1] * 2 * elat * cos(thetad);
    out[2] = pCost_[2] * 4 * POW3(sin((x[2] - thetad) / 2)) * (cos((x[2] - thetad) / 2) / 2);
}

void MobileRobotProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{
    // gradient of non-quadratic cost function
    out[0] = pCost_[3] * 4 * POW3(u[0] - 1.0);
    out[1] = pCost_[4] * 4 * POW3(u[1] - 1.0);
}

void MobileRobotProblemDescription::Vfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
{
    out[0] = 0;
}

void MobileRobotProblemDescription::dVdx(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
{
    out[0] = 0;
    out[1] = 0;
    out[2] = 0;
}

void MobileRobotProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{
    out[0] = pCon_[0] - x[0]; // min x-position
    out[1] = x[0] - pCon_[1]; // max x-position
    out[2] = pCon_[2] - x[1]; // min y-position
    out[3] = x[1] - pCon_[3]; // max y-position
}

void MobileRobotProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{
    out[0] = vec[1] - vec[0];
    out[1] = vec[3] - vec[2];
    out[2] = 0;
}

void MobileRobotProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{
    out[0] = 0;
    out[1] = 0;
}
