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
#ifndef MOBILE_ROBOT_PROBLEM_DESCRIPTION_HPP
#define MOBILE_ROBOT_PROBLEM_DESCRIPTION_HPP

#include "problem_description.hpp"

using namespace grampc;

class MobileRobotProblemDescription : public ProblemDescription
{
public:
    MobileRobotProblemDescription(const std::vector<typeRNum>& pSys, const std::vector<typeRNum>& pCost, const std::vector<typeRNum>& pCon);

    virtual void ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override;
    virtual void dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override;
    virtual void dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override;

    virtual void lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override;
    virtual void dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override;
    virtual void dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override;

    virtual void Vfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override;
    virtual void dVdx(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override;

    virtual void hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override;
    virtual void dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override;
    virtual void dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override;

private:
    // system parameters
    std::vector<typeRNum> pSys_;

    // cost parameters
    std::vector<typeRNum> pCost_;

    // constraint parameters
    std::vector<typeRNum> pCon_;
};

#endif // MOBILE_ROBOT_PROBLEM_DESCRIPTION_HPP