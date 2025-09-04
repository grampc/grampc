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
#ifndef PYGRAMPC_TYPES_HPP
#define PYGRAMPC_TYPES_HPP

extern "C" 
{
    #include "grampc.h"
}
#include <Eigen/Dense>
#include <vector>

namespace grampc
{
    typedef Eigen::Matrix<typeRNum, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Matrix<typeRNum, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Ref<Vector> VectorRef;
    typedef const Eigen::Ref<const Vector>& VectorConstRef;

    struct GrampcParam
    {
        GrampcParam();
        GrampcParam(const typeGRAMPCparam* param);
        void remap_memory(const typeGRAMPC* grampc);

        const int* Nx;
        const int* Nu;
        const int* Np;
        const int* Ng;
        const int* Nh;
        const int* NgT;
        const int* NhT;
        const int* Nc;

        Eigen::Map<Vector> x0;
        Eigen::Map<Vector> xdes;

        Eigen::Map<Vector> u0;
        Eigen::Map<Vector> udes;
        Eigen::Map<Vector> umax;
        Eigen::Map<Vector> umin;

        Eigen::Map<Vector> p0;
        Eigen::Map<Vector> pmax;
        Eigen::Map<Vector> pmin;

        const typeRNum* Thor;
        const typeRNum* Tmax;
        const typeRNum* Tmin;

        const typeRNum* dt;
        const typeRNum* t0;
    };

    
    /** Options for GRAMPC computations, which are part of the GRAMPC data struct.*/
    struct GrampcOpt
    {
        GrampcOpt();
        void remap_memory(const typeGRAMPC* grampc);

        int* Nhor;
        int* MaxGradIter;
        int* MaxMultIter;
        int* ShiftControl;

        int* TimeDiscretization;

        int* IntegralCost;
        int* TerminalCost;
        int* IntegratorCost;

        int* Integrator;
        typeRNum* IntegratorRelTol;
        typeRNum* IntegratorAbsTol;
        typeRNum* IntegratorMinStepSize;
        int*    IntegratorMaxSteps;
        std::vector<int> FlagsRodas;

        int* LineSearchType;
        int* LineSearchExpAutoFallback;
        typeRNum* LineSearchMax;
        typeRNum* LineSearchMin;
        typeRNum* LineSearchInit;
        typeRNum* LineSearchAdaptAbsTol;
        typeRNum* LineSearchAdaptFactor;
        typeRNum* LineSearchIntervalTol;
        typeRNum* LineSearchIntervalFactor;

        int* OptimControl;
        int* OptimParam;
        typeRNum* OptimParamLineSearchFactor;
        int* OptimTime;
        typeRNum* OptimTimeLineSearchFactor;

        int* ScaleProblem;
        Eigen::Map<Vector> xScale;
        Eigen::Map<Vector> xOffset;
        Eigen::Map<Vector> uScale;
        Eigen::Map<Vector> uOffset;
        Eigen::Map<Vector> pScale;
        Eigen::Map<Vector> pOffset;
        typeRNum* TScale;
        typeRNum* TOffset;
        typeRNum* JScale;
        Eigen::Map<Vector> cScale;

        int* EqualityConstraints;
        int* InequalityConstraints;
        int* TerminalEqualityConstraints;
        int* TerminalInequalityConstraints;
        int* ConstraintsHandling;
        Eigen::Map<Vector> ConstraintsAbsTol;

        typeRNum* MultiplierMax;
        typeRNum* MultiplierDampingFactor;
        typeRNum* PenaltyMax;
        typeRNum* PenaltyMin;
        typeRNum* PenaltyIncreaseFactor;
        typeRNum* PenaltyDecreaseFactor;
        typeRNum* PenaltyIncreaseThreshold;
        typeRNum* AugLagUpdateGradientRelTol;

        int* ConvergenceCheck;
        typeRNum* ConvergenceGradientRelTol;
    };

    /** Solution structure of GRAMPC computations, which is part of the GRAMPC data struct. */
    struct GrampcSol
    {
        GrampcSol();
        void remap_memory(const typeGRAMPC* grampc);

        Eigen::Map<Vector> xnext;
        Eigen::Map<Vector> unext;
        Eigen::Map<Vector> pnext;
        typeRNum* Tnext;
        Eigen::Map<Vector> J;
        typeRNum* cfct;
        typeRNum* pen;
        std::vector<int> iter;
        int status;
    };

    /** Real workspace structure of GRAMPC that holds intermediate results
    and computation results. */
    struct GrampcRWS
    {
        GrampcRWS();
        void remap_memory(const typeGRAMPC* grampc);

        Eigen::Map<Vector> t;
        Eigen::Map<Vector> tls;

        Eigen::Map<Matrix> x;
        Eigen::Map<Matrix> adj;
        Eigen::Map<Matrix> dcdx;

        Eigen::Map<Matrix> u;
        Eigen::Map<Matrix> uls;
        Eigen::Map<Matrix> uprev;
        Eigen::Map<Matrix> gradu;
        Eigen::Map<Matrix> graduprev;
        Eigen::Map<Matrix> dcdu;

        Eigen::Map<Vector> p;
        Eigen::Map<Vector> pls;
        Eigen::Map<Vector> pprev;
        Eigen::Map<Vector> gradp;
        Eigen::Map<Vector> gradpprev;
        Eigen::Map<Matrix> dcdp;

        typeRNum* T;
        typeRNum* Tprev;
        typeRNum* gradT;
        typeRNum* gradTprev;
        typeRNum* dcdt;

        Eigen::Map<Matrix> mult;
        Eigen::Map<Matrix> pen;
        Eigen::Map<Matrix> cfct;
        Eigen::Map<Matrix> cfctprev;
        Eigen::Map<Vector> cfctAbsTol;

        Eigen::Map<Vector> lsAdapt;
        Eigen::Map<Vector> lsExplicit;
        Eigen::Map<Vector> rwsScale;
        int* lrwsGeneral;
        Eigen::Map<Vector> rwsGeneral;

        int* lworkRodas;
        int* liworkRodas;
        Eigen::Map<Vector> rparRodas;
        std::vector<int> iparRodas;
        Eigen::Map<Vector> workRodas;
        std::vector<int> iworkRodas;
    };
}

#endif