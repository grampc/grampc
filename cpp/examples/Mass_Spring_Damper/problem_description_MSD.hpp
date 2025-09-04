/* This file is part of GRAMPC - (https://github.com/grampc/grampc)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright 2014-2025 by Knut Graichen, Andreas Voelz, Thore Wietzke,
 * Tobias Englert (<v2.3), Felix Mesmer (<v2.3), Soenke Rhein (<v2.3),
 * Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 *
 * GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 *
 *
 *
 *
 *
 *
 *
 * This probfct-file describes the mass-spring-damper problem from
 * Kapernick, B.: Gradient-Based Nonlinear Model Predictive Control With
 * Constraint Transformation for Fast Dynamical Systems. Dissertation,
 * Ulm University. Shaker, Aachen, Germany (2016)
 *
 *                                           _T
 *		                                    /
 *      min    J(u,p,T;x0) = V(T,x(T),p) + / l(t,x(t),u(t),p) dt
 *   u(.),p,T                            _/
 *                                      0
 *             .
 *      s.t.   x(t) = f(t0+t,x(t),u(t),p), x(0) = x0
 *             u_min <= u(t) <= u_max
 *
 */

#ifndef MSD_HPP
#define MSD_HPP

#include <vector>
#include "problem_description.hpp"

 /* square macro */
#define POW2(a) ((a)*(a))

class MassSpringDamper : public grampc::ProblemDescription
{
    public:
        MassSpringDamper(typeInt NN, typeRNum m, typeRNum c, typeRNum d, 
                const std::vector<typeRNum>& Q, const std::vector<typeRNum>& R, const std::vector<typeRNum>& P);

        virtual void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT) override;

        /** System function f(t,x,u,p,param,userparam) **/
        virtual void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param) override;
        /** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
        virtual void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param) override;
        /** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
        virtual void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param) override;
        
        /** Integral cost l(t,x(t),u(t),p,param,userparam) **/
        virtual void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param) override;
        /** Gradient dl/dx **/
        virtual void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param) override;
        /** Gradient dl/du **/
        virtual void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param) override;
        
        /** Terminal cost V(T,x,p) */
        virtual void Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param) override;
        /** Gradient dV/dx **/
        virtual void dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param) override;

    private:
        ctypeInt NN_; //Number of discretization points
        ctypeInt Nx_;
        ctypeInt Nu_;
        ctypeRNum m_;
        ctypeRNum c_;
        ctypeRNum d_;
        std::vector<typeRNum> Q_;
        std::vector<typeRNum> R_;
        std::vector<typeRNum> P_;
};

#endif