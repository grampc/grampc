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
#ifndef PYPROBLEM_DESCRIPTION_HPP
#define PYPROBLEM_DESCRIPTION_HPP

extern "C"
{
    #include "grampc.h"
}

#include "pygrampc_types.hpp"
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace grampc
{
    class PyProblemDescription
    {
        public:
            ctypeInt Nx;
            ctypeInt Nu;
            ctypeInt Np;
            ctypeInt Ng;
            ctypeInt Nh;
            ctypeInt NgT;
            ctypeInt NhT;
            typeInt Rodas_Jac = 0;
            typeInt Rodas_M = 0;

        public:
            PyProblemDescription(ctypeInt Nx, ctypeInt Nu, ctypeInt Np, ctypeInt Ng, ctypeInt Nh, ctypeInt NgT, ctypeInt NhT) 
             : Nx(Nx), Nu(Nu), Np(Np), Ng(Ng), Nh(Nh), NgT(NgT), NhT(NhT)
            {}

            virtual ~PyProblemDescription() {}

            /** System function f(t,x,u,p,userparam)
            ------------------------------------ **/
            virtual void ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) = 0;
            /** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
            virtual void dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) = 0;
            /** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
            virtual void dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}
            /** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
            virtual void dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}


            /** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam)
            -------------------------------------------------- **/
            virtual void lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) {}
            /** Gradient dl/dx **/
            virtual void dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) {}
            /** Gradient dl/du **/
            virtual void dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) {}
            /** Gradient dl/dp **/
            virtual void dldp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) {}


            /** Terminal cost V(T,x,p) */
            virtual void Vfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const GrampcParam& param) {}
            /** Gradient dV/dx **/
            virtual void dVdx(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const GrampcParam& param) {}
            /** Gradient dV/dp **/
            virtual void dVdp(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const GrampcParam& param) {}
            /** Gradient dV/dT **/
            virtual void dVdT(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const GrampcParam& param) {}


            /** Equality constraints g(t,x,u,p) = 0 */
            virtual void gfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) {}
            /** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
            virtual void dgdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}
            /** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
            virtual void dgdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}
            /** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
            virtual void dgdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}


            /** Inequality constraints h(t,x,u,p) < 0 */
            virtual void hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) {}
            /** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dh/dx) **/
            virtual void dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}
            /** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dh/du) **/
            virtual void dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}
            /** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dh/dp) **/
            virtual void dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}


            /** Terminal equality constraints gT(T,x,p) = 0 */
            virtual void gTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const GrampcParam& param) {}
            /** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
            virtual void dgTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}
            /** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
            virtual void dgTdp_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}
            /** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
            virtual void dgTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}


            /** Terminal inequality constraints hT(T,x,p) < 0 */
            virtual void hTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const GrampcParam& param) {}
            /** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
            virtual void dhTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}
            /** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
            virtual void dhTdp_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}
            /** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
            virtual void dhTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}


            /** Additional functions required for semi-implicit systems
            M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS
            ------------------------------------------------------- **/
            /** Jacobian df/dx in vector form (column-wise) **/
            virtual void dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) {}
            /** Jacobian df/dx in vector form (column-wise) **/
            virtual void dfdxtrans(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) {}
            /** Jacobian df/dt **/
            virtual void dfdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) {}
            /** Jacobian d(dH/dx)/dt  **/
            virtual void dHdxdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) {}
            /** Mass matrix in vector form (column-wise, either banded or full matrix) **/
            virtual void Mfct(VectorRef out, const GrampcParam& param) {}
            /** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
            virtual void Mtrans(VectorRef out, const GrampcParam& param) {}
    };
    
    // Alias
    typedef std::shared_ptr<PyProblemDescription> ProblemDescriptionPtr;
    typedef std::shared_ptr<const PyProblemDescription> ProblemDescriptionConstPtr;

    /* Template class needed for overriting C++ functions in python*/
    class PyProblem : public PyProblemDescription
    {
        public:
            /* Inherit the constructor*/
            using PyProblemDescription::PyProblemDescription;

            void ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE_PURE(
                    void, 
                    PyProblemDescription, 
                    ffct,
                    out, t, x, u, p, param);
            }

            void dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE_PURE(
                    void, 
                    PyProblemDescription, 
                    dfdx_vec,
                    out, t, x, u, p, vec, param);
            }

            void dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dfdu_vec,
                    out, t, x, u, p, vec, param);
            }

            void dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dfdp_vec,
                    out, t, x, u, p, vec, param);
            }

            void lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    lfct,
                    out, t, x, u, p, param);
            }

            void dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dldx,
                    out, t, x, u, p, param);
            }

            void dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dldu,
                    out, t, x, u, p, param);
            }
            
            void dldp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dldp,
                    out, t, x, u, p, param);
            }

            void Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    Vfct,
                    out, t, x, p, param);
            }

            void dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dVdx,
                    out, t, x, p, param);
            }

            void dVdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dVdp,
                    out, t, x, p, param);
            }
            
            void dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dVdT,
                    out, t, x, p, param);
            }

            void gfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    gfct,
                    out, t, x, u, p, param);
            }

            void dgdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dgdx_vec,
                    out, t, x, u, p, vec, param);
            }

            void dgdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dgdu_vec,
                    out, t, x, u, p, vec, param);
            }
            void dgdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dgdp_vec,
                    out, t, x, u, p, vec, param);
            }

            void hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    hfct,
                    out, t, x, u, p, param);
            }

            void dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dhdx_vec,
                    out, t, x, u, p, vec, param);
            }

            void dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dhdu_vec,
                    out, t, x, u, p, vec, param);
            }
            void dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dhdp_vec,
                    out, t, x, u, p, vec, param);
            }

            void gTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    gTfct,
                    out, t, x, p, param);
            }

            void dgTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dgTdx_vec,
                    out, t, x, p, vec, param);
            }

            void dgTdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dgTdp_vec,
                    out, t, x, p, vec, param);
            }

            void dgTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dgTdT_vec,
                    out, t, x, p, vec, param);
            }

            void hTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    hTfct,
                    out, t, x, p, param);
            }

            void dhTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dhTdx_vec,
                    out, t, x, p, vec, param);
            }

            void dhTdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dhTdp_vec,
                    out, t, x, p, vec, param);
            }

            void dhTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dhTdT_vec,
                    out, t, x, p, vec, param);
            }

            void dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dfdx,
                    out, t, x, u, p, param);
            }
            
            void dfdxtrans(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dfdxtrans,
                    out, t, x, u, p, param);
            }
            
            void dfdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dfdt,
                    out, t, x, u, p, param);
            }
            
            void dHdxdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    dHdxdt,
                    out, t, x, u, p, vec, param);
            }
            
            void Mfct(VectorRef out, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    Mfct,
                    out, param);
            }
            void Mtrans(VectorRef out, const GrampcParam& param) 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    PyProblemDescription, 
                    Mtrans,
                    out, param);
            }
    };
}

#endif
