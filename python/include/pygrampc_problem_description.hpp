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

#include "problem_description.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace grampc
{
    /* Template class needed for overriting C++ functions in python*/
    class PyProblem : public ProblemDescription
    {
        public:
            /* Inherit the constructor*/
            using ProblemDescription::ProblemDescription;

            void ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE_PURE(
                    void, 
                    ProblemDescription, 
                    ffct,
                    out, t, x, u, p, param);
            }

            void dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE_PURE(
                    void, 
                    ProblemDescription, 
                    dfdx_vec,
                    out, t, x, u, p, vec, param);
            }

            void dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dfdu_vec,
                    out, t, x, u, p, vec, param);
            }

            void dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dfdp_vec,
                    out, t, x, u, p, vec, param);
            }

            void lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    lfct,
                    out, t, x, u, p, param);
            }

            void dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dldx,
                    out, t, x, u, p, param);
            }

            void dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dldu,
                    out, t, x, u, p, param);
            }
            
            void dldp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dldp,
                    out, t, x, u, p, param);
            }

            void Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    Vfct,
                    out, t, x, p, param);
            }

            void dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dVdx,
                    out, t, x, p, param);
            }

            void dVdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dVdp,
                    out, t, x, p, param);
            }
            
            void dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dVdT,
                    out, t, x, p, param);
            }

            void gfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    gfct,
                    out, t, x, u, p, param);
            }

            void dgdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dgdx_vec,
                    out, t, x, u, p, vec, param);
            }

            void dgdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dgdu_vec,
                    out, t, x, u, p, vec, param);
            }
            void dgdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dgdp_vec,
                    out, t, x, u, p, vec, param);
            }

            void hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    hfct,
                    out, t, x, u, p, param);
            }

            void dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dhdx_vec,
                    out, t, x, u, p, vec, param);
            }

            void dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dhdu_vec,
                    out, t, x, u, p, vec, param);
            }
            void dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dhdp_vec,
                    out, t, x, u, p, vec, param);
            }

            void gTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    gTfct,
                    out, t, x, p, param);
            }

            void dgTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dgTdx_vec,
                    out, t, x, p, vec, param);
            }

            void dgTdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dgTdp_vec,
                    out, t, x, p, vec, param);
            }

            void dgTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dgTdT_vec,
                    out, t, x, p, vec, param);
            }

            void hTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    hTfct,
                    out, t, x, p, param);
            }

            void dhTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dhTdx_vec,
                    out, t, x, p, vec, param);
            }

            void dhTdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dhTdp_vec,
                    out, t, x, p, vec, param);
            }

            void dhTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dhTdT_vec,
                    out, t, x, p, vec, param);
            }

            void dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dfdx,
                    out, t, x, u, p, param);
            }
            
            void dfdxtrans(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dfdxtrans,
                    out, t, x, u, p, param);
            }
            
            void dfdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dfdt,
                    out, t, x, u, p, param);
            }
            
            void dHdxdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dHdxdt,
                    out, t, x, u, p, vec, param);
            }
            
            void Mfct(VectorRef out, const GrampcParam& param) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    Mfct,
                    out, param);
            }
            void Mtrans(VectorRef out, const GrampcParam& param) 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    Mtrans,
                    out, param);
            }
    };
}

#endif
