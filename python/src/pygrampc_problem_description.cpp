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
#include "pygrampc_problem_description.hpp"

extern "C"
{
    void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        *Nx = problem->Nx;
        *Nu = problem->Nu;
        *Np = problem->Np;
        *Ng = problem->Ng;
        *Nh = problem->Nh;
        *NgT = problem->NgT;
        *NhT = problem->NhT;
    }

    void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);

        problem->ffct(outMap, t, xMap, uMap, pMap, paramMap);
    }
    
    void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dfdx_vec(outMap, t, xMap, uMap, pMap, vecMap, paramMap);
    }

    void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nu);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dfdu_vec(outMap, t, xMap, uMap, pMap, vecMap, paramMap);
    }

    void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dfdp_vec(outMap, t, xMap, uMap, pMap, vecMap, paramMap);
    }

    void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, 1);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->lfct(outMap, t, xMap, uMap, pMap, paramMap);
    }

    void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dldx(outMap, t, xMap, uMap, pMap, paramMap);
    }

    void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nu);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dldu(outMap, t, xMap, uMap, pMap, paramMap);
    }

    void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dldp(outMap, t, xMap, uMap, pMap, paramMap);
    }

    void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, 1);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->Vfct(outMap, T, xMap, pMap, paramMap);
    }

    void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dVdx(outMap, T, xMap, pMap, paramMap);
    }

    void dVdp(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dVdp(outMap, T, xMap, pMap, paramMap);
    }

    void dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, 1);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dVdT(outMap, T, xMap, pMap, paramMap);
    }

    void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Ng);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->gfct(outMap, t, xMap, uMap, pMap, paramMap);
    }

    void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Ng);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dgdx_vec(outMap, t, xMap, uMap, pMap, vecMap, paramMap);
    }

    void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nu);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Ng);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dgdu_vec(outMap, t, xMap, uMap, pMap, vecMap, paramMap);
    }

    void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Ng);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dgdp_vec(outMap, t, xMap, uMap, pMap, vecMap, paramMap);
    }

    void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nh);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->hfct(outMap, t, xMap, uMap, pMap, paramMap);
    }

    void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nh);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dhdx_vec(outMap, t, xMap, uMap, pMap, vecMap, paramMap);
    }

    void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nu);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nh);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dhdu_vec(outMap, t, xMap, uMap, pMap, vecMap, paramMap);
    }

    void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nh);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dhdp_vec(outMap, t, xMap, uMap, pMap, vecMap, paramMap);
    }

    void gTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->NgT);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->gTfct(outMap, T, xMap, pMap, paramMap);
    }

    void dgTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->NgT);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dgTdx_vec(outMap, T, xMap, pMap, vecMap, paramMap);
    }

    void dgTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->NgT);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dgTdp_vec(outMap, T, xMap, pMap, vecMap, paramMap);
    }

    void dgTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, 1);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->NgT);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dgTdT_vec(outMap, T, xMap, pMap, vecMap, paramMap);
    }

    void hTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->NhT);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->hTfct(outMap, T, xMap, pMap, paramMap);
    }

    void dhTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->NhT);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dhTdx_vec(outMap, T, xMap, pMap, vecMap, paramMap);
    }

    void dhTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->NhT);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);

        problem->dhTdp_vec(outMap, T, xMap, pMap, vecMap, paramMap);
    }

    void dhTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, 1);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->NhT);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dhTdT_vec(outMap, T, xMap, pMap, vecMap, paramMap);
    }

    void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Rodas_Jac); //Nx * (MLJAC + MUJAC + 1)
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dfdx(outMap, t, xMap, uMap, pMap, paramMap);
    }

    void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Rodas_Jac); //Nx * (MLJAC + MUJAC + 1)
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dfdxtrans(outMap, t, xMap, uMap, pMap, paramMap);
    }

    void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dfdt(outMap, t, xMap, uMap, pMap, paramMap);
    }

    void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        grampc::GrampcParam paramMap(param);
        
        problem->dHdxdt(outMap, t, xMap, uMap, pMap, vecMap, paramMap);
    }

    void Mfct(typeRNum *out, const typeGRAMPCparam *param, typeUSERPARAM *userparam) // Auf Python Seite eine Matrix?
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Rodas_M); //Nx * (MLMAS + MUMAS + 1)
        grampc::GrampcParam paramMap(param);
        
        problem->Mfct(outMap, paramMap);
    }

    void Mtrans(typeRNum *out, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
    {
        grampc::PyProblemDescription* problem = (grampc::PyProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Rodas_M); //Nx * (MLMAS + MUMAS + 1)
        grampc::GrampcParam paramMap(param);
        
        problem->Mtrans(outMap, paramMap);
    }

}
