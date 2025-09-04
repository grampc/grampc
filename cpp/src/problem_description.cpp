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


#include "problem_description.hpp"


extern "C"
{
	void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->ocp_dim(Nx, Nu, Np, Ng, Nh, NgT, NhT);
	}


	void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->ffct(out, t, x, u, p, param);
	}
	void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dfdx_vec(out, t, x, u, p, vec, param);
	}
	void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dfdu_vec(out, t, x, u, p, vec, param);
	}
	void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dfdp_vec(out, t, x, u, p, vec, param);
	}


	void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->lfct(out, t, x, u, p, param);
	}
	void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dldx(out, t, x, u, p, param);
	}
	void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dldu(out, t, x, u, p, param);
	}
	void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dldp(out, t, x, u, p, param);
	}


	void Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->Vfct(out, t, x, p, param);
	}
	void dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dVdx(out, t, x, p, param);
	}
	void dVdp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dVdp(out, t, x, p, param);
	}
	void dVdT(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dVdT(out, t, x, p, param);
	}


	void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->gfct(out, t, x, u, p, param);
	}
	void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->dgdx_vec(out, t, x, u, p, vec, param);
	}
	void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dgdu_vec(out, t, x, u, p, vec, param);
	}
	void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dgdp_vec(out, t, x, u, p, vec, param);
	}


	void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->hfct(out, t, x, u, p, param);
	}
	void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhdx_vec(out, t, x, u, p, vec, param);
	}
	void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhdu_vec(out, t, x, u, p, vec, param);
	}
	void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhdp_vec(out, t, x, u, p, vec, param);
	}


	void gTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->gTfct(out, t, x, p, param);
	}
	void dgTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dgTdx_vec(out, t, x, p, vec, param);
	}
	void dgTdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dgTdp_vec(out, t, x, p, vec, param);
	}
	void dgTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dgTdT_vec(out, t, x, p, vec, param);
	}


	void hTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->hTfct(out, t, x, p, param);
	}
	void dhTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhTdx_vec(out, t, x, p, vec, param);
	}
	void dhTdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhTdp_vec(out, t, x, p, vec, param);
	}
	void dhTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhTdT_vec(out, t, x, p, vec, param);
	}


	void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhTdT_vec(out, t, x, u, p, param);
	}
	void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhTdT_vec(out, t, x, u, p, param);
	}
	void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhTdT_vec(out, t, x, u, p, param);
	}
	void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->dHdxdt(out, t, x, u, p, vec, param);
	}
	void Mfct(typeRNum *out, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->Mfct(out, param);
	}
	void Mtrans(typeRNum *out, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->Mtrans(out, param);
	}

}
