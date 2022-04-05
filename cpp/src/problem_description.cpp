/* This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright 2014-2019 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
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


	void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->ffct(out, t, x, u, p);
	}
	void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dfdx_vec(out, t, x, adj, u, p);
	}
	void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dfdu_vec(out, t, x, adj, u, p);
	}
	void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dfdp_vec(out, t, x, adj, u, p);
	}


	void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->lfct(out, t, x, u, p, xdes, udes);
	}
	void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dldx(out, t, x, u, p, xdes, udes);
	}
	void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dldu(out, t, x, u, p, xdes, udes);
	}
	void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dldp(out, t, x, u, p, xdes, udes);
	}


	void Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->Vfct(out, t, x, p, xdes);
	}
	void dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dVdx(out, t, x, p, xdes);
	}
	void dVdp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dVdp(out, t, x, p, xdes);
	}
	void dVdT(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dVdT(out, t, x, p, xdes);
	}


	void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->gfct(out, t, x, u, p);
	}
	void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->dgdx_vec(out, t, x, u, p, vec);
	}
	void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dgdu_vec(out, t, x, u, p, vec);
	}
	void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dgdp_vec(out, t, x, u, p, vec);
	}


	void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->hfct(out, t, x, u, p);
	}
	void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhdx_vec(out, t, x, u, p, vec);
	}
	void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhdu_vec(out, t, x, u, p, vec);
	}
	void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhdp_vec(out, t, x, u, p, vec);
	}


	void gTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->gTfct(out, t, x, p);
	}
	void dgTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dgTdx_vec(out, t, x, p, vec);
	}
	void dgTdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dgTdp_vec(out, t, x, p, vec);
	}
	void dgTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dgTdT_vec(out, t, x, p, vec);
	}


	void hTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->hTfct(out, t, x, p);
	}
	void dhTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhTdx_vec(out, t, x, p, vec);
	}
	void dhTdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhTdp_vec(out, t, x, p, vec);
	}
	void dhTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhTdT_vec(out, t, x, p, vec);
	}


	void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhTdT_vec(out, t, x, u, p);
	}
	void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhTdT_vec(out, t, x, u, p);
	}
	void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->dhTdT_vec(out, t, x, u, p);
	}
	void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *adj, ctypeRNum *p, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->dHdxdt(out, t, x, u, adj, p);
	}
	void Mfct(typeRNum *out, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->Mfct(out);
	}
	void Mtrans(typeRNum *out, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->Mtrans(out);
	}

}
