/* This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright (C) 2014-2018 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * Developed at the Institute of Measurement, Control, and Microtechnology,
 * Ulm University. All rights reserved.
 *
 * GRAMPC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * GRAMPC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>
 *
 */



#include "grampc.hpp"


namespace grampc
{


	Grampc::Grampc(ProblemDescription *problem_description)
		: problem_description_(problem_description)
	{
		grampc_init(&grampc_, problem_description_);
	}

	Grampc::~Grampc()
	{
		grampc_free(&grampc_);
	}


	void Grampc::setopt_real(const char *optName, ctypeRNum optValue)
	{
		grampc_setopt_real(grampc_, optName, optValue);
	}

	void Grampc::setopt_int(const char *optName, ctypeInt optValue)
	{
		grampc_setopt_int(grampc_, optName, optValue);
	}

	void Grampc::setopt_string(const char *optName, const char *optValue)
	{
		grampc_setopt_string(grampc_, optName, optValue);
	}

	void Grampc::setopt_real_vector(const char *optName, ctypeRNum *optValue)
	{
		grampc_setopt_real_vector(grampc_, optName, optValue);
	}

	void Grampc::setopt_int_vector(const char *optName, ctypeInt *optValue)
	{
		grampc_setopt_int_vector(grampc_, optName, optValue);
	}

	void Grampc::printopt() const
	{
		grampc_printopt(grampc_);
	}


	void Grampc::setparam_real(const char *paramName, ctypeRNum paramValue)
	{
		grampc_setparam_real(grampc_, paramName, paramValue);
	}

	void Grampc::setparam_real_vector(const char *paramName, ctypeRNum *paramValue)
	{
		grampc_setparam_real_vector(grampc_, paramName, paramValue);
	}

	void Grampc::printparam() const
	{
		grampc_printparam(grampc_);
	}


	void Grampc::run()
	{
		grampc_run(grampc_);
	}

	typeInt Grampc::printstatus(ctypeInt status, ctypeInt level)
	{
		return grampc_printstatus(status, level);
	}

	const typeGRAMPCopt* Grampc::getOptions() const
	{
		return grampc_->opt;
	}

	const typeGRAMPCparam* Grampc::getParameters() const
	{
		return grampc_->param;
	}

	const typeGRAMPCrws* Grampc::getWorkspace() const
	{
		return grampc_->rws;
	}

	const typeGRAMPCsol* Grampc::getSolution() const
	{
		return grampc_->sol;
	}

}
