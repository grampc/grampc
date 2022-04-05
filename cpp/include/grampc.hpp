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


#ifndef GRAMPC_HPP
#define GRAMPC_HPP

#include "problem_description.hpp"


namespace grampc
{

	/** C++ interface for GRAMPC solver */
	class Grampc
	{
	public:
		/** Create GRAMPC solver for problem description */
		Grampc(ProblemDescription *problem_description);

		/** Free allocated memory */
		~Grampc();

	private:
		Grampc(const Grampc&);

		Grampc& operator=(const Grampc&);

	public:
		/** Set option with float/double value */
		void setopt_real(const char* optName, ctypeRNum optValue);

		/** Set option with int value */
		void setopt_int(const char* optName, ctypeInt optValue);

		/** Set option with string value */
		void setopt_string(const char* optName, const char* optValue);

		/** Set option with float/double vector value */
		void setopt_real_vector(const char* optName, ctypeRNum* optValue);

		/** Set option with int vector value */
		void setopt_int_vector(const char* optName, ctypeInt* optValue);

		/** Print options */
		void printopt() const;

		/** Set parameter with float/double value */
		void setparam_real(const char* paramName, ctypeRNum paramValue);

		/** Set parameter with float/double vector value */
		void setparam_real_vector(const char* paramName, ctypeRNum* paramValue);

		/** Print parameters */
		void printparam() const;

        /** Estimate bound for minimal penalty parameter */
        typeInt estim_penmin(ctypeInt rungrampc);

		/** Run GRAMPC solver */
		void run();

		/** Print solver status */
		typeInt printstatus(ctypeInt status, ctypeInt level);

		/** Get access to options */
		const typeGRAMPCopt* getOptions() const;

		/** Get access to parameters */
		const typeGRAMPCparam* getParameters() const;

		/** Get access to workspace */
		const typeGRAMPCrws* getWorkspace() const;

		/** Get access to solution */
		const typeGRAMPCsol* getSolution() const;

	private:
#ifdef FIXEDSIZE
        typeGRAMPCparam grampc_param_;
        typeGRAMPCopt grampc_opt_;
        typeGRAMPCsol grampc_sol_;
        typeGRAMPCrws grampc_rws_;
        typeGRAMPC grampc_struct_;
#endif
        typeGRAMPC *grampc_;
		ProblemDescription *problem_description_;
	};

}

#endif // GRAMPC_HPP
