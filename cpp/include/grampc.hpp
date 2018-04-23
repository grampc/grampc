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
		typeGRAMPC * grampc_;
		ProblemDescription *problem_description_;
	};

}

#endif // GRAMPC_HPP
