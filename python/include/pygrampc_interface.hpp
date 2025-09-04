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
#ifndef PYGRAMPC_INTERFACE_HPP
#define PYGRAMPC_INTERFACE_HPP

extern "C" 
{
    #include "grampc.h"
}
#include "pygrampc_problem_description.hpp"
#include "pygrampc_types.hpp"
#include <string>
#include <vector>

namespace grampc
{
    class GrampcBinding
    {
        public:
            GrampcParam param;
            GrampcOpt opt;
            GrampcSol sol;
            GrampcRWS rws;

            ProblemDescriptionPtr problem_description;

        public:
            /** Create Python interface to the GRAMPC solver for problem description */
            GrampcBinding(ProblemDescriptionPtr problem);
            ~GrampcBinding();

            /** Calls grampc_run and returns the wall clock time of one call in milliseconds. */
            typeRNum run();

            /** Estimate bound for minimal penalty parameter */
            void estim_penmin(bool run_grampc);
            
            /** Check gradients of the supplied problem description */
            void check_gradients(ctypeRNum tolerance, ctypeRNum step_size);

            /** Get options and parameter from file **/
            void get_config_from_file(const std::string& fileName);
            

            /** Set parameter with float/double value */
            void set_param_real(const std::string& key, typeRNum value);

            /** Set parameter with float/double vector value */
            void set_param_real_vec(const std::string& key, const std::vector<typeRNum>& values);

		    
            /** Set option with float/double value */
            void set_opt_real(const std::string& key, typeRNum value);

            /** Set option with string value */
            void set_opt_str(const std::string& key, const std::string& vstr);

            /** Set option with int value */
            void set_opt_int(const std::string& key, int value);

            /** Set option with float/double vector value */
            void set_opt_real_vec(const std::string& key, const std::vector<typeRNum>& values);

            /** Set option with int vector value */
            void set_opt_int_vec(const std::string& key, const std::vector<int>& values);

            void set_rws_u(const Eigen::Ref<const Matrix>& u_new);
            void set_rws_multiplier(const Eigen::Ref<const Matrix>& multiplier_new);
            void set_rws_penalty(const Eigen::Ref<const Matrix>& penalty_new);

            void fill_rws_memory(const Eigen::Ref<Matrix> rws_matrix, const Eigen::Ref<const Matrix>& new_data);

            /** access the probfct like the Matlab interface */
            Vector ffct(ctypeRNum t, const Vector x, const Vector u, const Vector p, const GrampcParam& param);
            typeRNum lfct(ctypeRNum t, const Vector x, const Vector u, const Vector p, const GrampcParam& param);
            typeRNum Vfct(ctypeRNum t, const Vector x, const Vector p, const GrampcParam& param);
            Vector gfct(ctypeRNum t, const Vector x, const Vector u, const Vector p, const GrampcParam& param);
            Vector hfct(ctypeRNum t, const Vector x, const Vector u, const Vector p, const GrampcParam& param);
            Vector gTfct(ctypeRNum t, const Vector x, const Vector p, const GrampcParam& param);
            Vector hTfct(ctypeRNum t, const Vector x, const Vector p, const GrampcParam& param);

            /** Print parameters */
            void print_params();

            /** Print options */
            void print_opts();

            /** Print solver status */
            void print_status();

        private: 
            typeGRAMPC* grampc_;

            /** Copy and assignment constructors are private, to prevent duplicate ressource management of the grampc struct */
            GrampcBinding(const GrampcBinding&);
            GrampcBinding& operator=(const GrampcBinding&);
    };
}

#endif