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
#include "pygrampc_interface.hpp"
#include <chrono>
#include <stdexcept>
#include <memory>

namespace grampc
{
    GrampcBinding::GrampcBinding(ProblemDescriptionPtr problem)
        : problem_description(problem)
    {
        grampc_init(&grampc_, problem_description.get());
        param.remap_memory(grampc_);
        opt.remap_memory(grampc_);
        rws.remap_memory(grampc_);
        sol.remap_memory(grampc_);
    }

    GrampcBinding::~GrampcBinding()
    {
        grampc_free(&grampc_);
    }

    typeRNum GrampcBinding::run()
    {
        // check if dt and Thor are valid. This redundancy is necessary to throw an error for Python
        if (*param.dt <= 0.0)
        {
            throw std::runtime_error("Sampling time dt is not valid. Must be greater than zero");
        } 
        if (*param.Thor < *param.dt)
        {
            throw std::runtime_error("Horizon Thor is not valid. Must be greater than sampling time");
        }
        if (*opt.Integrator == INT_SYSDISC) {
            if (*param.Thor != (*opt.Nhor - 1) * *param.dt) {
                throw std::runtime_error(DISCRETE_NOT_VALID);
            }
            if (*opt.OptimTime) {
                throw std::runtime_error(OPTIM_TIME_NOT_VALID);
            }
        }

        auto begin = std::chrono::high_resolution_clock::now();
        grampc_run(grampc_);
        auto end = std::chrono::high_resolution_clock::now();
        // time in milliseconds
        return (typeRNum)std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() * 1e-6; 
    }

    void GrampcBinding::estim_penmin(bool run_grampc)
    {
        // check if dt and Thor are valid. This redundancy is necessary to throw an error for Python
        if (*param.dt <= 0.0 && run_grampc)
        {
            throw std::runtime_error("Sampling time dt is not valid. Must be greater than zero");
        } 
        else if (*param.Thor < *param.dt && run_grampc)
        {
            throw std::runtime_error("Horizon Thor is not valid. Must be greater than sampling time");
        }

        grampc_estim_penmin(grampc_, run_grampc);
    }

    void GrampcBinding::check_gradients(ctypeRNum tolerance, ctypeRNum step_size)
    {
        grampc_check_gradients(grampc_, tolerance, step_size);
    }

    void GrampcBinding::get_config_from_file(const std::string& fileName)
    {
        grampc_get_config_from_file(grampc_, fileName.c_str());
    }

    Vector GrampcBinding::ffct(ctypeRNum t, const Vector x, const Vector u, const Vector p, const GrampcParam& param)
    {
        if (x.size() != problem_description->Nx || u.size() != problem_description->Nu || p.size() != problem_description->Np) 
        {
            throw std::runtime_error("Size of x, u or p does not match the problem description");
        }
        Vector out(problem_description->Nx);
        problem_description->ffct(out, t, x, u, p, param);
        return out;
    }

    typeRNum GrampcBinding::lfct(ctypeRNum t, const Vector x, const Vector u, const Vector p, const GrampcParam& param)
    {
        if (x.size() != problem_description->Nx || u.size() != problem_description->Nu || p.size() != problem_description->Np) 
        {
            throw std::runtime_error("Size of x, u or p does not match the problem description");
        }
        Vector out(1);
        problem_description->lfct(out, t, x, u, p, param);
        return out[0];
    }

    typeRNum GrampcBinding::Vfct(ctypeRNum t, const Vector x, const Vector p, const GrampcParam& param)
    {
        if (x.size() != problem_description->Nx || p.size() != problem_description->Np) 
        {
            throw std::runtime_error("Size of x or p does not match the problem description");
        }
        Vector out(1);
        problem_description->Vfct(out, t, x, p, param);
        return out[0];
    }

    Vector GrampcBinding::gfct(ctypeRNum t, const Vector x, const Vector u, const Vector p, const GrampcParam& param)
    {
        if (x.size() != problem_description->Nx || u.size() != problem_description->Nu || p.size() != problem_description->Np) 
        {
            throw std::runtime_error("Size of x, u or p does not match the problem description");
        }
        Vector out(problem_description->Ng);
        problem_description->gfct(out, t, x, u, p, param);
        return out;
    }

    Vector GrampcBinding::hfct(ctypeRNum t, const Vector x, const Vector u, const Vector p, const GrampcParam& param)
    {
        if (x.size() != problem_description->Nx || u.size() != problem_description->Nu || p.size() != problem_description->Np) 
        {
            throw std::runtime_error("Size of x, u or p does not match the problem description");
        }
        Vector out(problem_description->Nh);
        problem_description->hfct(out, t, x, u, p, param);
        return out;
    }

    Vector GrampcBinding::gTfct(ctypeRNum t, const Vector x, const Vector p, const GrampcParam& param)
    {
        if (x.size() != problem_description->Nx || p.size() != problem_description->Np) 
        {
            throw std::runtime_error("Size of x or p does not match the problem description");
        }
        Vector out(problem_description->NgT);
        problem_description->gTfct(out, t, x, p, param);
        return out;
    }

    Vector GrampcBinding::hTfct(ctypeRNum t, const Vector x, const Vector p, const GrampcParam& param)
    {
        if (x.size() != problem_description->Nx || p.size() != problem_description->Np) 
        {
            throw std::runtime_error("Size of x or p does not match the problem description");
        }
        Vector out(problem_description->NhT);
        problem_description->hTfct(out, t, x, p, param);
        return out;
    }

    void GrampcBinding::set_param_real(const std::string &key, typeRNum value)
    {
        grampc_setparam_real(grampc_, key.c_str(), value);
    }

    void GrampcBinding::set_param_real_vec(const std::string &key, const std::vector<typeRNum> &values)
    {
        grampc_setparam_real_vector(grampc_, key.c_str(), values.data());
    }

    void GrampcBinding::set_opt_str(const std::string &key, const std::string &vstr)
    {
        grampc_setopt_string(grampc_, key.c_str(), vstr.c_str());
        if (key == "Integrator" || key == "LineSearchType" || key == "IntegratorCost")
        {    
            rws.remap_memory(grampc_);
        }
    }

    void GrampcBinding::set_opt_int(const std::string &key, int value)
    {
        grampc_setopt_int(grampc_, key.c_str(), value);
        if (key == "Nhor" || key == "MaxMultIter" || key == "MaxGradIter")
        {
            rws.remap_memory(grampc_);
            sol.remap_memory(grampc_);
        }
    }

    void GrampcBinding::set_opt_int_vec(const std::string &key, const std::vector<int> &values)
    {
        grampc_setopt_int_vector(grampc_, key.c_str(), values.data());
        if (values[4] < problem_description->Nx) // Check for MLJAC, if a banded structure is used
        {
            problem_description->Rodas_Jac = problem_description->Nx * (values[4] + values[5] + 1); // Nx * (MLJAC + MUJAC + 1)
        }
        else if (values[4] == problem_description->Nx)
        {
            problem_description->Rodas_Jac = problem_description->Nx * problem_description->Nx;
        }

        if (values[6] < problem_description->Nx) // Check for MLMAS, if a banded structure is used
        {
            problem_description->Rodas_M = problem_description->Nx * (values[6] + values[7] + 1); // Nx * (MLMAS + MUMAS + 1)
        }
        else if (values[6] == problem_description->Nx)
        {
            problem_description->Rodas_M = problem_description->Nx * problem_description->Nx;
        }
    }

    void GrampcBinding::set_opt_real(const std::string &key, typeRNum value)
    {
        grampc_setopt_real(grampc_, key.c_str(), value);
    }

    void GrampcBinding::set_opt_real_vec(const std::string &key, const std::vector<typeRNum> &values)
    {
        grampc_setopt_real_vector(grampc_, key.c_str(), values.data());
    }

    void GrampcBinding::fill_rws_memory(Eigen::Ref<Matrix> rws_matrix, const Eigen::Ref<const Matrix>& new_data)
    {
        if (rws_matrix.data() == nullptr)
        {
            throw std::runtime_error("Memory for matrix is not allocated");
        }
        if (rws_matrix.rows() != new_data.rows() || rws_matrix.cols() != new_data.cols())
        {
            // Formats an error message which outputs the expected dimensions and the dimensions of 
            // new_data. Results in an message like: 
            // "Wrong dimensions detected. Expected (2, 10), got (2, 15)."
            std::ostringstream error_message; 
            error_message << "Wrong dimensions detected. Expected (" << rws_matrix.rows() << ", " << rws_matrix.cols() << "), got ("
                << new_data.rows() << ", " << new_data.cols() << ")";
            throw std::length_error(error_message.str());
        }

        // copy the new data to the rws_matrix.
        rws_matrix = new_data;
    }

    void GrampcBinding::set_rws_u(const Eigen::Ref<const Matrix>& u_new)
    {
        fill_rws_memory(rws.u, u_new);
    }

    void GrampcBinding::set_rws_multiplier(const Eigen::Ref<const Matrix>& multiplier_new)
    {
        fill_rws_memory(rws.mult, multiplier_new);
    }

    void GrampcBinding::set_rws_penalty(const Eigen::Ref<const Matrix>& penalty_new)
    {
        fill_rws_memory(rws.pen, penalty_new);
    }

    void GrampcBinding::print_params()
    {
        grampc_printparam(grampc_);
    }

    void GrampcBinding::print_opts()
    {
        grampc_printopt(grampc_);
    }

    void GrampcBinding::print_status()
    {
        grampc_printstatus(sol.status, STATUS_LEVEL_DEBUG);
    }
}