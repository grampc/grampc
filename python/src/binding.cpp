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
#include <pygrampc_interface.hpp>
#include <pygrampc_problem_description.hpp>
#include <pygrampc_types.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

using namespace grampc;

PYBIND11_MODULE(_core, m)
{
    pybind11::class_<GrampcParam>(m, "GrampcParam")
        .def(pybind11::init<>())
        .def_readonly("Nx", &GrampcParam::Nx)
        .def_readonly("Nu", &GrampcParam::Nu)
        .def_readonly("Np", &GrampcParam::Np)
        .def_readonly("Ng", &GrampcParam::Ng)
        .def_readonly("Nh", &GrampcParam::Nh)
        .def_readonly("NgT", &GrampcParam::NgT)
        .def_readonly("NhT", &GrampcParam::NhT)
        .def_readonly("Nc", &GrampcParam::Nc)

        .def_readonly("x0", &GrampcParam::x0)
        .def_readonly("xdes", &GrampcParam::xdes)

        .def_readonly("u0", &GrampcParam::u0)
        .def_readonly("udes", &GrampcParam::udes)
        .def_readonly("umax", &GrampcParam::umax)
        .def_readonly("umin", &GrampcParam::umin)

        .def_readonly("p0", &GrampcParam::p0)
        .def_readonly("pmax", &GrampcParam::pmax)
        .def_readonly("pmin", &GrampcParam::pmin)

        .def_readonly("Thor", &GrampcParam::Thor)
        .def_readonly("Tmax", &GrampcParam::Tmax)
        .def_readonly("Tmin", &GrampcParam::Tmin)
        .def_readonly("dt", &GrampcParam::dt)
        .def_readonly("t0", &GrampcParam::t0);
    
    pybind11::class_<GrampcOpt>(m, "GrampcOpt")
        .def(pybind11::init<>())
        .def_readonly("Nhor", &GrampcOpt::Nhor)
        .def_readonly("MaxGradIter", &GrampcOpt::MaxGradIter)
        .def_readonly("MaxMultIter", &GrampcOpt::MaxMultIter)
        .def_readonly("ShiftControl", &GrampcOpt::ShiftControl)

        .def_readonly("TimeDiscretization", &GrampcOpt::TimeDiscretization)

        .def_readonly("IntegralCost", &GrampcOpt::IntegralCost)
        .def_readonly("TerminalCost", &GrampcOpt::TerminalCost)
        .def_readonly("IntegratorCost", &GrampcOpt::IntegratorCost)

        .def_readonly("Integrator", &GrampcOpt::Integrator)
        .def_readonly("IntegratorRelTol", &GrampcOpt::IntegratorRelTol)
        .def_readonly("IntegratorAbsTol", &GrampcOpt::IntegratorAbsTol)
        .def_readonly("IntegratorMinStepSize", &GrampcOpt::IntegratorMinStepSize)
        .def_readonly("IntegratorMaxSteps", &GrampcOpt::IntegratorMaxSteps)
        .def_readonly("FlagsRodas", &GrampcOpt::FlagsRodas)

        .def_readonly("LineSearchType", &GrampcOpt::LineSearchType)
        .def_readonly("LineSearchExpAutoFallback", &GrampcOpt::LineSearchExpAutoFallback)
        .def_readonly("LineSearchMax", &GrampcOpt::LineSearchMax)
        .def_readonly("LineSearchMin", &GrampcOpt::LineSearchMin)
        .def_readonly("LineSearchInit", &GrampcOpt::LineSearchInit)
        .def_readonly("LineSearchAdaptAbsTol", &GrampcOpt::LineSearchAdaptAbsTol)
        .def_readonly("LineSearchAdaptFactor", &GrampcOpt::LineSearchAdaptFactor)
        .def_readonly("LineSearchIntervalTol", &GrampcOpt::LineSearchIntervalTol)
        .def_readonly("LineSearchIntervalFactor", &GrampcOpt::LineSearchIntervalFactor)

        .def_readonly("OptimControl", &GrampcOpt::OptimControl)
        .def_readonly("OptimParam", &GrampcOpt::OptimParam)
        .def_readonly("OptimParamLineSearchFactor", &GrampcOpt::OptimParamLineSearchFactor)
        .def_readonly("OptimTime", &GrampcOpt::OptimTime)
        .def_readonly("OptimTimeLineSearchFactor", &GrampcOpt::OptimTimeLineSearchFactor)

        .def_readonly("ScaleProblem", &GrampcOpt::ScaleProblem)
        .def_readonly("xScale", &GrampcOpt::xScale)
        .def_readonly("xOffset", &GrampcOpt::xOffset)
        .def_readonly("uScale", &GrampcOpt::uScale)
        .def_readonly("uOffset", &GrampcOpt::uOffset)
        .def_readonly("pScale", &GrampcOpt::pScale)
        .def_readonly("pOffset", &GrampcOpt::pOffset)
        .def_readonly("TScale", &GrampcOpt::TScale)
        .def_readonly("TOffset", &GrampcOpt::TOffset)
        .def_readonly("JScale", &GrampcOpt::JScale)
        .def_readonly("cScale", &GrampcOpt::cScale)

        .def_readonly("EqualityConstraints", &GrampcOpt::EqualityConstraints)
        .def_readonly("InequalityConstraints", &GrampcOpt::InequalityConstraints)
        .def_readonly("TerminalEqualityConstraints", &GrampcOpt::TerminalEqualityConstraints)
        .def_readonly("TerminalInequalityConstraints", &GrampcOpt::TerminalInequalityConstraints)
        .def_readonly("ConstraintsHandling", &GrampcOpt::ConstraintsHandling)
        .def_readonly("ConstraintsAbsTol", &GrampcOpt::ConstraintsAbsTol)

        .def_readonly("MultiplierMax", &GrampcOpt::MultiplierMax)
        .def_readonly("MultiplierDampingFactor", &GrampcOpt::MultiplierDampingFactor)
        .def_readonly("PenaltyMax", &GrampcOpt::PenaltyMax)
        .def_readonly("PenaltyMin", &GrampcOpt::PenaltyMin)
        .def_readonly("PenaltyIncreaseFactor", &GrampcOpt::PenaltyIncreaseFactor)
        .def_readonly("PenaltyDecreaseFactor", &GrampcOpt::PenaltyDecreaseFactor)
        .def_readonly("PenaltyIncreaseThreshold", &GrampcOpt::PenaltyIncreaseThreshold)
        .def_readonly("AugLagUpdateGradientRelTol", &GrampcOpt::AugLagUpdateGradientRelTol)

        .def_readonly("ConvergenceCheck", &GrampcOpt::ConvergenceCheck)
        .def_readonly("ConvergenceGradientRelTol", &GrampcOpt::ConvergenceGradientRelTol);

    pybind11::class_<GrampcSol>(m, "GrampcSol")
        .def(pybind11::init<>())
        .def_readonly("xnext", &GrampcSol::xnext)
        .def_readonly("unext", &GrampcSol::unext)
        .def_readonly("pnext", &GrampcSol::pnext)
        .def_readonly("Tnext", &GrampcSol::Tnext)
        .def_readonly("J", &GrampcSol::J)
        .def_readonly("cfct", &GrampcSol::cfct)
        .def_readonly("pen", &GrampcSol::pen)
        .def_readonly("iter", &GrampcSol::iter)
        .def_readonly("status", &GrampcSol::status);

    pybind11::class_<GrampcRWS>(m, "GrampcRWS")
        .def(pybind11::init<>())
        .def_readonly("t", &GrampcRWS::t)
        .def_readonly("tls", &GrampcRWS::tls)

        .def_readonly("x", &GrampcRWS::x)
        .def_readonly("adj", &GrampcRWS::adj)
        .def_readonly("dcdx", &GrampcRWS::dcdx)

        .def_readonly("u", &GrampcRWS::u)
        .def_readonly("uls", &GrampcRWS::uls)
        .def_readonly("uprev", &GrampcRWS::uprev)
        .def_readonly("gradu", &GrampcRWS::gradu)
        .def_readonly("graduprev", &GrampcRWS::graduprev)
        .def_readonly("dcdu", &GrampcRWS::dcdu)

        .def_readonly("p", &GrampcRWS::p)
        .def_readonly("pls", &GrampcRWS::pls)
        .def_readonly("pprev", &GrampcRWS::pprev)
        .def_readonly("gradp", &GrampcRWS::gradp)
        .def_readonly("gradpprev", &GrampcRWS::gradpprev)
        .def_readonly("dcdp", &GrampcRWS::dcdp)

        .def_readonly("T", &GrampcRWS::T)
        .def_readonly("Tprev", &GrampcRWS::Tprev)
        .def_readonly("gradT", &GrampcRWS::gradT)
        .def_readonly("dcdt", &GrampcRWS::dcdt)

        .def_readonly("mult", &GrampcRWS::mult)
        .def_readonly("pen", &GrampcRWS::pen)
        .def_readonly("cfct", &GrampcRWS::cfct)
        .def_readonly("cfctprev", &GrampcRWS::cfctprev)
        .def_readonly("cfctAbsTol", &GrampcRWS::cfctAbsTol)

        .def_readonly("lsAdapt", &GrampcRWS::lsAdapt)
        .def_readonly("lsExplicit", &GrampcRWS::lsExplicit)
        .def_readonly("rwsScale", &GrampcRWS::rwsScale)
        .def_readonly("lrwsGeneral", &GrampcRWS::lrwsGeneral)
        .def_readonly("rwsGeneral", &GrampcRWS::rwsGeneral)

        .def_readonly("lworkRodas", &GrampcRWS::lworkRodas)
        .def_readonly("liworkRodas", &GrampcRWS::liworkRodas)
        .def_readonly("rparRodas", &GrampcRWS::rparRodas)
        .def_readonly("iparRodas", &GrampcRWS::iparRodas)
        .def_readonly("workRodas", &GrampcRWS::workRodas)
        .def_readonly("iworkRodas", &GrampcRWS::iworkRodas);

    pybind11::class_<GrampcBinding> (m, "GrampcBinding")
        .def(pybind11::init<ProblemDescriptionPtr>())
        .def_readonly("param", &GrampcBinding::param)
        .def_readonly("opt", &GrampcBinding::opt)
        .def_readonly("sol", &GrampcBinding::sol)
        .def_readonly("rws", &GrampcBinding::rws)
        .def_readwrite("problem", &GrampcBinding::problem_description)
        .def("run", &GrampcBinding::run)
        .def("_set_param_real", &GrampcBinding::set_param_real)
        .def("_set_param_real_vec", &GrampcBinding::set_param_real_vec)
        .def("_set_opt_str", &GrampcBinding::set_opt_str)
        .def("_set_opt_int", &GrampcBinding::set_opt_int)
        .def("_set_opt_int_vec", &GrampcBinding::set_opt_int_vec)
        .def("_set_opt_real", &GrampcBinding::set_opt_real)
        .def("_set_opt_real_vec", &GrampcBinding::set_opt_real_vec)
        .def("set_rws_u", &GrampcBinding::set_rws_u)
        .def("set_rws_multiplier", &GrampcBinding::set_rws_multiplier)
        .def("set_rws_penalty", &GrampcBinding::set_rws_penalty)
        .def("print_opts", &GrampcBinding::print_opts)
        .def("print_params", &GrampcBinding::print_params)
        .def("print_status", &GrampcBinding::print_status)

        .def("ffct", &GrampcBinding::ffct)
        .def("lfct", &GrampcBinding::lfct)
        .def("Vfct", &GrampcBinding::Vfct)
        .def("gfct", &GrampcBinding::gfct)
        .def("hfct", &GrampcBinding::hfct)
        .def("gTfct", &GrampcBinding::gTfct)
        .def("hTfct", &GrampcBinding::hTfct)
        .def("get_config_from_file", &GrampcBinding::get_config_from_file)
        .def("estim_penmin", &GrampcBinding::estim_penmin)
        .def("check_gradients", &GrampcBinding::check_gradients);

    pybind11::class_<PyProblemDescription, PyProblem, ProblemDescriptionPtr>(m, "ProblemDescription")
        .def(pybind11::init<ctypeInt,ctypeInt,ctypeInt,ctypeInt,ctypeInt,ctypeInt,ctypeInt>(), 
            pybind11::arg("Nx"), pybind11::arg("Nu"), pybind11::arg("Np"), pybind11::arg("Ng"), pybind11::arg("Nh"), pybind11::arg("NgT"), pybind11::arg("NhT"))
        .def_readonly("Nx", &PyProblemDescription::Nx)
        .def_readonly("Nu", &PyProblemDescription::Nu)
        .def_readonly("Np", &PyProblemDescription::Np)
        .def_readonly("Ng", &PyProblemDescription::Ng)
        .def_readonly("Nh", &PyProblemDescription::Nh)
        .def_readonly("NgT", &PyProblemDescription::NgT)
        .def_readonly("NhT", &PyProblemDescription::NhT)

        .def("ffct", &PyProblemDescription::ffct)
        .def("dfdx_vec", &PyProblemDescription::dfdx_vec)
        .def("dfdu_vec", &PyProblemDescription::dfdu_vec)
        .def("dfdp_vec", &PyProblemDescription::dfdp_vec)

        .def("lfct", &PyProblemDescription::lfct)
        .def("dldx", &PyProblemDescription::dldx)
        .def("dldu", &PyProblemDescription::dldu)
        .def("dldp", &PyProblemDescription::dldp)

        .def("Vfct", &PyProblemDescription::Vfct)
        .def("dVdx", &PyProblemDescription::dVdx)
        .def("dVdp", &PyProblemDescription::dVdp)
        .def("dVdT", &PyProblemDescription::dVdT)

        .def("gfct", &PyProblemDescription::gfct)
        .def("dgdx_vec", &PyProblemDescription::dgdx_vec)
        .def("dgdu_vec", &PyProblemDescription::dgdu_vec)
        .def("dgdp_vec", &PyProblemDescription::dgdp_vec)

        .def("hfct", &PyProblemDescription::hfct)
        .def("dhdx_vec", &PyProblemDescription::dhdx_vec)
        .def("dhdu_vec", &PyProblemDescription::dhdu_vec)
        .def("dhdp_vec", &PyProblemDescription::dhdp_vec)

        .def("gTfct", &PyProblemDescription::gTfct)
        .def("dgTdx_vec", &PyProblemDescription::dgTdx_vec)
        .def("dgTdp_vec", &PyProblemDescription::dgTdp_vec)
        .def("dgTdT_vec", &PyProblemDescription::dgTdT_vec)

        .def("hTfct", &PyProblemDescription::hTfct)
        .def("dhTdx_vec", &PyProblemDescription::dhTdx_vec)
        .def("dhTdp_vec", &PyProblemDescription::dhTdp_vec)
        .def("dhTdT_vec", &PyProblemDescription::dhTdT_vec)

        .def("dfdx", &PyProblemDescription::dfdx)
        .def("dfdxtrans", &PyProblemDescription::dfdxtrans)
        .def("dfdt", &PyProblemDescription::dfdt)
        .def("dHdxdt", &PyProblemDescription::dHdxdt)
        .def("Mfct", &PyProblemDescription::Mfct)
        .def("Mtrans", &PyProblemDescription::Mtrans);
}