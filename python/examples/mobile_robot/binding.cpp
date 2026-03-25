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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "mobile_robot_problem_description.hpp"

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

template <typename Class>
void problem_methods(py::class_<Class, grampc::ProblemDescription, std::shared_ptr<Class>> &cls)
{
    cls.def("ffct", &Class::ffct);
    cls.def("dfdx_vec", &Class::dfdx_vec);
    cls.def("dfdu_vec", &Class::dfdu_vec);
    cls.def("dfdp_vec", &Class::dfdp_vec);

    cls.def("lfct", &Class::lfct);
    cls.def("dldx", &Class::dldx);
    cls.def("dldu", &Class::dldu);
    cls.def("dldp", &Class::dldp);

    cls.def("Vfct", &Class::Vfct);
    cls.def("dVdx", &Class::dVdx);
    cls.def("dVdp", &Class::dVdp);
    cls.def("dVdT", &Class::dVdT);

    cls.def("gfct", &Class::gfct);
    cls.def("dgdx_vec", &Class::dgdx_vec);
    cls.def("dgdu_vec", &Class::dgdu_vec);
    cls.def("dgdp_vec", &Class::dgdp_vec);

    cls.def("hfct", &Class::hfct);
    cls.def("dhdx_vec", &Class::dhdx_vec);
    cls.def("dhdu_vec", &Class::dhdu_vec);
    cls.def("dhdp_vec", &Class::dhdp_vec);

    cls.def("gTfct", &Class::gTfct);
    cls.def("dgTdx_vec", &Class::dgTdx_vec);
    cls.def("dgTdp_vec", &Class::dgTdp_vec);
    cls.def("dgTdT_vec", &Class::dgTdT_vec);

    cls.def("hTfct", &Class::hTfct);
    cls.def("dhTdx_vec", &Class::dhTdx_vec);
    cls.def("dhTdp_vec", &Class::dhTdp_vec);
    cls.def("dhTdT_vec", &Class::dhTdT_vec);
}

PYBIND11_MODULE(grampc_mobile_robot, m)
{
    py::module_::import("pygrampc");

    py::class_<MobileRobotProblemDescription, grampc::ProblemDescription, std::shared_ptr<MobileRobotProblemDescription>> mobile_robot_problem(m, "MobileRobotProblem");
    mobile_robot_problem.def(py::init<const std::vector<typeRNum>&, const std::vector<typeRNum>&, const std::vector<typeRNum>&>());

    problem_methods<MobileRobotProblemDescription>(mobile_robot_problem);
}