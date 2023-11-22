# GRAMPC
GRAMPC is a nonlinear MPC framework that is suitable for dynamical systems with sampling times in the (sub)millisecond range and that allows for an efficient implementation on embedded hardware. The algorithm is based on an augmented Lagrangian formulation with a tailored gradient method for the inner minimization problem. GRAMPC is implemented in plain C with an additional interface to C++ and MATLAB/Simulink.

The basic structure and usage of GRAMPC are described in the [documentation](doc/manual.pdf) that comes along with the source files. More details about the algorithm and its performance can be found in the corresponding article published in Optimization and Engeneering. The article is available online under open access at: https://doi.org/10.1007/s11081-018-9417-2.

Please cite the paper when you are using results obtained with GRAMPC.

## Support
You can find a discussion board on the [project's sourceforge page](https://sourceforge.net/p/grampc/discussion/general).

## Related projects
GRAMPC is used as the underlying solver in the distributed model predictive control framework [GRAMPC-D](https://github.com/grampc/grampc-d).
