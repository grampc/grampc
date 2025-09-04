The purpose of this example is to compare the continuous-time formulation to the discrete-time formulation.
As dynamics the example of a double integrator is chosen.
Calling

    [vec, grampc] = startMPC(1, false, 2);

compiles the toolbox as well as both problem variants and runs the continuous-time formulation.
Afterwards it suffices to call

    [vec, grampc] = startMPC(1, false);

and

    [vec, grampc] = startMPC(11, true);
    
to run the continuous-time and discrete-time MPC, respectively.
By using different figure numbers, the results can be compared.
In addition, the combination of model predictive control with moving horizon estimation can be compared using

    [vec, grampcMPC, grampcMHE] = startMHE(1, false, 2)

and similar arguments for figure number, continuous or discrete and the compilation.
