#include "finite_diff.h"

void Vfct_wrapper(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *dummy, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    Vfct(out, T, x, p, param, userparam);
}

void gTfct_wrapper(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *dummy, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    gTfct(out, T, x, p, param, userparam);
}

void hTfct_wrapper(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *dummy, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    hTfct(out, T, x, p, param, userparam);
}

void grampc_check_gradients(typeGRAMPC *grampc, ctypeRNum tolerance, ctypeRNum step_size)
{
    typeInt i, num_rows, max_dim, size;
    typeInt ml = grampc->opt->FlagsRodas[4];
    typeInt mu = grampc->opt->FlagsRodas[5];
    typeRNum *tmp = 0;
    typeRNum *fun, *fun_fd, *vec, *memory;

    // maximum of dimensions max(Nx, Nu, Np, Ng, NgT, Nh, NhT)
    max_dim = MAX(grampc->param->Nx, MAX(grampc->param->Nu, MAX(grampc->param->Np, MAX(grampc->param->Ng, MAX(grampc->param->NgT, MAX(grampc->param->Nh, grampc->param->NhT))))));
    // dimension of dfdx is Nx * num_rows depending on banded or full structure
    if (ml != grampc->param->Nx || mu != grampc->param->Nx) {
        num_rows = ml + 1 + mu;
    }
    else {
        num_rows = grampc->param->Nx;
    }
    size = MAX(max_dim, grampc->param->Nx * num_rows);

    tmp = (typeRNum*) malloc((2*size + 4*max_dim) * sizeof(typeRNum));
    if (tmp == NULL)
    {
        printf("Allocation of temporary arrays failed in gradient checker\n");
        return;
    }
    fun = tmp; // size elements
    fun_fd = fun + size; // size elements
    vec = fun_fd + size; // max_dim elements
    memory = vec + max_dim; // 3*max_dim elements

    MatSetScalar(vec, 0.0, max_dim, 1);

    /* Check gradient of ffct */
    // dfdx_vec
    for (i = 0; i < grampc->param->Nx; i++) {
        vec[i] = 1.0;
        MatSetScalar(fun, 0.0, size, 1);
        MatSetScalar(fun_fd, 0.0, size, 1);
        MatSetScalar(memory, 0.0, 3*max_dim, 1);
        finite_diff_vec(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DX, grampc->param->Nx, &ffct);
        dfdx_vec(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam);
        check_tolerance("dfdx_vec", fun_fd, fun, tolerance, i, grampc->param->Nx);
        vec[i] = 0.0;
    }

    // dfdu_vec
    if (grampc->param->Nu > 0) {
        for (typeInt i = 0; i < grampc->param->Nx; i++) {
            vec[i] = 1.0;
            MatSetScalar(fun, 0.0, size, 1);
            MatSetScalar(fun_fd, 0.0, size, 1);
            MatSetScalar(memory, 0.0, 3*max_dim, 1);
            finite_diff_vec(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DU, grampc->param->Nx, &ffct);
            dfdu_vec(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam);
            check_tolerance("dfdu_vec", fun_fd, fun, tolerance, i, grampc->param->Nu);
            vec[i] = 0.0;
        }
    }

    // dfdp_vec
    if (grampc->param->Np > 0) {
        for (i = 0; i < grampc->param->Nx; i++) {
            vec[i] = 1.0;
            MatSetScalar(fun, 0.0, size, 1);
            MatSetScalar(fun_fd, 0.0, size, 1);
            MatSetScalar(memory, 0.0, 3*max_dim, 1);
            finite_diff_vec(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DP, grampc->param->Nx, &ffct);
            dfdp_vec(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam);
            check_tolerance("dfdp_vec", fun_fd, fun, tolerance, i, grampc->param->Np);
            vec[i] = 0.0;
        }
    }

    /* Check gradient of lfct */
    // dldx
    vec[0] = 1.0;
    MatSetScalar(fun, 0.0, size, 1);
    MatSetScalar(fun_fd, 0.0, size, 1);
    MatSetScalar(memory, 0.0, 3*max_dim, 1);
    finite_diff_vec(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DX, 1, &lfct);
    dldx(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, grampc->param, grampc->userparam);
    check_tolerance("dldx", fun_fd, fun, tolerance, 0, grampc->param->Nx);
    vec[0] = 0.0;

    // dldu
    if (grampc->param->Nu > 0) {
        vec[0] = 1.0;
        MatSetScalar(fun, 0.0, size, 1);
        MatSetScalar(fun_fd, 0.0, size, 1);
        MatSetScalar(memory, 0.0, 3*max_dim, 1);
        finite_diff_vec(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DU, 1, &lfct);
        dldu(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, grampc->param, grampc->userparam);
        check_tolerance("dldu", fun_fd, fun, tolerance, 0, grampc->param->Nu);
        vec[0] = 0.0;
    }

    // dldp
    if (grampc->param->Np > 0) {
        vec[0] = 1.0;
        MatSetScalar(fun, 0.0, size, 1);
        MatSetScalar(fun_fd, 0.0, size, 1);
        MatSetScalar(memory, 0.0, 3*max_dim, 1);
        finite_diff_vec(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DP, 1, &lfct);
        dldp(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, grampc->param, grampc->userparam);
        check_tolerance("dldp", fun_fd, fun, tolerance, 0, grampc->param->Np);
        vec[0] = 0.0;
    }

    /* Check gradient of Vfct */
    // dVdx
    vec[0] = 1.0;
    MatSetScalar(fun, 0.0, size, 1);
    MatSetScalar(fun_fd, 0.0, size, 1);
    MatSetScalar(memory, 0.0, 3*max_dim, 1);
    finite_diff_vec(fun_fd, grampc->param->Thor, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DX, 1, &Vfct_wrapper);
    dVdx(fun, grampc->param->Thor, grampc->param->x0, grampc->param->p0, grampc->param, grampc->userparam);
    check_tolerance("dVdx", fun_fd, fun, tolerance, 0, grampc->param->Nx);
    vec[0] = 0.0;

    // dVdp
    if (grampc->param->Np > 0) {
        vec[0] = 1.0;
        MatSetScalar(fun, 0.0, size, 1);
        MatSetScalar(fun_fd, 0.0, size, 1);
        MatSetScalar(memory, 0.0, 3*max_dim, 1);
        finite_diff_vec(fun_fd, grampc->param->Thor, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DP, 1, &Vfct_wrapper);
        dVdp(fun, grampc->param->Thor, grampc->param->x0, grampc->param->p0, grampc->param, grampc->userparam);
        check_tolerance("dVdp", fun_fd, fun, tolerance, 0, grampc->param->Np);
        vec[0] = 0.0;
    }

    // dVdT
    vec[0] = 1.0;
    MatSetScalar(fun, 0.0, size, 1);
    MatSetScalar(fun_fd, 0.0, size, 1);
    MatSetScalar(memory, 0.0, 3*max_dim, 1);
    finite_diff_vec(fun_fd, grampc->param->Thor, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DT, 1, &Vfct_wrapper);
    dVdT(fun, grampc->param->Thor, grampc->param->x0, grampc->param->p0, grampc->param, grampc->userparam);
    check_tolerance("dVdT", fun_fd, fun, tolerance, 0, 1);
    vec[0] = 0.0;

    /* Check gradient of gfct */
    if (grampc->param->Ng > 0) {
        // dgdx_vec
        for (i = 0; i < grampc->param->Ng; i++) {
            vec[i] = 1.0;
            MatSetScalar(fun, 0.0, size, 1);
            MatSetScalar(fun_fd, 0.0, size, 1);
            MatSetScalar(memory, 0.0, 3*max_dim, 1);
            finite_diff_vec(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DX, grampc->param->Ng, &gfct);
            dgdx_vec(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam);
            check_tolerance("dgdx_vec", fun_fd, fun, tolerance, i, grampc->param->Nx);
            vec[i] = 0.0;
        }

        // dgdu_vec
        if (grampc->param->Nu > 0) {
            for (i = 0; i < grampc->param->Ng; i++) {
                vec[i] = 1.0;
                MatSetScalar(fun, 0.0, size, 1);
                MatSetScalar(fun_fd, 0.0, size, 1);
                MatSetScalar(memory, 0.0, 3*max_dim, 1);
                finite_diff_vec(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DU, grampc->param->Ng, &gfct);
                dgdu_vec(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam);
                check_tolerance("dgdu_vec", fun_fd, fun, tolerance, i, grampc->param->Nu);
                vec[i] = 0.0;
            }
        }

        // dgdp_vec
        if (grampc->param->Np > 0) {
            for (i = 0; i < grampc->param->Ng; i++) {
                vec[i] = 1.0;
                MatSetScalar(fun, 0.0, size, 1);
                MatSetScalar(fun_fd, 0.0, size, 1);
                MatSetScalar(memory, 0.0, 3*max_dim, 1);
                finite_diff_vec(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DP, grampc->param->Ng, &gfct);
                dgdp_vec(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam);
                check_tolerance("dgdp_vec", fun_fd, fun, tolerance, i, grampc->param->Np);
                vec[i] = 0.0;
            }
        }
    }

    /* Check gradient of hfct */
    if (grampc->param->Nh > 0) {
        // dhdx_vec
        for (i = 0; i < grampc->param->Nh; i++) {
            vec[i] = 1.0;
            MatSetScalar(fun, 0.0, size, 1);
            MatSetScalar(fun_fd, 0.0, size, 1);
            MatSetScalar(memory, 0.0, 3*max_dim, 1);
            finite_diff_vec(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DX, grampc->param->Nh, &hfct);
            dhdx_vec(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam);
            check_tolerance("dhdx_vec", fun_fd, fun, tolerance, i, grampc->param->Nx);
            vec[i] = 0.0;
        }

        // dhdu_vec
        if (grampc->param->Nu > 0) {
            for (i = 0; i < grampc->param->Nh; i++) {
                vec[i] = 1.0;
                MatSetScalar(fun, 0.0, size, 1);
                MatSetScalar(fun_fd, 0.0, size, 1);
                MatSetScalar(memory, 0.0, 3*max_dim, 1);
                finite_diff_vec(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DU, grampc->param->Nh, &hfct);
                dhdu_vec(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam);
                check_tolerance("dhdu_vec", fun_fd, fun, tolerance, i, grampc->param->Nu);
                vec[i] = 0.0;
            }
        }

        // dhdp_vec
        if (grampc->param->Np > 0) {
            for (i = 0; i < grampc->param->Nh; i++) {
                vec[i] = 1.0;
                MatSetScalar(fun, 0.0, size, 1);
                MatSetScalar(fun_fd, 0.0, size, 1);
                MatSetScalar(memory, 0.0, 3*max_dim, 1);
                finite_diff_vec(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DP, grampc->param->Nh, &hfct);
                dhdp_vec(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam);
                check_tolerance("dhdp_vec", fun_fd, fun, tolerance, i, grampc->param->Np);
                vec[i] = 0.0;
            }
        }
    }

    /* Check gradient of gTfct */
    if (grampc->param->NgT > 0) {
        // dgTdx_vec
        for (i = 0; i < grampc->param->NgT; i++) {
            vec[i] = 1.0;
            MatSetScalar(fun, 0.0, size, 1);
            MatSetScalar(fun_fd, 0.0, size, 1);
            MatSetScalar(memory, 0.0, 3*max_dim, 1);
            finite_diff_vec(fun_fd, grampc->param->Thor, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DX, grampc->param->NgT, &gTfct_wrapper);
            dgTdx_vec(fun, grampc->param->Thor, grampc->param->x0, grampc->param->p0, vec, grampc->param, grampc->userparam);
            check_tolerance("dgTdx_vec", fun_fd, fun, tolerance, i, grampc->param->Nx);
            vec[i] = 0.0;
        }

        // dgTdp_vec
        if (grampc->param->Np > 0) {
            for (i = 0; i < grampc->param->NgT; i++) {
                vec[i] = 1.0;
                MatSetScalar(fun, 0.0, size, 1);
                MatSetScalar(fun_fd, 0.0, size, 1);
                MatSetScalar(memory, 0.0, 3*max_dim, 1);
                finite_diff_vec(fun_fd, grampc->param->Thor, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DP, grampc->param->NgT, &gTfct_wrapper);
                dgTdp_vec(fun, grampc->param->Thor, grampc->param->x0, grampc->param->p0, vec, grampc->param, grampc->userparam);
                check_tolerance("dgTdp_vec", fun_fd, fun, tolerance, i, grampc->param->Np);
                vec[i] = 0.0;
            }
        }

        // dgTdT_vec
        vec[0] = 1.0;
        MatSetScalar(fun, 0.0, size, 1);
        MatSetScalar(fun_fd, 0.0, size, 1);
        MatSetScalar(memory, 0.0, 3*max_dim, 1);
        finite_diff_vec(fun_fd, grampc->param->Thor, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DT, grampc->param->NgT, &gTfct_wrapper);
        dgTdT_vec(fun, grampc->param->Thor, grampc->param->x0, grampc->param->p0, vec, grampc->param, grampc->userparam);
        check_tolerance("dgTdT_vec", fun_fd, fun, tolerance, 0, 1);
        vec[0] = 0.0;
    }

    /* Check gradient of hTfct */
    if (grampc->param->NhT > 0) {
        // dhTdx_vec
        for (i = 0; i < grampc->param->NhT; i++) {
            vec[i] = 1.0;
            MatSetScalar(fun, 0.0, size, 1);
            MatSetScalar(fun_fd, 0.0, size, 1);
            MatSetScalar(memory, 0.0, 3*max_dim, 1);
            finite_diff_vec(fun_fd, grampc->param->Thor, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DX, grampc->param->NhT, &hTfct_wrapper);
            dhTdx_vec(fun, grampc->param->Thor, grampc->param->x0, grampc->param->p0, vec, grampc->param, grampc->userparam);
            check_tolerance("dhTdx_vec", fun_fd, fun, tolerance, i, grampc->param->Nx);
            vec[i] = 0.0;
        }

        // dhTdp_vec
        if (grampc->param->Np > 0) {
            for (i = 0; i < grampc->param->NhT; i++) {
                vec[i] = 1.0;
                MatSetScalar(fun, 0.0, size, 1);
                MatSetScalar(fun_fd, 0.0, size, 1);
                MatSetScalar(memory, 0.0, 3*max_dim, 1);
                finite_diff_vec(fun_fd, grampc->param->Thor, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DP, grampc->param->NhT, &hTfct_wrapper);
                dhTdp_vec(fun, grampc->param->Thor, grampc->param->x0, grampc->param->p0, vec, grampc->param, grampc->userparam);
                check_tolerance("dhTdp_vec", fun_fd, fun, tolerance, i, grampc->param->Np);
                vec[i] = 0.0;
            }
        }

        // dhTdT_vec
        vec[0] = 1.0;
        MatSetScalar(fun, 0.0, size, 1);
        MatSetScalar(fun_fd, 0.0, size, 1);
        MatSetScalar(memory, 0.0, 3*max_dim, 1);
        finite_diff_vec(fun_fd, grampc->param->Thor, grampc->param->x0, grampc->param->u0, grampc->param->p0, vec, grampc->param, grampc->userparam, memory, step_size, DT, grampc->param->NhT, &hTfct_wrapper);
        dhTdT_vec(fun, grampc->param->Thor, grampc->param->x0, grampc->param->p0, vec, grampc->param, grampc->userparam);
        check_tolerance("dhTdT_vec", fun_fd, fun, tolerance, 0, 1);
        vec[0] = 0.0;
    }

    /* Check gradient of ffct for implicit integrators */
    if (grampc->opt->Integrator == INT_RODAS) {
        // dfdx
        MatSetScalar(fun, 0.0, size, 1);
        MatSetScalar(fun_fd, 0.0, size, 1);
        MatSetScalar(memory, 0.0, 3*max_dim, 1);
        finite_diff_dfdx(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, grampc->param, grampc->userparam,
                         memory, step_size, ml, mu, 0);
        dfdx(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, grampc->param, grampc->userparam);
        for(i = 0; i < grampc->param->Nx; ++i) {
            check_tolerance("dfdx", fun_fd + i * num_rows, fun + i * num_rows, tolerance, i, num_rows);
        }

        // dfdxtrans
        MatSetScalar(fun, 0.0, size, 1);
        MatSetScalar(fun_fd, 0.0, size, 1);
        MatSetScalar(memory, 0.0, 3*max_dim, 1);
        finite_diff_dfdx(fun_fd, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, grampc->param, grampc->userparam,
                         memory, step_size, ml, mu, 1);
        dfdxtrans(fun, grampc->param->t0, grampc->param->x0, grampc->param->u0, grampc->param->p0, grampc->param, grampc->userparam);
        for(i = 0; i < grampc->param->Nx; ++i) {
            check_tolerance("dfdxtrans", fun_fd + i * num_rows, fun + i * num_rows, tolerance, i, num_rows);
        }
    }

    free(tmp);
}

/* Checks the relative deviation of the gradient like in
https://de.mathworks.com/help/optim/ug/checkgradients.html#mw_2fc8fe2e-5e32-4697-aa66-b007be316c38 */
void check_tolerance(const char *function, ctypeRNum *fun_fd, ctypeRNum *fun, ctypeRNum tolerance, ctypeInt col, ctypeInt num_rows)
{
    typeRNum rel_diff;
    for (typeInt row = 0; row < num_rows; row++) {
        rel_diff = ABS((fun_fd[row] - fun[row]) / MAX(1.0, ABS(fun[row])));
        if (rel_diff > tolerance) {
            printf("%s: element (%u,%u) exceeds supplied tolerance with %e\n", function, row, col, rel_diff);
        }
    }
}

void finite_diff_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam,
                     typeRNum *memory, ctypeRNum step_size, const typeFiniteDiffTarget target, ctypeInt func_out_size, const typeFiniteDiffFctPtr func)
{
    typeInt i, j, size;
    typeRNum *fun = memory; // size 2*func_out_size + size(=Nx/Nu/Np/1)
    typeRNum *delta_fun = fun + func_out_size;
    typeRNum *delta_arg = delta_fun + func_out_size;

    /* select the derivative target x, u, p or T */
    switch(target) {
    case DX:
        size = param->Nx;
        MatCopy(delta_arg, x, size, 1);
        break;
    case DU:
        size = param->Nu;
        MatCopy(delta_arg, u, size, 1);
        break;
    case DP:
        size = param->Np;
        MatCopy(delta_arg, p, size, 1);
        break;
    case DT:
        size = 1;
        *delta_arg = t;
        break;
    default:
        size = 0; // no-op for this function as the main loop does not execute
        break;
    }

    (*func)(fun, t, x, u, p, param, userparam);
    for(i = 0; i < size; i++) {
        delta_arg[i] += step_size;
        switch(target) {
        case DX:
            (*func)(delta_fun, t, delta_arg, u, p, param, userparam);
            break;
        case DU:
            (*func)(delta_fun, t, x, delta_arg, p, param, userparam);
            break;
        case DP:
            (*func)(delta_fun, t, x, u, delta_arg, param, userparam);
            break;
        case DT:
            (*func)(delta_fun, delta_arg[0], x, u, p, param, userparam);
            break;                                                                                  
        }
        delta_arg[i] -= step_size;

        out[i] = 0.0;
        for(j = 0; j < func_out_size; j++) {
            out[i] += vec[j] * (delta_fun[j] - fun[j]) / step_size;
        }
    }
}

void finite_diff_dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam,
                      typeRNum *memory, ctypeRNum step_size, ctypeInt ml, ctypeInt mu, typeBoolean transpose)
{
    typeInt i, jmin, jmax;
    typeInt m = ml + 1 + mu;
    typeInt size = param->Nx;
    typeRNum *fun = memory; // size 3*Nx
    typeRNum *delta_fun = fun + size;
    typeRNum *delta_arg = delta_fun + size;

    MatCopy(delta_arg, x, size, 1);
    ffct(fun, t, x, u, p, param, userparam);
    for(i = 0; i < size; i++) {
        delta_arg[i] += step_size;
        ffct(delta_fun, t, delta_arg, u, p, param, userparam);
        delta_arg[i] -= step_size;

        if (ml != param->Nx || mu != param->Nx) {
            // banded matrix
            if (transpose) {
                // dfdxtrans
                jmin = MAX(0, i - ml); // shifts the start index by the correct number of zeros at the beginning
                jmax = MIN(size, i + 1 + mu); // shifts the end index by the correct number of zeros at the end
                for(typeInt j = jmin; j < jmax; j++) {
                    out[j * m + (i - j + ml)] = (delta_fun[j] - fun[j]) / step_size;
                }
            }
            else {
                // dfdx
                jmin = MAX(0, i - mu); // shifts the start index by the correct number of zeros at the beginning
                jmax = MIN(size, i + 1 + ml); // shifts the end index by the correct number of zeros at the end
                for(typeInt j = jmin; j < jmax; j++) {
                    out[i * m + (j - i + mu)] = (delta_fun[j] - fun[j]) / step_size;
                }
            }
        }
        else {
            // full matrix
            if (transpose) {
                // dfdxtrans
                for(typeInt j = 0; j < param->Nx; j++) {
                    out[j * param->Nx + i] = (delta_fun[j] - fun[j]) / step_size;
                }
            }
            else {
                // dfdx
                for(typeInt j = 0; j < param->Nx; j++) {
                    out[i * param->Nx + j] = (delta_fun[j] - fun[j]) / step_size;
                }
            }
        }
    }
}
