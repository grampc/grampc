#ifndef FINITE_DIFF_H
#define FINITE_DIFF_H

#include "grampc.h"

/* Enum specifying with respect to which argument the derivatives are computed */
typedef enum {
    DX,
    DU,
    DP,
    DT
} typeFiniteDiffTarget;

/* Function signature used for finite differences */
typedef void(*typeFiniteDiffFctPtr)(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam);

void Vfct_wrapper(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *dummy, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam);
void gTfct_wrapper(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *dummy, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam);
void hTfct_wrapper(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *dummy, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam);

/* Check gradients of probfct by evaluating relative deviation w.r.t. finite differences comparable to
https://de.mathworks.com/help/optim/ug/checkgradients.html#mw_2fc8fe2e-5e32-4697-aa66-b007be316c38 */
void grampc_check_gradients(typeGRAMPC *grampc, ctypeRNum tolerance, ctypeRNum step_size);

void check_tolerance(const char *function, ctypeRNum *fun_fd, ctypeRNum *fun, ctypeRNum tolerance, ctypeInt col, ctypeInt num_rows);

/* Approximate gradient of function using finite differences */
void finite_diff_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam,
                     typeRNum *memory, ctypeRNum step_size, const typeFiniteDiffTarget target, ctypeInt size_vec, const typeFiniteDiffFctPtr func);

/* Approximate dfdx and dfdxtrans using finite differences */
void finite_diff_dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam,
                      typeRNum *memory, ctypeRNum step_size, ctypeInt ml, ctypeInt mu, typeBoolean transpose);

#endif 