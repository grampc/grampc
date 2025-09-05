.. _sec:AlgOpt:LineSearch:

Line search
-----------

The projected gradient method as part of the GRAMPC algorithm
:prf:ref:`alg:AlgOpt:GrampcAlgorithm`
requires the solution of the line search
problem :footcite:`Englert:OE:2019`

.. math::
    :label: eq:AlgOpt:LineSearchProblem

    \min_{\alpha} \bar{J} \left( 
        \mb{\psi}_{\mb{u}} \left( \mb{u}^{i|j} - \alpha \mb{d}_{\mb{u}}^{i|j} \right),
        \mb{\psi}_{\mb{p}} \left( \mb{p}^{i|j} - \gamma_{\mb{p}} \alpha \mb{d}_{\mb{p}}^{i|j} \right),
        \psi_{T} \left( T^{i|j} - \gamma_{T} \alpha d_{T}^{i|j} \right);
        \mb{\bar \mu}, \mb{\bar c}, \mb{x}_0
    \right)\,.

GRAMPC implements two efficient strategies to solve this problem in an
approximate manner. The following options apply to both line search
methods that are detailed below
(:ref:`sec:AlgOpt:LineSearchAdaptive` and
:ref:`sec:AlgOpt:LineSearchExplicit`):

-  ``LineSearchType``: This option selects either the adaptive line search strategy (value
   ``adaptive``) or the explicit approach (value ``explicit1`` or ``explicit2``).

-  ``LineSearchExpAutoFallback``: If this option is activated, the automatic fallback strategy is
   used in the case that the explicit formulas result in negative step sizes.

-  ``LineSearchMax``: This option sets the maximum value :math:`\alpha_{\max}` of the
   step size :math:`\alpha`.

-  ``LineSearchMin``: This option sets the minimum value :math:`\alpha_{\min}` of the
   step size :math:`\alpha`.

-  ``LineSearchInit``: Indicates the initial value :math:`\alpha_{\text{init}}>0` for the
   step size :math:`\alpha`. If the adaptive line search is used, the
   sample point :math:`\alpha_2` is set
   to :math:`\alpha_2 = \alpha_{\text{init}}`.

-  ``OptimParamLineSearchFactor``: This option sets the adaptation factor
   :math:`\gamma_{\mb{p}}` that weights the update of the
   parameter vector :math:`\mb{p}` against the update of the
   control :math:`\mb{u}`.

-  ``OptimTimeLineSearchFactor``: This option sets the adaptation factor :math:`\gamma_{T}` that
   weights the update of the end time :math:`T` against the update of
   the control :math:`\mb{u}`.

.. _sec:AlgOpt:LineSearchAdaptive:

Adaptive line search
~~~~~~~~~~~~~~~~~~~~

An appropriate way to determine the step size is the adaptive line
search approach from :footcite:`InTech_GraichenKaepernick2012`,
where a polynomial approximation of the cost :math:`\bar J` is used and an adaptation of
the search intervals is performed. More precisely, the cost functional
is evaluated at three sample points
:math:`\alpha_1 < \alpha_2 < \alpha_3` with
:math:`\alpha_2=\frac{1}{2}\left(\alpha_1+\alpha_3\right)`, which are
used to construct a quadratic polynomial of the cost according to

.. math::
    :label: eq:AlgOpt:ls_adapt_approx

    \bar{J}\left(
        \mb{\psi}_{\mb{u}} \left(\mb{u} ^{i|j} - \alpha \mb{d}_{\mb{u}} ^{i|j}\right),
        \mb{\psi}_{\mb{p}} \left(\mb{p} ^{i|j} - \gamma_{\mb{p}}\alpha \mb{d}_{\mb{p}} ^{i|j}\right),
        \psi_{T} \left(T ^{i|j} - \gamma_{T}\alpha   d_{T} ^{i|j} \right);
        \mb{\bar \mu}, \mb{\bar c}, \mb{x}_0
    \right)\\
    \approx \Phi(\alpha) = p_0 + p_1 \alpha + p_2 \alpha^2 \,.

.. figure:: ../img/ls_adapt.*
    :name: fig:AlgOpt:LS_Adapt

    Adaptation of line search interval.

Subsequently, a step size :math:`\alpha^{j}` is computed by minimizing
the cost approximation :math:numref:`eq:AlgOpt:ls_adapt_approx`. If
necessary, the interval :math:`[\alpha_1,\alpha_3]` is adapted for the
next gradient iteration in the following way

.. math::
    :label: eq:AlgOpt:IntervalAdaptation

    [\alpha_1,\alpha_3] & \leftarrow 
    \begin{cases}
    \hfill \kappa\,[\alpha_1,\alpha_3] & \text{if} \ \alpha \geq
    \alpha_3 - \varepsilon_{\alpha}(\alpha_3-\alpha_1) \ \text{and} \
    \alpha_3 \leq \alpha_{\max} \ \text{and} \
    |\Phi(\alpha_1)-\Phi(\alpha_3)| >  \varepsilon_{\phi}
    \\[1.5ex]
    \hfill \frac{1}{\kappa} \,[\alpha_1,\alpha_3] & \text{if} \ \alpha \leq
    \alpha_1 + \varepsilon_{\alpha}(\alpha_3-\alpha_1) \ \text{and} \
    \alpha_1 \geq \alpha_{\min}\ \text{and} \
    |\Phi(\alpha_1)-\Phi(\alpha_3)| > \varepsilon_{\phi}
    \\[1.5ex]
    \hfill[\alpha_1,\alpha_3] & \text{otherwise}
    \end{cases} \\[1.5ex]
    %
    \alpha_2 & \leftarrow \frac{1}{2}\left(\alpha_1+\alpha_3\right)

with the adaptation factor :math:`\kappa > 1`, the interval tolerance
:math:`\varepsilon_{\alpha} \in (0,0.5)`, the absolute cost tolerance
:math:`{\varepsilon_{\phi}\in[0,\infty)}` for adapting the interval and
the interval bounds :math:`\alpha_{\max}>\alpha_{\min}>0`.
The modification :math:numref:`eq:AlgOpt:IntervalAdaptation` of
the line search interval tracks the minimum point :math:`\alpha^{j}` of
the line search problem in the case when :math:`\alpha^{j}` is either
outside of the interval :math:`[\alpha_1,\alpha_3]` or close to one of
the outer bounds :math:`\alpha_1`, :math:`\alpha_3`, as illustrated in
:numref:`fig:AlgOpt:LS_Adapt`. The adaptation factor
:math:`\kappa` accounts for scaling as well as shifting of the interval
:math:`[\alpha_1,\alpha_3]` in the next gradient iteration, if
:math:`\alpha^{j}` lies in the vicinity of the interval bounds
:math:`[\alpha_1,\alpha_3]` as specified by the interval tolerance
:math:`\varepsilon_\alpha`. This adaptive strategy allows one to track
the minimum of the line search problem :math:numref:`eq:AlgOpt:LineSearchProblem` over
the gradient iterations :math:`j` and MPC steps :math:`k`, while
guaranteeing a fixed number of operations in view of a real-time MPC
implementation. The absolute tolerance :math:`\varepsilon_{\phi}` of the
difference in the (scaled) costs at the interval bounds
:math:`|\Phi(\alpha_1)-\Phi(\alpha_3)|` avoids oscillations of the
interval width in regions where the cost function :math:`\bar{J}` is
almost constant.

The following options apply specifically to the adaptive line search
strategy:

-  ``LineSearchAdaptAbsTol``: This option sets the absolute tolerance :math:`\varepsilon_{\phi}`
   of the difference in costs at the interval bounds :math:`\alpha_1`
   and :math:`\alpha_2`. If the difference in the (scaled) costs on
   these bounds falls below :math:`\varepsilon_{\phi}`, the adaption of
   the interval is stopped in order to avoid oscillations.

-  ``LineSearchAdaptFactor``: This option sets the adaptation factor :math:`\kappa > 1` in
   :math:numref:`eq:AlgOpt:IntervalAdaptation`
   that determines how much the line search interval can be adapted from
   one gradient iteration to the next.

-  ``LineSearchIntervalTol``: This option sets the interval tolerance
   :math:`\varepsilon_{\alpha} \in (0,0.5)` in
   :math:numref:`eq:AlgOpt:IntervalAdaptation`
   that determines for which values of :math:`\alpha` the adaption is
   performed.

-  ``LineSearchIntervalFactor``: This option sets the interval factor :math:`\beta \in (0,1)` that
   specifies the interval bounds :math:`[\alpha_1,\alpha_3]` according
   to :math:`\alpha_1 = \alpha_2 (1 - \beta)` and
   :math:`\alpha_3 = \alpha_2 (1 + \beta)`, whereby the mid sample point
   is initialized as :math:`\alpha_2 = \alpha_\text{init}`.

.. _sec:AlgOpt:LineSearchExplicit:

Explicit line search
~~~~~~~~~~~~~~~~~~~~

An alternative way to determine the step size in order to further reduce
the computational effort for time-critical problems is the explicit line
search approach originally discussed in :footcite:`Barzilai1988`
and adapted in :footcite:`Kaepernick2013` for the optimal
control case. The motivation is to minimize the difference between two
consecutive control updates :math:`u_k^{i|j}(\tau)` and
:math:`u_k^{i|j+1}(\tau)` in the unconstrained case and additionally
assuming the same step size :math:`\alpha^{i|j}`, i.e.

.. math::
    :label: eq:AlgOpt:ls_expl_prob

    \alpha^{i|j} &= \underset{\alpha > 0}{\arg\min} \; \Bigl\| {\mb{u}}^{i|j+1} - {\mb{u}}^{i|j} \Bigr\|^2_{L^2_m[0,T]} + \Bigl\| {\mb{p}}^{i|j+1} - {\mb{p}}^{i|j} \Bigr\|^2_2 + \Bigl| {T}^{i|j+1} - {T}^{i|j} \Bigr| \\
    &= \underset{\alpha > 0}{\arg\min} \; \Bigl\| \underbrace{{\mb{u}}^{i|j} - {\mb{u}}^{j-1}}_{=:\Delta {\mb{u}}^{i|j}}- \alpha \underbrace{\left(\mb{d}^{i|j}_{\mb{u}}-\mb{d}^{j-1}_{\mb{u}}\right)}_{=:\Delta \mb{d}^{i|j}_{\mb{u}}} \Bigr\|^2_{L^2_m[0,T]} + \Bigl\| \underbrace{{\mb{p}}^{i|j} - {\mb{p}}^{j-1}}_{=:\Delta {\mb{p}}^{i|j}}- \gamma_{\mb{p}} \alpha \underbrace{\left(\mb{d}^{i|j}_{\mb{p}}-\mb{d}^{j-1}_{\mb{p}}\right)}_{=:\Delta \mb{d}^{i|j}_{\mb{p}}} \Bigr\|^2_2 \\
    & \qquad\qquad +\Bigl\| \underbrace{{T}^{i|j} - {T}^{j-1}}_{=:\Delta {T}^{i|j}}- \gamma_{T} \alpha \underbrace{\left({d}^{i|j}_{T}-{d}^{j-1}_{T}\right)}_{=:\Delta {d}^{i|j}_{T}} \Bigr\|^2_2 

with :math:`\| z \|^2_{L^2_m[0,T]} = \langle z,z \rangle :=\int_0^T z^\mathsf{T}(t) z(t) \, {\rm d}t`.

.. figure:: ../img/ls_explicit.*
    :name: fig:AlgOpt:LS_Explicit

    Motivation for the explicit line search strategy.

:numref:`fig:AlgOpt:LS_Explicit` illustrates the general idea
behind :math:numref:`eq:AlgOpt:ls_expl_prob`. To solve
:math:numref:`eq:AlgOpt:ls_expl_prob`, consider the
following function

.. math::
    :label: eq:AlgOpt:q_alpha

    q(\alpha) : & = 
    \Bigl\| \Delta u_k^{j} - \alpha \Delta d_k^{j} \Bigr\|^2_{L^2_m[0,T]}
    \\
    & = \int_{0}^{T} \left( \Delta u_k^{j} - \alpha \Delta d_k^{j} \right)^\mathsf{T}
    \left( \Delta u_k^{j} - \alpha \Delta d_k^{j} \right) \, {\rm d}t 
    \\
    & = \int_{0}^{T} \left(\Delta u_k^{j}\right)^\mathsf{T}\Delta u_k^{j} \, {\rm d}t
    + \alpha^2 \int_{0}^{T} \left(\Delta d_k^{j}\right)^\mathsf{T}\Delta d_k^{j} \, {\rm d}t 
    - 2 \alpha \int_{0}^{T} \left(\Delta u_k^{j}\right)^\mathsf{T}\Delta d_k^{j} \, {\rm d}t.

The minimum has to satisfy the stationarity condition

.. math::

    \frac{\partial^{} q(\alpha)}{\partial \alpha^{}} = 
    2 \alpha \int_{0}^{T} \left(\Delta d_k^{j}\right)^\mathsf{T}\Delta d_k^{j} \, {\rm d}t 
    - 2 \int_{0}^{T} \left(\Delta u_k^{j} \right)^\mathsf{T}\Delta d_k^{j} \, {\rm d}t 
    = 0 .

A suitable step size :math:`\alpha^{j}` then follows to

.. math::

    \alpha^{j} = 
    \frac{\int_{0}^{T} \left( \Delta u_k^{j}\right)^\mathsf{T}\Delta d_k^{j} \, {\rm d}t}
    {\int_{0}^{T} \left( \Delta d_k^{j}\right)^\mathsf{T}\Delta d_k^{j} \, {\rm d}t} 
    = \frac{\langle \Delta u_k^{j},\Delta d_k^{j}\rangle}
    {\langle \Delta d_k^{j},\Delta d_k^{j} \rangle} \,.

Another way to compute an appropriate step size for the control update
can be achieved by reformulating
:math:numref:`eq:AlgOpt:q_alpha` in the following way:

.. math::

   \label{eq:AlgOpt:OptStepSizeQbar}
   q(\alpha) = 
   \Bigl\| \Delta u_k^{j} - \alpha \Delta d_k^{j} \Bigr\|^2_{L^2_m[0,T]} = 
   \alpha^2 \Bigl\| \frac{1}{\alpha} \Delta u_k^{j} - \Delta d_k^{j} \Bigr\|^2_{L^2_m[0,T]} 
   =: \alpha^2 \bar q(\alpha).

In the subsequent, the new function :math:`\bar q(\alpha)` is minimized
w.r.t. the step size leading to a similar solution

.. math::

   \label{eq:AlgOpt:OptStepSize2}
   \alpha^{j} = \frac{\langle\Delta u^{j}, \Delta u^{j} \rangle}
   {\langle \Delta u^{j}, \Delta d_k^{j} \rangle} \,.

For the original problem :math:numref:`eq:AlgOpt:ls_expl_prob`, the solution

.. math::
    :label: eq:AlgOpt:ls_expl1

    \alpha^{i|j} = 
    \frac{\langle \Delta \mb{u}^{i|j},\Delta\mb{d}^{i|j}_{\mb{u}}\rangle 
        +\gamma_{\mb{p}} \langle \Delta\mb{p}^{i|j},\Delta\mb{d}^{i|j}_{\mb{p}}\rangle 
        +\gamma_{T} \Delta T^{i|j}\Delta{d}^{i|j}_{T}}
    {\langle \Delta\mb{d}^{i|j}_{\mb{u}},\Delta \mb{d}^{i|j}_{\mb{u}} \rangle 
        + \gamma_{\mb{p}}^2 \langle{\Delta\mb{d}^{i|j}_{\mb{p}}}^\mathsf{T}\Delta \mb{d}^{i|j}_{\mb{p}}\rangle 
        + \gamma_{T}^2\big( \Delta {d}^{i|j}_{T}\big)^2}

.. math::
    :label: eq:AlgOpt:ls_expl2
    
    \alpha^{i|j} = 
    \frac{\langle \Delta \mb{u}^{i|j}, \Delta\mb{u}^{i|j}\rangle 
        + \gamma_{\mb{p}} \langle\Delta\mb{p}^{i|j}, \Delta\mb{p}^{i|j} \rangle
        + \gamma_{T} \big(\Delta T^{i|j}\big)^2}
    {\langle \Delta\mb{u}^{i|j}, \Delta \mb{d}^{i|j}_{\mb{u}} \rangle 
        + \gamma_{\mb{p}}^2 \langle\Delta\mb{p}^{i|j}, \Delta \mb{d}^{i|j}_{\mb{p}} \rangle 
        + \gamma_{T}^2 \Delta {T}^{i|j} \Delta {d}^{i|j}_{T}} 

follows.

In the GRAMPC implementation, both approaches :math:numref:`eq:AlgOpt:ls_expl1` and
:math:numref:`eq:AlgOpt:ls_expl2` are available. 
In addition, the step size :math:`\alpha^{j}` is bounded by the upper and
lower values :math:`\alpha_{\max} > \alpha_{\min} > 0`. However, if the
originally computed step size :math:`\alpha^{j}` is less than zero [1]_,
either the initial step size :math:`\alpha^{j} = \alpha_\text{init}` or
the automatic fallback strategy that is detailed in the next subsection
is used in order to achieve a valid step size. The fallback strategy is
set with the following option (only available for the explicit line
search strategies):

-  ``LineSearchExpAutoFallback``: If this option is activated, the automatic fallback strategy is
   used in the case that the explicit formulas result in negative step sizes.

.. _sec:AlgOpt:LineSearchFallback:

Fallback strategy for explicit line search
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While the initial step size can be used as fallback solution if the
explicit step size computation yields negative values, it often requires
problem-specific tuning of :math:`\alpha_\text{init}` for achieving
optimal performance. As alternative, GRAMPC implements an automatic
fallback strategy that is based on the idea of using at most 1% of the
control range defined by ``umax`` and ``umin``. For this purpose, the maximum absolute
value
:math:`\mb{d}^{i|j}_{\mb{u},\text{max}} = \| \mb{d}_{\mb{u}} ^{i|j}(t) \|_{L^\infty}`
of the search direction
:math:`\mb{d}_{\mb{u}} ^{i|j}(t)` over the horizon is
determined. Subsequently, the step size

.. math::

   	\alpha^{i|j} = \frac{1}{100} \cdot \min_{k\in\{1,\dots,N_{\mb{u}}\}}\left\lbrace  
    \frac{u_{\text{max},k}- u_{\text{min},k}} {d_{\mb{u},\text{max},k}^{i|j}}  \right\rbrace

follows as the minimal step size required to perform a step of 1% with
respect to the range of at least one control in at least one time step.
Additionally, the step size is limited to 10% of the maximum step size
:math:`\alpha_{\max}`. Since this strategy requires reasonable limits
for the controls, it is only executed if ``LineSearchExpAutoFallback`` is activated and if these
limits are defined by the user. Furthermore, this strategy can only be
used if ``OptimControl`` is switched ``on``. In all other case, the initial step size
:math:`\alpha^{j}=\alpha_\text{init}` will be used as fallback solution.

.. rubric:: Footnotes

.. [1] It can be shown that the step size :math:`\alpha` is negative if the cost function is locally non-convex.

.. footbibliography::