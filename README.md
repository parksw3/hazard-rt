## Theory

There are broadly two ways of modeling changes in infection incidence. First, we can assume that epidemic dynamics are driven by constant-strength interventions (basically what everyone else does): 
$$
i(t) = R_0 P(t) \int_0^\infty g(\tau) i(t-\tau) d\tau.
$$
In this case, we can estimate the instantaneous reproduction number $R(t)$ using the intrinsic generation-interval distribution $g(\tau)$ which does not change over time:
$$
R(t) = R_0 P(t) = \frac{i(t)}{\int_0^\infty g(\tau) i(t-\tau) d\tau}.
$$
Instead, we can also assume that epidemic dynamics are driven by constant-speed interventions:
$$
i(t) = R_0 \int_0^\infty\left[ \exp\left(-\int_0^\tau h(t-\sigma) d\sigma \right) g(\tau) i(t-\tau) \right] d\tau,
$$
where $h(t)$ represents the hazard of isolation at time $t$. For example, an individual infected at time $t-\tau$ will experience the hazard from $h(t-\tau)$ to $h(t)$ at time $t$ and are therefore less likely to transmit. In this case, calculating the instantaneous reproduction number $R(t)$ is similar to before:
$$
\begin{align}
R(t) &= R_0 \int_0^\infty\left[ \exp\left(-\int_0^\tau h(t-\sigma) d\sigma \right) g(\tau) \right] d\tau\\
&=  \frac{i(t)}{\int_0^\infty g_t(\tau) i(t-\tau) d\tau}
\end{align}
$$
But we need a time-varying instantaneous generation-interval distribution:
$$
g_t(\tau) = \frac{ \exp\left(-\int_0^\tau h(t-\sigma) d\sigma \right) g(\tau)}{\int_0^\infty\left[ \exp\left(-\int_0^\tau h(t-\sigma) d\sigma \right) g(\tau) \right] d\tau}.
$$
More generally, given a kernel $K(t, \tau)$ and an associated intervention function $I(t, \tau)$, we have
$$
\begin{aligned}
K(t, \tau) &= R_0 I(t, \tau) g(\tau) \\
i(t) &= \int_0^\infty K(t, \tau) i(t-\tau) d\tau\\
R(t) &= \int_0^\infty K(t, \tau) d\tau\\
g(t, \tau) &= K(t,\tau)/R(t)\\
R(t) &= \frac{i(t)}{\int_0^\infty g(t,\tau) i(t-\tau) d\tau}
\end{aligned}
$$
Constant-strength $I(t, \tau)  = P(t)$ and constant-speed $I(t, \tau) = \exp\left(-\int_0^\tau h(t-\sigma) d\sigma \right)$ are just special cases.

## Application

First, some discretization steps:
$$
i(t) = \sum_{\tau=1}^t \exp\left(-\sum_{\sigma=1}^\tau h(t-\sigma) \right) g(\tau) i(t-\tau).
$$

