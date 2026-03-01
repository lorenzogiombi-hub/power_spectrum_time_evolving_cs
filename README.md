# power_spectrum_time_evolving_cs

This module allows to compute the gravitational wave power spectrum from sound waves according to the formula
```math
\mathcal{P}_{\rm{gw}}=  3(1+\omega)^2\mathcal{H}_*^2 \left(\frac{a_*}{a_{\text{r}}}\right)^{\frac{2\nu}{1+\nu}}\frac{k^3}{2\pi^2} \iint_{\eta_*}^{\eta_{\text{end}}} d\eta_1  d\eta_2 \left(\frac{\eta_*^2}{\eta_1\eta_2}\right)^{1-\nu} G_{k}^\prime(\eta, \eta_1)G_{k}^\prime(\eta, \eta_2) U_{\mathcal{S}}(k, \eta_1, \eta_2),
```
The module operates a 4-dimensional integration with the simpson rule.
