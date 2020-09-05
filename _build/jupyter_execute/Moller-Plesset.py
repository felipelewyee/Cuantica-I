# Moller-Plesset

Las ecuaciones de Moller-Plesset surgen de tomar un sistema representado por Hartree-Fock como sistema conocido ($\mathcal{H}_0$) y aplicarle una perturbación ($\mathcal{V}$) para convertirlo en el sistema con correlación ($\mathcal{H}$).
\begin{equation}
\mathcal{H} = \mathcal{H}_0 + \mathcal{V}
\end{equation}

Recordando de teoría de perturbaciones, es posible realizar correcciones de n-orden a la energía $E_i^{(n)}$ y a la función de onda $\psi_i^{(n)}$. La corrección a la energía del estado basal ($i=0$) toma la forma

\begin{equation}
E_0 = E_0^{(0)} + E_0^{(1)} + E_0^{(2)} + E_0^{(3)} + ...
\end{equation}

y en Moller-Plesset estos términos son:

- $E_0^{(0)} = \sum_a {\varepsilon_a}$
- $E_0^{(1)} = -\frac{1}{2} \sum_{ab} [aa|bb] - [ab|ba]$
- $E_0^{(2)} = 2 \sum_{abrs}^{N/2} \frac{[ar|bs][ra|sb]}{\varepsilon_{a}+\varepsilon_{b}-\varepsilon_{r}-\varepsilon_{s}} - \sum_{abrs}^{N/2} \frac{[ar|bs][rb|sa]}{\varepsilon_{a}+\varepsilon_{b}-\varepsilon_{r}-\varepsilon_{s}}$

donde los índices $a$, $b$ refieren a orbitales moleculares ocupados, y $r$, $s$ a orbitales moleculares desocupados, los términos $\varepsilon_a$, $\varepsilon_b$, $\varepsilon_r$, $\varepsilon_s$ refieren a sus energías orbitales, y los términos entre paréntesis son integrales de repulsión electrónica moleculares, que se calculan mediante 
\begin{equation}
[pq|rs] = \sum_{\mu\nu\sigma\lambda}^N C_{\mu p} C_{\nu q} C_{\sigma r} C_{\lambda s} [\mu \nu | \sigma \lambda]
\end{equation}

donde $\mu$, $\nu$, $\sigma$ y $\lambda$ se refieren a orbitales atómicos, y
\begin{equation}
[\mu \nu | \sigma \lambda] = \int \int \frac{\mu(r_1) \nu(r_1) \sigma(r_2) \lambda(r_2)}{ |r_1 - r_2| } dr_1 dr_2
\end{equation}

En química cuántica es de especial interés la corrección a la energía, por lo que dependiendo de su orden de corrección la metodología se denomina MPn.

La energía de Hartree-Fock está dada por
\begin{equation}
E_{HF} = E_0^{(0)} + E_0^{(1)}
\end{equation}

por lo que la correlación electrónica es
\begin{equation}
E_{corr} = E_0^{(2)} + E_0^{(3)} + ...
\end{equation}

