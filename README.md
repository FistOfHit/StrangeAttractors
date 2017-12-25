# UROP

MATLAB code developed to conduct the Undergreaduate research project in Computational Dynamical systems at Imperial College London,
under the supervision of Dr. P. Ramsden. (Summer 2017). Project attempted to, and successfully found evidence for the existence of non-
chaotic, strange attractors in the quasi-periodically driven Duffing oscillator. 

License: MIT

## Usage ##
<img src="https://latex.codecogs.com/gif.latex?\begin{align*} 
y &= \dot{x} \\
\dot{y} &= \alpha y + \beta x + \delta x^3 + \gamma_1 \cos (t) + \gamma_2 \cos (\tau)
\end{align*}" /> 

The code essentially revolves around solving the quais-periodically driven Duffing oscillator (shown above) using Runge-Kutta Fourth 
numerical integration. Many of the functions are dedicated to the generation and analysis of Poincare sections made from the "snapshots"
of the flow. Analysis involves multiple methods of finding the fractal dimension of the attractors (Poincare sections) generated using 
three methods mentioned below, as well as Fourier transforms and Lyapunov exponents.



