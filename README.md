# UROP

MATLAB code developed to conduct the Undergreaduate research project in Computational Dynamical systems at Imperial College London,
under the supervision of Dr. P. Ramsden. (Summer 2017). Project attempted to, and successfully found evidence for the existence of non-
chaotic, strange attractors in the quasi-periodically driven Duffing oscillator. 

License: MIT

## Usage ##
<a href="https://imgur.com/8l69u7w"><img src="https://i.imgur.com/8l69u7w.jpg" title="source: imgur.com" /></a>

The code essentially revolves around solving the quais-periodically driven Duffing oscillator (shown above) using Runge-Kutta Fourth 
numerical integration. Many of the functions are dedicated to the generation and analysis of Poincare sections made from the "snapshots"
of the flow. Analysis involves multiple methods of finding the fractal dimension of the attractors (Poincare sections) generated using 
three methods mentioned below, as well as Fourier transforms and Lyapunov exponents.

After checking out the wikipedia page on Duffing oscillators and perhaps reading one or two papers (almost any) you should have a good 
idea of what parameters to input into the functions to get the rough kind of solution you want. Using Plot_solution, you can quite easily and quickly see what the actual wave and phase plot look like as you change parameters:

<a href="https://imgur.com/sFb7oRy"><img src="https://i.imgur.com/sFb7oRy.jpg" title="source: imgur.com" /></a>

This example is with the parameters: -0.2, 1, -1.63, 0.3, 0.3, (sqrt(5) - 1)/2, [1, 0], which results in a chaotic system (huge changes
in end result with small changes in initial condition). You can then test this by changing the inital condtions from [1, 0] to perhaps 
[1.0001, 0] and see if the solution is wildly different. If so, you've just found a chaotic system! Feel free to mess around with
parameters, but be aware that some regions of parameter space may result in systems that die off immediately or break the program with
infinite values resulting. Please read up on the physical meanings of the variable to avoid this.

If you want to create your very own beautiful attractor sets, then use Create_attractor with the exact same parameters (with the exception of the initial conditon array being replaced by number of cycles). Here's the kind of amazing sets you can create:

<a href="https://imgur.com/3DCGdyc"><img src="https://i.imgur.com/3DCGdyc.jpg" title="source: imgur.com" /></a>

This was created with the same parameters, with the last changed to number of driving cycles (time in multiples of *2pi*) for the
driving force: -0.2, 1, -1.63, 0.3, 0.3, (sqrt(5) - 1) / 2, 10^5.

Additionally, you can take Fourier transforms to get frequency spectrums, with similar parameter inputs


