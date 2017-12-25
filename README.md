# UROP #

MATLAB code developed to conduct the Undergreaduate research project in Computational Dynamical systems at Imperial College London,
under the supervision of Dr. P. Ramsden. (Summer 2017). Project attempted to, and successfully found evidence for the existence of non-
chaotic, strange attractors in the quasi-periodically driven Duffing oscillator. 

## Summary of abilities ##
- Solves Quasi-periodically driven Duffing oscillator (can adapt to almost any similar ODE really)
- Find frequency spectra of solutions to the oscillator
- Determines if a system with a certain set of parameters is chaotic or not via Lyapunov exponents
- Creates (estimates of) attractors by taking Poincare sections in 2 and 3 dimensions (2 spacial and 1 temporal)
- Calculates fractal dimensions of such attractors with 3 distinct and useful methods.

## Requirements ##

- Matlab 2017a or newer (older versions may work, check their documentation on functions such as "unique", and especially on logical
arrays)

## Description and Usage ##
<a href="https://imgur.com/8l69u7w"><img src="https://i.imgur.com/8l69u7w.jpg" title="source: imgur.com" /></a>

The code essentially revolves around solving the Quasi-periodically driven Duffing oscillator (shown above) using Runge-Kutta Fourth
order numerical integration. Many of the functions are dedicated to the generation and analysis of Poincare sections made from the
"snapshots" of the flow. Analysis involves multiple methods of finding the fractal dimension of the attractors (Poincare sections)
generated using three methods mentioned below, as well as Fourier transforms and Lyapunov exponents.

After checking out the wikipedia page on Duffing oscillators and perhaps reading one or two papers (almost any) you should have a good 
idea of what parameters to input into the functions to get the rough kind of solution you want. Using Plot_solution, you can quite
easily and quickly see what the actual wave and phase plot look like as you change parameters:

<a href="https://imgur.com/sFb7oRy"><img src="https://i.imgur.com/sFb7oRy.jpg" title="source: imgur.com" /></a>

This example is with the parameters: -0.2, 1, -1.63, 0.3, 0.3, (sqrt(5) - 1)/2, [1, 0], which results in a chaotic system (huge changes
in end result with small changes in initial condition). You can then test this by changing the inital condtions from [1, 0] to perhaps 
[1.0001, 0] and see if the solution is wildly different. If so, you've just found a chaotic system! Feel free to mess around with
parameters, but be aware that some regions of parameter space may result in systems that die off immediately or break the program with
infinite values resulting. Please read up on the physical meanings of the variables to avoid this.

If you want to create your very own beautiful attractor sets, then use Create_attractor with the exact same parameters (with the
exception of the initial conditon array being replaced by number of cycles). Here's the kind of amazing sets you can create:

<a href="https://imgur.com/3DCGdyc"><img src="https://i.imgur.com/3DCGdyc.jpg" title="source: imgur.com" /></a>

This was created with the same parameters, with the last changed to number of driving cycles (time in multiples of *2pi*) for the
driving force: -0.2, 1, -1.63, 0.3, 0.3, (sqrt(5) - 1) / 2, 10^5.

Additionally, you can take Fourier transforms to get frequency spectrums, with similar parameter inputs:

<a href="https://imgur.com/okoHHyQ"><img src="https://i.imgur.com/okoHHyQ.jpg" title="source: imgur.com" /></a>

(Inputs used: -0.2, 1, -1.63, 0.3, 0.3, (sqrt(5) - 1) / 2, [1, 0])

The whole codebase has been designed for you to use almost exactly the same inputs with minimal variation for all the functions to
perform analysis. This dosent just make it easier to remeber inputs layouts (Trust me, you memorise entire sets of parameters that are
interesting very quickly) but also makes putting the functions into scripts quite a lot easier.

After you've explored the dynamical system a bit and have a feel for whats going on as you change parameters, start the analysis by
finding the Lyapunov exponents and fractal dimensions of the attractors you create. This will allow you too to find non-chaotic, strange
attractors like I and my supervisor did! The methods are explained in much more detail in both the function descriptions and the final
report in the repo above. (which I strongly recommend you have a quick read through atleast) However, to summarise, you can find the
Lyapunov exponent which is a measure of rate of divergence of the flows which tells you whether the system is chaotic or not, and the
box counting, correlation and Lyapunov dimensions of the attractors which can tell you whether the attractor is geometrically strange or
not. As per the title, your aim is to find sets of parameters which have a negative Lyapunov exponent in the x (position) dimension 
and have the following pattern of fractal dimensions:

<a href="https://imgur.com/dX2Gb2f"><img src="https://i.imgur.com/dX2Gb2f.jpg" title="source: imgur.com" /></a>

However, the Box counting dimension will likely be very close to but greater than 1 and the correlation dimension will be the same
usually. This is explained in the final report, so dont worry too much but as long as your results look something like the above, you've
quite likely found a non-chaotic strange attractor! 

If you'd like to take part in the search for SNA's (strange non-chaotic attractors, not a species of alien sadly) then feel free to use
all my code and trawl through infinite parameter space to your hearts content. If you do come across strong evidence for one, then
please let me know about it! (seriously, email me at hitesh.kumar15@imperial.ac.uk)

## Credits ##
I worked on this project under the supervision of Dr. P. Ramsden, who supported me through many difficult parts of the project and
taught me a lot about dynamical systems, fourier analysis and research in general. Credits for that (and especially for the Box_count
function) go to him.

Additionally, I would like to thank the EPSRC for trusting me and this project and selecting me as one of the few to get their VERY
generous funding. 

## License ##

MIT license, feel free to use and develop for your own investigations!
