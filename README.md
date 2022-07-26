# SK-$\tau$-ROCK methods
The purpose of this code is to implement and test the explicit stabilized $\tau$-leap methods, called PSK-$\tau$-ROCK and SK-$\tau$-ROCK, proposed in
- Abdulle, A., Gander, L., & Rosilho de Souza, G. Optimal explicit stabilized postprocessed $\tau$-leap method for the simulation of chemical kinetics. [arXiv: math.NA/2106.09339](https://arxiv.org/abs/2106.09339), 2021.

## Implemented $\tau$-leap methods
In addition to the SK-$\tau$-ROCK schemes, other popular $\tau$-leap methods are implemented for comparison purposes. 
The following methods are available:
- the exact, up to statistical error, Stochastic Simulation Algorithm (`ssa`) [1],
- the standard explicit $\tau$-leap method (`tl`) [2],
- the SK-$\tau$-ROCK method with or without postprocessing (`str`) [3],
- the $\tau$-ROCK method and reversed $\tau$-ROCK methods (`tr` and `rtr`) [4],
- the implicit $\tau$-leap method with or without postprocessing (`itl`) [5],
- the trapezoidal $\tau$-leap method (`ttl`) [6],
- the split-step implicit $\tau$-leap method (`ssitl`) [7].

## Problems already implemented
The following problems are already implemented in the code and can be solved by any of the above methods (see src/ProblemsList.cpp for the details):
 1) Reversible Isomerization, a simple linear scalar model problem [6],
 2) Nonlinear Reversible Reaction [4],
 3) Genetic Positive Feedback Loop [8],
 4) Michaelis-Menten,
 5) Schlogl Reaction [9],
 6) Decaying Dimerizing [2],
 7) E. Coli [10].
 
 ## Installation
 _Prerequisites_: CMake, git and a C++ compiler.  
 
 The code can be compiled with the following commands:
 ```
 chmod u+x configure.sh build.sh clear.sh
 ./configure.sh
 ./build.sh
 ```
 Upon configuration, the script will use git to automatically download the required external libraries: `Eigen` and `GetPot`.
 
 The compiled executable is found in the ``install/`` directory.
 
 ## Running the code
 
To run the code the following options are available:
- `-prob`: the example to solve, the argument n=1,2,...,7 corresponds to the number in the above list.
- `-solver`: the tau-leap method to use, must be one of: `ssa, tl, str, tr, rtr, itl, ttl, ssitl`.
- `-mc`: number of Monte Carlo samples. Choose 1 if you want to simulate only one path and save it. By default is 1.
- `-tau`: the step size. When the solver is `tl`, tau might be reduced during integration due to stability restrictions.
- `-nout`: When mc is 1 it indicates, approximately, the number of outputs.
- `-ofile`: name of output file, by default is _sol_.
- `-pp`: must be 0 or 1 to indicate if postprocessing of str or itl is activated or not. By default is 1.
- `-refsol`: if given, it must correspond to a previously computed solution. The result of the new simulation will be compared the provided one.
- `-Ntol`: tolerance in Newton algorithm of implicit solvers. By default it is 1e-2.
- `-s`: if given, the explicit stabilized methods will not choose the number of stages automatically according to their stability constraint formula but $s$ will be fixed to the given value.
- `-s_add`: add some stages to the minimum required for stability. May be necessary when the stiffness increases within one time step and the algorithm becomes unstable. By default it is 1.
- `-damping`: changes the damping parameter, for `tr` and `rtr` this option is considered only when `-s` is given as well. By default, for `str` damping is 0.05, for `tr` and `rtr` it is chosen according to a formula depending on $s$.

### Examples 
Before running the examples, move into the `install/` folder.

- to run the Genetic Positive Feedback Loop problem, using 1e4 Monte Carlo iterations with the Stochastic Simulation Algorithm and name the output file _ssa_sol_, run:
```
 ./TauLeapMethods -prob 3 -mc 1e4 -solver ssa -ofile ssa_sol
 ```
 Simulation results are found in the `install/GeneticPositiveFeedbackLoop` folder.
- to run the same problem, with tau=0.05, 1e4 Monte Carlo iterations, using SK-$\tau$-ROCK with postprocessing and using 4 additional stages, naming the output file _sk-tau-rock_sol_ and finally compare the result with the previous simulation, run:
```
./TauLeapMethods -prob 3 -tau 0.05 -mc 1e4 -solver str -pp 1 -s_add 4 -ofile sk-tau-rock_sol -refsol ssa_sol
```
- to compute one sample of the Michaelis-Menten problem, using SK-$\tau$-ROCK without postprocessing, with tau=0.001 and 1e3 output points, naming the output file _MM_str_, run:
```
./TauLeapMethods -prob 4 -tau 0.001 -nout 1e3 -mc 1 -solver str -pp 0 -ofile MM_str
```
The solution can be displayed with the `Plot_path.m` script.

## License
See `LICENSE.txt` for licensing information.

## References
- <sub>[1] D. T. Gillespie. Exact stochastic simulation of coupled chemical reactions. Journal of Physical Chemistry, 81(25), 1977.</sub>
- <sub>[2] D. T. Gillespie. Approximate accelerated stochastic simulation of chemically reacting systems. Journal of Chemical Physics, 115(4), 2001.</sub>
- <sub>[3] Abdulle, A., Gander, L., & Rosilho de Souza, G. Optimal explicit stabilized postprocessed tau-leap method for the simulation of chemical kinetics. arXiv: math.NA/2106.09339, 2021.</sub>
- <sub>[4] A. Abdulle, Y. Hu, and T. Li. Chebyshev Methods with Discrete Noise: the τ-ROCK Methods. Journal of Computational Mathematics, 28(2), 2010.</sub>
- <sub>[5] M. Rathinam, L. R. Petzold, Y. Cao, and D. T. Gillespie. Stiffness in stochastic chemically reacting systems: The implicit tau-leaping method. Journal of Chemical Physics, 119(24), 2003.</sub>
- <sub>[6] Y. Cao, L. R. Petzold, M. Rathinam, and D. T. Gillespie. The numerical stability of leaping methods for stochastic simulation of chemically reacting systems. Journal of Chemical Physics, 121(24), 2004.</sub>
- <sub>[7] C. B. Hammouda, A. Moraes, and R. Tempone. Multilevel hybrid split-step implicit tau-leap. Numerical Algorithms, 74(2), 2017.</sub>
- <sub>[8] Y. Yang, M. Rathinam, and J. Shen. Integral tau methods for stiff stochastic chemical systems. Journal of Chemical Physics, 134(4), 2011.</sub>
- <sub>[9] M. Rathinam, L. R. Petzold, Y. Cao, and D. T. Gillespie. Consistency and Stability of Tau-Leaping Schemes. Multiscale Model. Simul., 4(3), 2005.</sub>
- <sub>[10] Hu, Y., Abdulle, A., & Li, T.. Boosted hybrid method for solving chemical reaction systems with multiple scales in time and population size. Communications in Computational Physics, 12(4), 2012.</sub>
