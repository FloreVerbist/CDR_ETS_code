# CDR-ETS integration code


This code belongs to the paper:
F. Verbist, J. Meus, J. A. Moncada, P. Valkering, and E. Delarue, “Carbon removals meet Emission Trading System design: A precautionary path towards integration,” Energy Economics, vol. 145, p. 108389, May 2025, doi: 10.1016/j.eneco.2025.108389.

The code contains the 8 different CDR-ETS integration cases. They are modelled as Mixed Complementarity Problems (MCP) and solved using an ADMM-based optimisation algorithm in Julia, JuMP. 


---

## Repository Structure

```The code needs to be accessed and ran via the 'main.jl' file. This will load all parameters first in the parameters.jl file. A user can change the parameters to customise the code.
The case_processing.jl file will initialise the optimisation models. The Cnr-case_name.jl allow case by case to solve the CDR-ETS optimisation framework.
Results and plots are generated using the results.jl script file and corresponding functions. 

.
├── main.jl
├── parameters.jl
├── case_processing.jl
├── Cnr-case_name.jl
└── Results.jl  

---

## Requirements

```This code uses the Gurobi solver to solve each MCP. Using: 

using Pkg 
Pkg.activate("./Environment_CDR")

allows to download the required julia packages. 

