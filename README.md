# WaveSim
1D and 2D wave simulation from scratch. For more details, please check out this [YouTube video](https://youtu.be/4IL8n8yYNjw).

![thumbnail](https://raw.githubusercontent.com/tgautam03/WaveSim/master/thumbnail.png)


The repository is setup as a package with source code in the folder named `src`. 

> Before running the demo files, please add this package to your Julia global environment.

- To run the 1D wave equation demo, please open `app_1d.jl` file in Pluto notebook.

- To run the 2D wave equation demo, please open `app_2d.jl` file in Pluto notebook.

### Quick instructions for adding `WaveSim` package to the global environment
> WARNING: The instructions below might not work, so please do check out the official Julia programming language webpage for instructions.

- Open Julia in the terminal and make sure that the working directory is `WaveSim` folder.
- Move to the Pkg manager by pressing `]` key.
- Execute `dev .` command. This should add `WaveSim` to the global environment.
- Go back to the Julia prompt and open Pluto notebooks.


## Topics
- **1D Wave Equation**
    - Finite Difference Method
    - Simulation Accuracy and Stability
        - Numerical Dispersion
        - Fine vs. coarse grid
        - Higher-order Finite Difference approximations
        - Von Neumann Analysis
    - Boundary Conditions
        - Dirichlet boundary condition
        - Neumann boundary conditions
        - Absorbing boundary conditions
- **2D Wave Equation**
    - Finite difference method for 2D problem
    - Heterogeneous medium
    - Simulation Accuracy and Stability
        - Numerical Dispersion
        - Fine vs. coarse grid
        - Higher-order Finite Difference approximations
        - Von Neumann Analysis
    - Boundary Conditions
        - Dirichlet boundary condition
        - Neumann boundary conditions
        - Absorbing boundary conditions

## References
- [FDM 5-point stencil](https://math.stackexchange.com/questions/262701/how-to-obtain-prove-5-stencil-formula-for-2nd-derivative)
- Excellent resources on Von Neumann analysis:
    - [MIT notes](https://math.mit.edu/classes/18.300/Notes/Notes_vNSA.pdf)
    - [Jeff Chasnov's YouTube video](https://www.youtube.com/watch?v=QUiUGNwNNmo)
    - [LMU Seismology](https://www.youtube.com/watch?v=5_VtrWGaEGM)
- [Absorbing boundary condition](https://www.jpier.org/ac_api/download.php?id=0506213)