var documenterSearchIndex = {"docs":
[{"location":"2d_wave/2d_FDM/#Finite-Difference-Method","page":"Finite Difference Method","title":"Finite Difference Method","text":"","category":"section"},{"location":"1d_wave/1d_FDM/#Finite-Difference-Method","page":"Finite Difference Method","title":"Finite Difference Method","text":"","category":"section"},{"location":"#WaveSim.jl","page":"Home","title":"WaveSim.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Earthquakes are one of the most powerful and destructive natural phenomena on Earth. They are caused by the sudden movement or shifting of tectonic plates deep within the Earth's crust, which releases enormous amounts of energy in the form of seismic waves. These waves travel through the ground, causing it to shake and tremble violently.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is my attempt at simulating these waves (accurately) which can be modeled using the Wave Equation","category":"page"},{"location":"","page":"Home","title":"Home","text":"fracpartial^2 p(xyzt)partial t^2 = c^2(xyz) biggfracpartial^2 p(xyzt)partial x^2 + fracpartial^2 p(xyzt)partial y^2 + fracpartial^2 p(xyzt)partial z^2 bigg + s(xyzt)","category":"page"},{"location":"#1D-Wave-Equation","page":"Home","title":"1D Wave Equation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"I will start with the simplest case, i.e. 1D Wave Equation which can be seen as a model for wave on a string (homogenous material)","category":"page"},{"location":"","page":"Home","title":"Home","text":"fracpartial^2 u(xt)partial t^2 = c^2 biggfracpartial^2 u(xt)partial x^2 bigg + s(xt)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The variable u(xt) describes the displacement of string (in the vertical direction) at position x and time t. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Figure 1: Wave on a string. @fig:1)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The constant c defines the speed of wave propagation. In this case, it is dependent on the material property. The higher the value of c, the faster will the wave propagate.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Figure 2: Wave on a string with speed 334 m/s. @fig:2)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Figure 3: Wave on a string with speed 668 m/s. @fig:3)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Using this as an example, I will explain","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finite Difference Method for solving PDEs\nSimulation Accuracy: \nNumerical Dispersion\nHow the number of points per wavelength affect the accuracy?\nHow much can Higher Order Finite Difference Schemes improve the accuracy?\nSimulation Stability:\nWhat exactly is stability (or a stable solution)?\nVon Neumann Analysis\nBoundary Conditions:\nWhy are they so important? (Hint: they ensure a unique solution!)\nUnderstanding different boundary conditions:\nDirichlet boundary condition (Zero/fixed/clamped boundary conditions)\nNeumann boundary conditions (Free boundary conditions)\nAbsorbing boundary conditions","category":"page"},{"location":"#2D-Wave-Equation-(NOT-FINALIZED-YET)","page":"Home","title":"2D Wave Equation (NOT FINALIZED YET)","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"After finishing up the 1D case, I will delve into the 2D simulation, which is modeled using the 2D Wave Equation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"fracpartial^2 p(xyt)partial t^2 = c^2(xy) biggfracpartial^2 p(xyt)partial x^2 + fracpartial^2 p(xyt)partial y^2 bigg + s(xyt)","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is where things get interesting, and I will be discussing","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finite Difference Method to solve 2D Wave Equation in a Homogeneous medium \nSimulation Accuracy:\nNumerical Anisotropy\nFinite Difference Method to solve 2D Wave Equation in a Heterogeneous medium ","category":"page"}]
}
