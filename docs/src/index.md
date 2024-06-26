# WaveSim.jl
Earthquakes are one of the most powerful and destructive natural phenomena on Earth. They are caused by the sudden movement or shifting of tectonic plates deep within the Earth's crust, which releases enormous amounts of energy in the form of **seismic waves**. These waves travel through the ground, causing it to shake and tremble violently.

This is my attempt at simulating these waves (accurately) which can be modeled using the Wave Equation

$$\frac{\partial^2 p(x,y,z,t)}{\partial t^2} = c^2(x,y,z) \bigg[\frac{\partial^2 p(x,y,z,t)}{\partial x^2} + \frac{\partial^2 p(x,y,z,t)}{\partial y^2} + \frac{\partial^2 p(x,y,z,t)}{\partial z^2} \bigg] + s(x,y,z,t)$$

## 1D Wave Equation
### Introduction
I will start with the simplest case, i.e. 1D Wave Equation which can be seen as a model for wave on a string (homogenous material)

$$\frac{\partial^2 u(x,t)}{\partial t^2} = c^2 \bigg[\frac{\partial^2 u(x,t)}{\partial x^2} \bigg] + s(x,t)$$

The variable $u(x,t)$ describes the displacement of string (in the vertical direction) at position $x$ and time $t$. 

---
```@raw html
<div class="imgcap">
  <img src="https://raw.githubusercontent.com/tgautam03/WaveSim/master/docs/src/img/index/1d_wave.png" alt="this slowpoke moves"  width="800"/>
  <div class="thecap">Figure 1: Wave on a string. </div>
</div>
```
---
\
The constant $c$ defines the speed of wave propagation. In this case, it is dependent on the material property. The higher the value of $c$, the faster will the wave propagate.

---
```@raw html
<div class="imgcap">
  <img src="https://raw.githubusercontent.com/tgautam03/WaveSim/master/docs/src/img/index/1D_c1.gif" alt="this slowpoke moves"  width="800"/>
  <div class="thecap">Figure 2: Wave on a string with speed 334 m/s. </div>
</div>
```
---
---

```@raw html
<div class="imgcap">
  <img src="https://raw.githubusercontent.com/tgautam03/WaveSim/master/docs/src/img/index/1D_c2.gif" alt="this slowpoke moves"  width="800"/>
  <div class="thecap">Figure 3: Wave on a string with speed 668 m/s. </div>
</div>
```
---
\
Finally, $s(x,t)$ is the source function and it is responsible for initiating the wave at location $x$ and time $t$. The definition of source function also dictates the shape of the resulting wave (see Figures 4 and 5). 

---
```@raw html
<div class="imgcap">
  <img src="https://raw.githubusercontent.com/tgautam03/WaveSim/master/docs/src/img/index/1D_deriv_gauss_src.gif" alt="this slowpoke moves"  width="800"/>
  <div class="thecap">Figure 4: Wave on a string using the derivative of Gaussian as a source function. </div>
</div>
```
---
---
```@raw html
<div class="imgcap">
  <img src="https://raw.githubusercontent.com/tgautam03/WaveSim/master/docs/src/img/index/1D_gauss_src.gif" alt="this slowpoke moves"  width="800"/>
  <div class="thecap">Figure 5: Wave on a string using the Gaussian as a source function. </div>
</div>
```
---
\
One thing to remember is that partial differential equations (PDEs) are typically defined for the interior of a domain, not the boundary. This holds true for both spatial and temporal domains where we have to explicitly define 
- Boundary Conditions: $u(x,t)$ at $x=0$ and $x=x_{max}$ for **all $t$ values!**
- Initial Condition: $u(x,t)$ at $t=0$ for **all $x$ values!**

For example, if our spatial and temporal domains range from 0 to 1, and we can say that 
- Boundary Conditions: $u(x=0,t)=u(x=1,t)=0$, 
- Initial Condition: $u(x,t=0)=0$ 

then the problem formulation looks something like
- For $t=0$
  - The solution is $u(x,0) = 0$
- For $t>0$
  - if $x=0$ or $x=1$
    - The solution is $u(0,t)=u(1,t)=0$
  - if $0<x<1$
    - We get the solution by solving $\frac{\partial^2 u(x,t)}{\partial t^2} = c^2 \bigg[\frac{\partial^2 u(x,t)}{\partial x^2} \bigg] + s(x,t)$

### Table of Contents
Using this as an example, I have explained in detail
- Finite Difference Method for solving PDEs
- Simulation Accuracy: 
    - Numerical Dispersion
    - How the number of points per wavelength affect the accuracy?
    - How much can Higher Order Finite Difference Schemes improve the accuracy?
- Simulation Stability:
    - What exactly is stability (or a stable solution)?
    - Von Neumann Analysis
- Boundary Conditions:
    - Why are they so important? (Hint: they ensure a unique solution!)
    - Understanding different boundary conditions:
        - Dirichlet boundary condition (Zero/fixed/clamped boundary conditions)
        - Neumann boundary conditions (Free boundary conditions)
        - Absorbing boundary conditions

## 2D Wave Equation (NOT FINALIZED YET)
After finishing up the 1D case, I will delve into the 2D simulation, which is modeled using the 2D Wave Equation.

$$\frac{\partial^2 p(x,y,t)}{\partial t^2} = c^2(x,y) \bigg[\frac{\partial^2 p(x,y,t)}{\partial x^2} + \frac{\partial^2 p(x,y,t)}{\partial y^2} \bigg] + s(x,y,t)$$

This is where things get interesting, and I will be discussing
- Finite Difference Method to solve 2D Wave Equation in a Homogeneous medium 
- Simulation Accuracy:
    - Numerical Anisotropy
- Finite Difference Method to solve 2D Wave Equation in a Heterogeneous medium 