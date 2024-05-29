# Simulation Accuracy
We some assumptions (mainly discretization) while implementing the Finite Difference Method, so how do we know the accuracy of the numerical solution? 

The Wave Equation with no source term (and some non-zero initial condition) has an analytical solution such that the wave at $t=0$ (i.e. initial condition) propagates as it is along the length with time.

$$\frac{\partial^2 u(x,t)}{\partial t^2} = c^2 \bigg[\frac{\partial^2 u(x,t)}{\partial x^2} \bigg] \tag{1a}$$

Initial Conditions:

$$u(x,t=0)=u_0 \\ \frac{\partial u}{\partial t}\big(x,t=0\big)=0 \tag{1b}$$

Then the final solution is given as

$$u(x,t)=\frac{1}{2}\bigg[u_0(ct-x) + u_0(ct+x)\bigg]$$

However, if there is a source term, we need to use the concept of Green's Function and solve the equation with a delta function (in space time) as the source. I won't go in any more details but just tell you that after that analysis, we observe that the solution is just the integral of the source function propagating with time.

