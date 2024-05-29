# Finite Difference Method
Given the Wave Equation, we want to evaluate $u(x,t)$ for $0 \le x \le x_n$ and $0 \le t \le t_n$. For simplicity, let's set $x_n=1$ and $t_n=1$. Right now the domain is continuous and solving the Wave Equation analytically would be quite challenging. This is where we could use Numerical Methods like Finite Difference Method. The very first step is to discretize the domain into finite number of points. 

To keep things simple, I have picked 11 points (each) to represent the $x$ and $t$ dimensions. This transforms the problem where we now have to just evaluate $u$ at these 121 discrete points. Remember that we have to pre-define boundary and initial conditions, which leaves the number of unknown points to 90 (see Figure 1 for details). We will satisfy the Wave Equation at these points and get the value of $u$.

---
```@raw html
<div class="imgcap">
  <img src="https://raw.githubusercontent.com/tgautam03/WaveSim/master/docs/src/img/1d_wave/grid.png" alt="this slowpoke moves"  width="800"/>
  <div class="thecap">Figure 1: Discretized grid with green shaded points showing locations where value of u is known. </div>
</div>
```
---

## Differential to Difference Equation
Partial Differential Equations are represented in the continuous domain, but we just discretized our domain into finite number of points. This is where we need to transform our differential equation into a difference equation (which is the discrete version of a differential equation). 

Taylor Series is a way of approximating a function around some value $x$ using the derivatives of that function 

$$f(x+dx)=f(x)+dx \frac{df}{dx}\big(x\big)+\frac{dx^2}{2!} \frac{d^2f}{dx^2}\big(x\big)+\frac{dx^3}{3!} \frac{d^3f}{dx^3}\big(x\big)+\cdots \tag{1}$$

Looking at the Wave Equation, we can see that it has 2nd order derivatives in time and space.

$$\frac{\partial^2 u(x,t)}{\partial t^2} = c^2 \bigg[\frac{\partial^2 u(x,t)}{\partial x^2} \bigg] + s(x,t) \tag{2}$$

With Taylor Series we can approximate both $\frac{\partial^2 u(x,t)}{\partial t^2}$ and $\frac{\partial^2 u(x,t)}{\partial x^2}$ as follows

$$u(x+dx, t)=u(x, t)+dx \frac{du}{dx}\big(x,t\big)+\frac{dx^2}{2!} \frac{d^2u}{dx^2}\big(x,t\big)+\frac{dx^3}{3!} \frac{d^3u}{dx^3}\big(x,t\big)+\cdots \tag{3a}$$

$$u(x-dx, t)=u(x, t)-dx \frac{du}{dx}\big(x, t\big)+\frac{dx^2}{2!} \frac{d^2u}{dx^2}\big(x, t\big)-\frac{dx^3}{3!} \frac{d^3u}{dx^3}\big(x, t\big)+\cdots \tag{3b}$$

Adding equations $(3a)$ and $(3b)$, we get

$$\frac{u(x+dx, t) - 2u(x, t)+ u(x-dx, t)}{dx^2}= \frac{d^2u}{dx^2}\big(x, t\big)+2\frac{dx^2}{4!} \frac{d^4u}{dx^4}\big(x, t\big)+\cdots \tag{3}$$

Similarly for partial derivative with respect to $t$, we can write

$$\frac{u(x, t+dt) - 2u(x, t)+ u(x, t-dt)}{dt^2}= \frac{d^2u}{dt^2}\big(x,t\big)+2\frac{dt^2}{4!} \frac{d^4u}{dt^4}\big(x, t\big)+\cdots \tag{4}$$

We can ignore the 4th derivative term in equations $(3)$ and $(4)$, as long as we acknowledge that our approximation is $\mathcal{O}(dx^2)$ and $\mathcal{O}(dt^2)$ accurate respectively (commonly known as 2nd order accuracy). Plugging the approximation to the derivatives into the differential equation, we get a difference equation

$$\frac{u(x, t+dt) - 2u(x, t)+ u(x, t-dt)}{dt^2} = c^2 \bigg[\frac{u(x+dx, t) - 2u(x, t)+ u(x-dx, t)}{dx^2} \bigg] + s(x,t) \tag{5a}$$

In Figure 2, we can see how this ties back to our discretized grid. Difference equation $(5a)$ can also be written more succintly as

$$\frac{u_{i}^{j+1} - 2u_{i}^{j}+ u_{i}^{j-1}}{dt^2} = c^2 \bigg[\frac{u_{i+1}^{j} - 2u_i^j+ u_{i-1}^j}{dx^2} \bigg] + s_i^j \tag{5}$$

> Note that I have represented the space dimension in the subscript and the time dimension in the superscript.

---
```@raw html
<div class="imgcap">
  <img src="https://raw.githubusercontent.com/tgautam03/WaveSim/master/docs/src/img/1d_wave/Difference_eq_pts.png" alt="this slowpoke moves"  width="800"/>
  <div class="thecap">Figure 2: Connecting Discretized grid to the difference equations. </div>
</div>
```
---

## Algorithm
If we look carefully, we can solve equation $(5)$ as an extrapolation problem. All we have to do is isolate the future $u$ i.e. $u_i^{j+1}$ on the left hand side and get it's value by evaluating the remaining right hand side expression.

$$u_{i}^{j+1} = c^2 \frac{dt^2}{dx^2} \bigg[u_{i+1}^{j} - 2u_i^j+ u_{i-1}^j \bigg] + 2u_{i}^{j} - u_{i}^{j-1} + s_i^j \tag{5}$$



```julia
# Pressure Field at p(x,t)
p = zeros(nx)
# Pressure Field at p(x,t-dt)
p_prev = zeros(nx)
# Pressure Field at p(x,t+dt)
p_next = zeros(nx)
# Solution at each step
p_sols = zeros(nt,nx)

  # Looping over time
  for it in range(start=2, stop=nt)
      # Looping over Space
      for ix in range(start=2, stop=nx-1)
          # Evaluating 2nd derivative wrt x
          d2p_dx2 = p[ix+1] - 2*p[ix] + p[ix-1]
          # Updating Solution
          if ix == isrc
              p_next[ix] = (c*dt/dx)^2 * d2p_dx2 + 2*p[ix] - p_prev[ix] + dt^2 * src[it]
          else
              p_next[ix] = (c*dt/dx)^2 * d2p_dx2 + 2*p[ix] - p_prev[ix]
          end
      end
  
      # Current Sol becomes Previous Sol
      p_prev[:] = p[:]
      # Next Sol becomes Current Sol
      p[:] = p_next[:]
  
      # Storing solutions at each time step
      p_sols[it,:] = p
  end

  return p_sols
```