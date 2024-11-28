```@raw html
<body>
<div>
  \(
  \definecolor{florange}{rgb}{0.832	,0.367 ,0}
  \definecolor{flblue}{rgb}{0	,0.445	,0.695}
  \definecolor{flgreen}{rgb}{0	,0.617	,0.449}
  \definecolor{flred}{rgb}{1	,0  	,0    }
  \)
</div>
```


# Tests cases

## About testing
cf. [this tutorial](https://erikexplores.substack.com/p/julia-testing-best-pratice)  
cf. [this tutorial](https://erikexplores.substack.com/p/julia-testing-best-pratice)  

## Manufactured solutions


See [`set_poisson`](@ref)

### Poisson equation inside a circle

See 2.4.4 in [`Rodriguez 2024`](https://theses.fr/s384455)

See [`(Rodriguez et al. 2024)`](https://link.springer.com/article/10.1007/s00707-024-04133-4)

``p (x, y) = cos (π^2 x y) sin (π^2 x y)``

### Poisson equation inside a square
```@raw html
<a name="tagPoisson"></a> 
```
The Poisson equation is solved inside a square domain: ``\nabla \cdot \nabla p = f~\text{in}~\Omega``



In this test case, we impose:
* A non-homogeneous Neumann boundary condition on the left wall
* A homogeneous Neumann boundary condition on the top and bottom walls
* A zero Dirichlet boundary condition on the right wall


Recalling [Example away from interfaces](@ref), the divergence in 2D in a cell ``(i,j)`` away from the wall, without an interface, with a regular mesh of constant spacings ``h_x`` and ``h_y`` is:

```math
\begin{aligned}
{\mathrm{div}(q^ω, q^{γ})}_i &= \color{florange}{A^x_{i+\frac{1}{2}} q_{i+\frac{1}{2}}^\omega-A^x_{i-\frac{1}{2}} q_{i-\frac{1}{2}}^\omega} \\
&+\color{flblue}{(B_i^x -A^x_{i+\frac{1}{2}}) q_{i+\frac{1}{2}}^\gamma +(A^x_{i-\frac{1}{2}}-B_i^x ) q_{i-\frac{1}{2}}^\gamma } \\
&+\color{flgreen}{A^y_{j+\frac{1}{2}} q_{j+\frac{1}{2}}^\omega -A^y_{j-\frac{1}{2}} q_{j-\frac{1}{2}}^\omega } \\
&+\color{flblue}{(B_j^y -A^y_{j+\frac{1}{2}}) q_{j+\frac{1}{2}}^\gamma +(A^y_{j-\frac{1}{2}}-B_j^y ) q_{j-\frac{1}{2}}^\gamma } \\
&= \color{flblue}{h_y (\frac{c_{i+1}^\omega-c_i^\omega}{h_x}) - h_y (\frac{c_{i}^\omega-c_{i-1}^\omega}{h_x})} \\
&+ \color{flgreen}{h_x (\frac{c_{j+1}^\omega-c_j^\omega}{h_y}) - h_x (\frac{c_{j}^\omega-c_{j-1}^\omega}{h_y})}
\end{aligned}
```
With ``h_x = h_y``, this simplifies to:

```math
{\mathrm{div}(q^ω, q^{γ})}_{i,j} = -4 c_{i,j} + c_{i-1,j} + c_{i+1,j} + c_{i,j-1} + c_{i,j+1}
```

!!! todo "TODO"
    Dirichlet, Neumann

!!! todo "TODO"
    corners

At the left wall, we have:

```math
\begin{aligned}
{\mathrm{div}(q^ω, q^{γ})}_i &= \color{florange}{A^x_{i+\frac{1}{2}} q_{i+\frac{1}{2}}^\omega-A^x_{i-\frac{1}{2}} q_{i-\frac{1}{2}}^\omega} \\
&+\color{flblue}{(B_i^x -A^x_{i+\frac{1}{2}}) q_{i+\frac{1}{2}}^\gamma +(A^x_{i-\frac{1}{2}}-B_i^x ) q_{i-\frac{1}{2}}^\gamma } \\
&+\color{flgreen}{A^y_{j+\frac{1}{2}} q_{j+\frac{1}{2}}^\omega -A^y_{j-\frac{1}{2}} q_{j-\frac{1}{2}}^\omega } \\
&+\color{flblue}{(B_j^y -A^y_{j+\frac{1}{2}}) q_{j+\frac{1}{2}}^\gamma +(A^y_{j-\frac{1}{2}}-B_j^y ) q_{j-\frac{1}{2}}^\gamma } \\
&= \color{flblue}{h_y (\frac{c_{i+1}^\omega-c_i^\omega}{h_x}) - h_y (\frac{c_{i}^\omega-c_{i-1}^\omega}{\frac{h_x}{2}})} \\
&+ \color{flgreen}{h_x (\frac{c_{j+1}^\omega-c_j^\omega}{h_y}) - h_x (\frac{c_{j}^\omega-c_{j-1}^\omega}{h_y})}
\end{aligned}
```
With ``h_x = h_y``, this simplifies to:

```math
\mathrm{div} = -5 c_{i,j} + 2 c_{i-1,j} + c_{i+1,j} + c_{i,j-1} + c_{i,j+1} 
```

With:

```julia
@time @inbounds @threads for i in 1:A.m
   @inbounds A[i,i] += 1e-10
end
```
the coefficients are no longer exact.

We study the relative errors in the discrete ``l_2`` and ``l_\infty`` norms in 2D:

```math
\begin{aligned}
   l_2&=\sqrt{\frac{\sum{S_i\left( p_i- p_{i}^e\right)}^2}{\sum{S_i\left( p_{i}^e\right)}^2}}
\end{aligned}
```

```math
\begin{aligned}
   l_\infty&=\frac{\max\lvert p_i- p_{i}^e\rvert}{\lvert\max( p_{i}^e)\rvert}
\end{aligned}
```

THe boundary conditions are:
* left: Neumann: +10
* right: p = 0
* top and bottom: homogeneous Neumann


The function is ``p(x,y) = -10(\frac{L}{2}-x)`` with ``L = 2.5`` the domain length.

!!! todo "TODO"
    As you can see in [`set_borders!`](@ref), in the current implementation needs to set the value of the boundary condition times the scalar product ``n \cdot e_x`` for the left and right faces and ``n \cdot e_y`` for the bottom and top faces.


```julia
 @inline function normf_2(field, pos, cap, h)
        AVG = 0.
        RMS = 0.
        VOLUME = 0.
        MAX = 0.
        dv = 0.
        @inbounds for II in pos
            v = abs(field[II])
            dv = cap[II]*h^2
            if dv > 0.
                VOLUME += dv
                AVG += dv*v
                RMS += dv*v^2
                if (v > MAX) MAX = v end
            end
        end
        return AVG/VOLUME, sqrt(RMS/VOLUME), MAX
    end
    
    @inline function normf_2(field, pos, h)
        AVG = 0.
        RMS = 0.
        VOLUME = 0.
        MAX = 0.
        @inbounds for II in pos
            v = abs(field[II])
            VOLUME += h^2
            AVG += (h^2)*v
            RMS += (h^2)*v^2
            if (v > MAX) MAX = v end
        end
        return AVG/VOLUME, sqrt(RMS/VOLUME), MAX
    end
```


### Poisson equation inside a square with circular interface


### Poisson equation inside a square with circular interface at wall

### Gradient



## Poiseuille (TODO)

### Epsilon to prevent NaN




with epsilon in capacities, ... cf `num.epsilon_mode`: coefficients not exact 

after removing eps : Laplacian 1 1 -4 1 1
1 1 -5 1 1 2: exact coefficients (like with finite differences) at machine error


### Taylor expansion
The mesh spacings are h/2 then h at left border. 

With ``v_{i+1}`` ``v_i`` and ``v_{i-1}`` with h and h/2 spacings (we can assimilate the BC to ``v_{i-1}``)


```math
\begin{aligned}
&v_{i+1}=v(x)+\frac{h_1}{1!}v'(x)+\frac{{h_1}^2}{2!}v''(x)+\frac{{h_1}^3}{3!}v^{(3)}(x)+\frac{{h_1}^4}{4!}v^{(4)}(x)\\
&v_i\\
&v_{i-1}=v(x)-\frac{h_2}{1!}v'(x)+\frac{{h_2}^2}{2!}v''(x)-\frac{{h_2}^3}{3!}v^{(3)}(x)+\frac{{h_2}^4}{4!}v^{(4)}(x)\\
\end{aligned}
```

 then I get that the coefficients for the second order derivative are ``4/3(1*v_{i+1} - 3v_{i} +2v_{i-1})``. 
 
 In Flower the coefficients at the BC, away from the wall are 1 -5 1 1 and 2 for the BC (like 1 -2 1 for y and 1 -3 2 for x). There is a factor 4/3 missing and in fact when I multiply the result given by Flower by 4/3, I have the exact value, like I have away from the wall: rho/Re*Laplacian*v = imposed pressure gradient
```-3.3882972000e+02 * 4/3 = -4.5177296000e+02```.

The scheme is not second order at the wall because of the irregular spacing. 

Here it is exact . 
There is a factor ``\frac{4}{3}`` missing at the wall (bulk and border) without interface to get the exact Laplacian (because the Poiseuille profile is a second-order polynomial).

For the interface, there would be a factor ``\frac{2h}{h+dist}``.




when I do a Taylor expansion with ``v_{i+1}`` ``v_i`` and ``v_{i-1}`` with h and h/2 spacings (if I understood correctly, we can assimilate the BC to ``v_{i-1}``, then I get that the coefficients for the second order derivative are ``4/3(1*v_{i+1} - 3v_{i} +2v_{i-1})``. In Flower the coefficients at the BC, away from the wall are 1 -5 1 1 and 2 for the BC (like 1 -2 1 for y and 1 -3 2 for x). I think there is a factor 4/3 missing and in fact when I multiply the result given by Flower by 4/3, I have the exact value, like I have away from the wall: rho/Re*Laplacian*v=-3.3882972000e+02 * 4/3 = -4.5177296000e+02.



Normally, we have a b c coefficients such that ``a*h1 -c*h2 = 0`` in ``a*v_{i+1} + b*v_i +c*v_{i-1}``.

At the moment, we get ``3/4*L +\mathcal{O}(h)``  at the corners and ``L + \mathcal{O}(h^2)`` in the bulk when there is no interface


