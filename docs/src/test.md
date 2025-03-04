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


## Summary

The following test cases are planned:
```@raw html
<div w3-include-html="./assets/test_table.html"></div> 
```


```@raw html
<table class="styled-table">

    <thead>
        <tr>
            <th> Name </th>
            <!-- <td>Poisson_square</td> --> 
            <!-- or use td if you do not want bold -->
            <th>Poisson_square (steady) </th>
            <th>Poisson_square_circle (steady) </th>
            <th>Poisson_square_circle_arc (steady)</th>
            <th>Poisson_square_circle_arrow (unsteady)</th>
            <th>Poisson_square_circle_arc_arrow (unsteady)</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <th> </th>
            <td> <a href="test.html#Poisson-equation-inside-a-square"><img class="mid" src="./assets/Poisson_square.svg" width="75vh"> </a> </td>
            <td><a href="test.html#Poisson-equation-inside-a-square-with-circular-interface"><img src="./assets/Poisson_square_circle.svg" width="75vh"></a></td>
            <td><a href="test.html#Poisson-equation-inside-a-square-with-circular-interface-at-wall"><img src="./assets/Poisson_square_circle_arc.svg" width="75vh"></a></td>
            <td><img src="./assets/Poisson_square_circle_arrow.svg" width="75vh" ></td>
            <td><img src="./assets/Poisson_square_circle_arc_arrow.svg" width="75vh" ></td>
        </tr>
 
        <tr>

        <th>Poisson: $ \nabla \cdot (\kappa \nabla \phi)=0$ </th>
        <td> <a href="https://github.com/flnt/Flower.jl/blob/electrolysis/test/poisson_no_interface.jl"> ✔️ </a> </td>
        <td>  </td>
        <td> </td>
        <td> </td>
        <td> </td>
        </tr>

        <tr>
        <th>Diffusion</th>
        <td> </td>
        <td> </td>
        <td> </td>
        <td> </td>
        <td> </td>
        </tr>

        <tr>
        <th>Convection</th>
        <td> </td>
        <td> </td>
        <td> </td>
        <td> </td>
        <td> </td>
        </tr>

        <!-- and so on... ✔️ ❌ • • •-->
    </tbody>
</table>
```


!!! info "About testing in Julia"
    cf. [this tutorial](https://erikexplores.substack.com/p/julia-testing-best-pratice)  

## Unit tests Orientation, operators

* normal defined by ``\alpha``: pointing towards domain 

* [Orientation](https://github.com/flnt/Flower.jl/blob/electrolysis/test/orientation.jl): tests gradient, interpolation, divergence (work in progress)

* [Poisson](https://github.com/flnt/Flower.jl/blob/electrolysis/test/poisson_no_interface.jl): tests matrix coefficients at bulk, left wall away from corners, bottom left corner, computes convergence order (TODO test value of slope, level of error)

* [Poisson](https://github.com/flnt/Flower.jl/blob/electrolysis/test/poisson_no_interface_right_Neumann.jl): tests matrix coefficients at bulk, left wall away from corners, bottom left corner, computes convergence order (TODO test value of slope, level of error), Neumann at right

```@raw html
See this equation for the signs given to the g function at the
<a href="documentation.html#tagdivleftwall">left wall</a> and the <a href="documentation.html#tagdivbubbleinterface"> bubble interface</a>
```
For a regular mesh with constant mesh spacing, with ``h_x=h_y``
* Test that the matrix coefficients away from the interface for a Laplacian are -4, 1, 1, 1, 1
* Test that the matrix coefficients at the left wall, away from corners for a Laplacian are -5, ...
* Test that the matrix coefficients at the bottom left corner for a Laplacian are -6, ...
* Test that the matrix coefficients at the bottom left corner for a Laplacian in the BC system with Neumann at left wall are 2 for bulk, -2 for interfacial value
* @testset "interpolate_scalar_to_staggered_u_v_grids_at_border!": test the interpolation of electrical conductivity at wall

### Interpolation of conductivity
We have at the wall an average conductivity, see test "interpolate_scalar_to_staggered_u_v_grids_at_border!" in [Orientation](https://github.com/flnt/Flower.jl/blob/electrolysis/test/orientation.jl)


### Interpolation
To check [`interpolate_scalar!`](@ref)
```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/test/runtests_operators.jl
```


Example

```@raw html
As seen in <a href="documentation.html#Coefficients-in-a-simple-configuration"> </a>:
```

```@raw html
<table class="styled-table">

    <thead>
        <tr>
        <th> BC </th>
        <th>Border </th>
        <th>Bulk </th>
        </tr>
    </thead>

    <tbody>

    <tr>
        <th>Left</th>
        <td> +2</td>
        <td> -2</td>
    </tr>

    <tr>
        <th>Right</th>
        <td> +2</td>
        <td> -2</td>
    </tr>

    <tr>
        <th>Bottom</th>
        <td> +2</td>
        <td> -2</td>
    </tr>

    <tr>
        <th>Top</th>
        <td> +2</td>
        <td> -2</td>
    </tr>

    </tbody>
</table>
```
```@raw html
The coefficients of the matrix for the discretization of a Neumann boundary condition are tested in<a href="https://github.com/flnt/Flower.jl/blob/electrolysis/test/poisson_no_interface.jl"> poisson_no_interface.jl</a> and <a href="https://github.com/flnt/Flower.jl/blob/electrolysis/test/poisson_no_interface_right_Neumann.jl"> poisson_no_interface_right_Neumann.jl </a> with tmp the index of a bulk value and i_corner_vecb_bottom the index of the boundary value at the bottom:
```


```julia
@testset "A[i_corner_vecb_bottom,tmp]" begin
    @test A[i_corner_vecb_bottom,tmp] ≈ -2.0 atol = test_tolerance
end

@testset "A[i_corner_vecb_bottom,i_corner_vecb_bottom]" begin
    @test A[i_corner_vecb_bottom,i_corner_vecb_bottom] ≈ 2.0 atol = test_tolerance
end
```



!!! todo "Orientation"

    vector given by ``(\cos(\alpha),\sin(\alpha))`` is the interior normal so ``\alpha+\pi`` should be used to have the exterior normal
    ```julia
    function neumann_bcs!(gp, N)
    @unpack x, y, dx, dy, LS, ind = gp

    @inbounds @threads for II in ind.inside
        x_bc = LS[1].mid_point[II].x * dx[II] + x[II]
        y_bc = LS[1].mid_point[II].y * dy[II] + y[II]
        Nx = ∇fx(x_bc, y_bc)
        Ny = ∇fy(x_bc, y_bc)

        N[II] = Nx * cos(LS[1].α[II]+π) + Ny * sin(LS[1].α[II]+π)
    end


    return nothing
    end
    ```

## Definition of errors

We study the relative errors in the discrete $l_1$, $l_2$ and $l_\infty$ norms in 2D, with $S_i$ the surface of the cell and $p_i^e$ the analytical solution:

```math
\begin{equation}
\begin{aligned}
l_1&=\frac{\sum{S_i \lvert p_i- p_{i}^e \rvert } }{\sum{S_i\lvert p_{i}^e \rvert}}
\end{aligned}
\end{equation}
```

```math
\begin{equation}
\begin{aligned}
l_2&=\sqrt{\frac{\sum{S_i\left( p_i- p_{i}^e\right)}^2}{\sum{S_i\left( p_{i}^e\right)}^2}}
\end{aligned}
\end{equation}
```

```math
\begin{equation}
\begin{aligned}
l_\infty&=\frac{\max\lvert p_i- p_{i}^e\rvert}{\max( \lvert p_{i}^e \rvert )}
\end{aligned}
\end{equation}
```



## Manufactured solutions


See [`set_poisson`](@ref)

## Poisson equation inside a circle

See 2.4.4 in [`Rodriguez 2024`](https://theses.fr/s384455)

See [`(Rodriguez et al. 2024)`](https://link.springer.com/article/10.1007/s00707-024-04133-4)

``p (x, y) = cos (π^2 x y) sin (π^2 x y)``

### Axisymmetric case

``p(x, y) = x^2 + y^2``

## Poisson equation inside a square
```@raw html
<a name="tagPoisson"></a> 
```
The Poisson equation is solved inside a square domain: ``-\nabla \cdot \nabla p = f~\text{in}~\Omega``.

With:
```
function f(x, y)
    return -10 * (1.25-x) #L=2.5 from -L/2 to L/2
end
```

The analytical function is: ``g = -f``.



See [`solve_poisson_variable_coeff!`](@ref) in [Electrical potential](@ref) and 
```@raw html
 <a href="documentation.html#Poisson-equation"> Poisson equation section </a>.
```

The boundary conditions are:
* left: Neumann: +10
* right: p = 0
* top and bottom: homogeneous Neumann

```@raw html
<figure>
    <a name="Poisson_square"></a> 
    <img src="./assets/test_case.svg" alt="Poisson" title="Poisson">
    <figcaption>"Poisson in a square domain"</figcaption>
</figure>
```

The function is ``p(x,y) = -10(\frac{L}{2}-x)`` with ``L = 2.5`` the domain length.

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




```@raw html
Recalling <a href="documentation.html#Poisson-equation"> Poisson equation section </a>
```


```math
\begin{equation}
    \left [ \begin{array}{>{\centering\arraybackslash$} p{2.0cm} <{$} >{\centering\arraybackslash$} p{3.2cm} <{$}}
    G ^ \top W ^ \dagger G & G ^ \top W ^ \dagger H \\
    I _ b H ^ \top W ^ \dagger G & I _ b H ^ \top W ^ \dagger H + I _ a I _ \Gamma
    \end{array} \right ] \left [ \begin{array}{c}
    p ^ \omega \\
    p ^ \gamma
    \end{array} \right ] \simeq \left [ \begin{array}{c}
    V f ^ \omega \\
     I _ \Gamma g ^ \gamma
    \end{array} \right ],
\label{eq:roblapmat}
\end{equation}
```


"
Appendix One-dimensional Poisson's equation from [`(Rodriguez et al. 2024)`](https://link.springer.com/article/10.1007/s00707-024-04133-4)

 

```@raw html
The discretized  <a href="documentation.html#tagPoisson"> Poisson's equation </a> in terms of the cut-cell operators reads:
```

 
```math
\begin{equation}
    G ^ \top W ^ \dagger G p ^ \omega + G ^ \top W ^ \dagger H p ^ \gamma = V f ^ \omega,
\end{equation}
```
with the one-dimensional version expanding to 
```math
\begin{equation}
    - \left [ B _ x D _ x ^ + W _ x ^ \dagger D _ x ^ - B _ x p ^ \omega + B _ x D _ x ^ + W _ x ^ \dagger \left ( A _ x D _ x ^ - - D _ x ^ - B _ x \right ) p ^ \gamma \right ] = V f ^ \omega
\end{equation}
```
in terms of the elementary discrete operators, if the x-direction is considered. We can now further expand this expression for a given cell i, which yields the following expression
```math
\begin{multline}
    - \mathcal{B} _ {x, i} \left [ \mathcal{W} _ {x, i + 1} \left ( \mathcal{B} _ {x, i + 1} p _ {i + 1} ^ \omega - \mathcal{B} _ {x, i} p _ {i} ^ \omega \right ) - \mathcal{W} _ {x, i}  \left ( \mathcal{B} _ {x, i} p _ {i} ^ \omega - \mathcal{B} _ {x, i - 1} p _ {i - 1} ^ \omega \right ) \right ] - \\
    \mathcal{B} _ {x, i} \left \{ \mathcal{W} _ {x, i + 1} \left [ \left ( \mathcal{B} _ {x, i + 1} - \mathcal{A} _ {x, i + 1} \right )  p _ {i + 1} ^ \gamma + \left ( \mathcal{A} _ {x, i + 1} - \mathcal{B} _ {x, i} \right ) p _ {i} ^ \gamma \right ] - \right . \\
    \left . \mathcal{W} _ {x, i} \left [ \left ( \mathcal{B} _ {x, i} - \mathcal{A} _ {x, i} \right ) p _ {i} ^ \gamma + \left ( \mathcal{A} _ {x, i} - \mathcal{B} _ {x, i - 1} \right ) p _ {i - 1} ^ \gamma \right ] \right \} = \mathcal{V} _ i f _ i ^ \omega.
\end{multline}
```
As seen in the latter equation, the discrete Laplacian yields a 3-point stencil in both ``p ^ \omega`` and ``p ^ \gamma``.
"

#### Equations


By printing the coefficients we see we have -4 when there is no interface. What is discretized is ``\nabla \cdot \nabla p = f``.

!!! todo "TODO"
    ```julia
    A[1:ni,end-nb+1:end] = bc_L_b
    ```
!!! todo "TODO"
    Hx Hb...



### Implementation

!!! info "Changes"
    * Replaced a1=-1 (i.e. a in Robin BC) by a1=1
    * Replaced -chi * a by +chi * a
    
```@docs
set_borders_poisson!
solve_poisson
```


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

### Convergence study
L = 2.5
```math
p_e(x,y) = cos((1.25+x)*(π /5.0))
```


#### n = 16
number_small_cells_for_error 000 
0.0% of mixed cells
ALL: (0.0008035776793641828, 6.457370867720527e-7, 0.0008035776793634693)
MIXED: (NaN, NaN, NaN)
FULL: (0.0008035776793641828, 6.457370867720527e-7, 0.0008035776793634693)

#### n =32
number_small_cells_for_error 000 
0.0% of mixed cells
ALL: (0.00020082180972763426, 4.0329399262224025e-8, 0.0002008218097276315)
MIXED: (NaN, NaN, NaN)
FULL: (0.00020082180972763426, 4.0329399262224025e-8, 0.0002008218097276315)

#### n = 64

number_small_cells_for_error 000 
0.0% of mixed cells
ALL: (5.0200915977696063e-5, 2.520131964877533e-9, 5.0200915972993364e-5)
MIXED: (NaN, NaN, NaN)
FULL: (5.0200915977696063e-5, 2.520131964877533e-9, 5.0200915972993364e-5)

Order 2 for a cos, but careful with the function, Taylor development : ``1 - x^2 +\mathcal{O}(h^4)`` near zero   




```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/test/runtests_manufactured_solutions.jl

python3 -c "import convergence_study; convergence_study.plot_errors_from_h5()" ../Flower.jl/test/poisson_square_solve_poisson.yml convergence_study.h5
```

!!! todo "RECHECK values"


```@raw html
<!-- <table border="1" class="dataframe"> -->
<table class="styled-table">
  <thead>
    <tr style="text-align: right;">
      <th>nx_list</th>
      <th>l1_rel_error</th>
      <th>l2_rel_error</th>
      <th>linfty_rel_error</th>
      <th>1/n</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>16</td>
      <td>1.44e-02</td>
      <td>1.42e-02</td>
      <td>1.47e-02</td>
      <td>6.25e-02</td>
    </tr>
    <tr>
      <td>32</td>
      <td>3.59e-03</td>
      <td>3.53e-03</td>
      <td>3.59e-03</td>
      <td>3.12e-02</td>
    </tr>
    <tr>
      <td>64</td>
      <td>8.96e-04</td>
      <td>8.80e-04</td>
      <td>8.93e-04</td>
      <td>1.56e-02</td>
    </tr>
    <tr>
      <td>128</td>
      <td>2.24e-04</td>
      <td>2.20e-04</td>
      <td>2.23e-04</td>
      <td>7.81e-03</td>
    </tr>
    <tr>
      <td>256</td>
      <td>5.60e-05</td>
      <td>5.50e-05</td>
      <td>5.57e-05</td>
      <td>3.91e-03</td>
    </tr>
  </tbody>
</table>
```

```@raw html
<figure>
    <a name="Poisson_square_results"></a> 
    <img src="./assets/Poisson_square_errors.svg" alt="Poisson" title="Poisson">
    <figcaption>"Poisson in a square domain"</figcaption>
</figure>
```


## Poisson equation inside a square with circular interface

In the current implementation of the test, we did not add the tests for matrix coefficients, they are tested in the case without interface. With an interface, depending on the posittion of the interface and the resolution, the coefficients may no longer be the same as in classical finite differences (-4,1,1,1,1), ...

Like in [`(Rodriguez et al. 2024)`](https://link.springer.com/article/10.1007/s00707-024-04133-4), the analytical solution is imposed at the centroids of the interface.

```@raw html
<figure>
    <a name="Poisson_square"></a> 
    <img src="./assets/test_case_circle.svg" alt="Poisson" title="Poisson">
    <figcaption>"Poisson in a square domain"</figcaption>
</figure>
```


### Dirichlet boundary condition on interface



### Neumann BC


!!! todo "Not symmetrical"
    check why not symmetrical
    check sign BC interface

Relative error in ``l_\infty`` norm on interface: ``\approx \mathcal{O}(10^{?})``


A unit test should be added to check the signs: 

[Robin boundary conditions](https://en.wikipedia.org/wiki/Robin_boundary_condition)


## Poisson equation inside a square with circular interface at wall

```bash
python3 -c "import convergence_study; convergence_study.plot_errors_from_h5()" ../Flower.jl/test/poisson_circular_interface_wall_Dirichlet.yml poisson_circular_interface_wall_Dirichlet.h5
```

```bash
python3 -c "import plot_flower; plot_flower.plot_all_fig_func()" ../Flower.jl/test/poisson_circular_interface_wall_Neumann_256.yml poisson_circular_interface_wall_Neumann_00000256.h5
```

## Gradient

!!! todo "TODO"
    Test gradient, need to recheck localization/ error levels 




```@raw html
<figure>
    <a name="Poisson_square_wall"></a> 
    <img src="./assets/Poisson_square_circle_arc.svg" alt="Poisson" title="Poisson">
    <figcaption>"Poisson in a square domain"</figcaption>
</figure>
```

### Dirichlet on interface and wall

### Neumann on interface and wall

### Dirichlet on interface and  Neumann at wall

### Neumann on interface and Dirichlet at wall



## Butler-Volmer
We use a similar approach as in [`example for sinh`](https://math.stackexchange.com/questions/3472880/solving-sinh-x-kx). We assume a constant conductivity ``\sigma = 1`` and a simplified Butler-Volmer equation.
The boundary conditions are:
* left: simplified Butler-Volmer 

```math
2i_0 \mathrm{sinh}\left( \right)
```

```math
i=i_0\left[\exp{\left(\frac{\alpha_aF\eta}{R_\mu T}\right)}
-\exp{\left(\frac{\alpha_c F\eta}{R_\mu T}\right)}\right]
```

* right: ``\Phi = 0``
* bottom and top ``\frac{\partial \Phi}{\partial n} = 0``

[`example for sinh`](https://math.stackexchange.com/questions/3472880/solving-sinh-x-kx)

!!! todo "Needs to be checked" 
    Poisson with constant conductivity, electroneutrality: ``\Delta \Phi = 0`` so ``\Phi(x) = ax+b``.
    At ``x = L``, ``\Phi = 0`` so we have ``\Phi(x) = a(x-L)`` 
    
!!! todo ""
    with ``a>0`` because 

!!! todo
    At ``0^+``: ``\Phi(0^+) = -aL`` and ``\frac{\partial\Phi}{\partial x }(0^+) = a`` so
    ``a = 2 i_0 \mathrm{sinh}(\frac{\alpha F}{RT}(-0.6-(-aL)) )``

```@raw html
The results of <a href="https://github.com/flnt/Flower.jl/blob/electrolysis/test/butler.jl"> butler.jl</a> run by <a href="https://github.com/flnt/Flower.jl/blob/electrolysis/test/runtest_butler.jl"> runtest_butler.jl </a> can be compared with <a href="https://github.com/flnt/Flower.jl/blob/electrolysis/example/solve_potential_1D.py"> solve_potential_1D.py</a>.
```

In Flower, the following potential is obtained at the electrode: 

```@raw html
<table class="styled-table">
<!-- <table border="1" class="dataframe"> -->
    <table>
        <thead>
            <tr>
                <th>Index</th>
                <th>Value</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td>0</td>
                <td>0</td>
            </tr>
            <tr>
                <td>1</td>
                <td>-0.01411703</td>
            </tr>
            <tr>
                <td>2</td>
                <td>-0.01119343</td>
            </tr>
            <tr>
                <td>3</td>
                <td>-0.0117445</td>
            </tr>
            <tr>
                <td>4</td>
                <td>-0.01163859</td>
            </tr>
            <tr>
                <td>5</td>
                <td>-0.01165887</td>
            </tr>
            <tr>
                <td>6</td>
                <td>-0.01165499</td>
            </tr>
            <tr>
                <td>7</td>
                <td>-0.01165573</td>
            </tr>
        </tbody>
    </table>
```

!!! todo "Proof for only one fixed point in 
\citet{olmoThermodynamicDerivationNernst2019}"


Newton's method with python 1D

```@raw html
<table class="styled-table">
<!-- <table border="1" class="dataframe"> -->
  <thead>
    <tr style="text-align: right;">
      <th>k</th>
      <th>phi_wall</th>
      <th>Relative error on slope</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>0</td>
      <td>-1.146e-02</td>
      <td>1.694e-02</td>
    </tr>
    <tr>
      <td>1</td>
      <td>-1.166e-02</td>
      <td>4.429e-06</td>
    </tr>
    <tr>
      <td>2</td>
      <td>-1.166e-02</td>
      <td>3.003e-13</td>
    </tr>
    <tr>
      <td>3</td>
      <td>-1.166e-02</td>
      <td>7.315e-16</td>
    </tr>
  </tbody>
</table>
```

Successive substitutions: with Flower or python in 1D:
```@raw html
<!-- <table border="1" class="dataframe"> -->
<table class="styled-table">
  <thead>
    <tr style="text-align: right;">
      <th>k</th>
      <th>phi_wall</th>
      <th>Relative error on slope</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>0</td>
      <td>-1.412e-02</td>
      <td>2.112e-01</td>
    </tr>
    <tr>
      <td>1</td>
      <td>-1.119e-02</td>
      <td>3.965e-02</td>
    </tr>
    <tr>
      <td>2</td>
      <td>-1.174e-02</td>
      <td>7.626e-03</td>
    </tr>
    <tr>
      <td>3</td>
      <td>-1.164e-02</td>
      <td>1.460e-03</td>
    </tr>
    <tr>
      <td>4</td>
      <td>-1.166e-02</td>
      <td>2.798e-04</td>
    </tr>
    <tr>
      <td>5</td>
      <td>-1.165e-02</td>
      <td>5.360e-05</td>
    </tr>
    <tr>
      <td>6</td>
      <td>-1.166e-02</td>
      <td>1.027e-05</td>
    </tr>
    <tr>
      <td>7</td>
      <td>-1.166e-02</td>
      <td>1.968e-06</td>
    </tr>
    <tr>
      <td>8</td>
      <td>-1.166e-02</td>
      <td>3.770e-07</td>
    </tr>
    <tr>
      <td>9</td>
      <td>-1.166e-02</td>
      <td>7.223e-08</td>
    </tr>
  </tbody>
</table>
```

Successive substitutions

```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/examples/convergence_Butler.jl ../Flower.jl/examples/convergence_Butler.yml
```

Newton-Raphson

```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/examples/convergence_Butler.jl ../Flower.jl/examples/convergence_Butler_Newton.yml
```

```bash
python3 -c "import convergence_study; convergence_study.plot_convergence_study_func()" ../Flower.jl/examples/convergence_Butler_Newton.yml poisson_00000032_0000000*
```

```@raw html
<table class="styled-table">
  <thead>
    <tr style="text-align: right;">
      <th>Iteration</th>
      <th>Residual</th>
      <th>Variation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>1</td>
      <td>2.013e-02</td>
      <td>1.146e-02</td>
    </tr>
    <tr>
      <td>2</td>
      <td>5.278e-06</td>
      <td>1.974e-04</td>
    </tr>
    <tr>
      <td>3</td>
      <td>3.701e-13</td>
      <td>5.163e-08</td>
    </tr>
    <tr>
      <td>4</td>
      <td>7.739e-15</td>
      <td>3.508e-15</td>
    </tr>
  </tbody>
</table>
```


Succesive substitutions

```@raw html
<table class="styled-table">
   <thead>
    <tr style="text-align: right;">
      <th>Iteration</th>
      <th>Residual</th>
      <th>Variation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>1</td>
      <td>2.612e-01</td>
      <td>1.412e-02</td>
    </tr>
    <tr>
      <td>2</td>
      <td>4.692e-02</td>
      <td>2.924e-03</td>
    </tr>
    <tr>
      <td>3</td>
      <td>9.100e-03</td>
      <td>5.511e-04</td>
    </tr>
    <tr>
      <td>4</td>
      <td>1.739e-03</td>
      <td>1.059e-04</td>
    </tr>
    <tr>
      <td>5</td>
      <td>3.334e-04</td>
      <td>2.028e-05</td>
    </tr>
    <tr>
      <td>6</td>
      <td>6.387e-05</td>
      <td>3.886e-06</td>
    </tr>
    <tr>
      <td>7</td>
      <td>1.224e-05</td>
      <td>7.445e-07</td>
    </tr>
    <tr>
      <td>8</td>
      <td>2.345e-06</td>
      <td>1.426e-07</td>
    </tr>
    <tr>
      <td>9</td>
      <td>4.492e-07</td>
      <td>2.733e-08</td>
    </tr>
  </tbody>
</table>
```

With a bubble

``` 
<table class="styled-table">
<thead>
    <tr style="text-align: right;">
      <th>Iteration</th>
      <th>Residual</th>
      <th>Variation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>1</td>
      <td>2.921e-01</td>
      <td>1.564e-02</td>
    </tr>
    <tr>
      <td>2</td>
      <td>5.613e-02</td>
      <td>3.525e-03</td>
    </tr>
    <tr>
      <td>3</td>
      <td>1.182e-02</td>
      <td>7.171e-04</td>
    </tr>
    <tr>
      <td>4</td>
      <td>2.442e-03</td>
      <td>1.492e-04</td>
    </tr>
    <tr>
      <td>5</td>
      <td>5.066e-04</td>
      <td>3.091e-05</td>
    </tr>
    <tr>
      <td>6</td>
      <td>1.050e-04</td>
      <td>6.409e-06</td>
    </tr>
    <tr>
      <td>7</td>
      <td>2.177e-05</td>
      <td>1.328e-06</td>
    </tr>
    <tr>
      <td>8</td>
      <td>4.512e-06</td>
      <td>2.754e-07</td>
    </tr>
    <tr>
      <td>9</td>
      <td>9.354e-07</td>
      <td>5.709e-08</td>
    </tr>
    <tr>
      <td>10</td>
      <td>1.939e-07</td>
      <td>1.183e-08</td>
    </tr>
    <tr>
      <td>11</td>
      <td>4.020e-08</td>
      <td>2.453e-09</td>
    </tr>
    <tr>
      <td>12</td>
      <td>8.333e-09</td>
      <td>5.085e-10</td>
    </tr>
    <tr>
      <td>13</td>
      <td>1.727e-09</td>
      <td>1.054e-10</td>
    </tr>
    <tr>
      <td>14</td>
      <td>3.581e-10</td>
      <td>2.185e-11</td>
    </tr>
    <tr>
      <td>15</td>
      <td>7.423e-11</td>
      <td>4.530e-12</td>
    </tr>
  </tbody>
</table>
```

With bubble Newton

```@raw html
<table class="styled-table">
  <thead>
    <tr style="text-align: right;">
      <th>Iteration</th>
      <th>Residual</th>
      <th>Variation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>1</td>
      <td>2.400e-02</td>
      <td>1.246e-02</td>
    </tr>
    <tr>
      <td>2</td>
      <td>8.490e-06</td>
      <td>2.507e-04</td>
    </tr>
    <tr>
      <td>3</td>
      <td>1.052e-12</td>
      <td>8.810e-08</td>
    </tr>
    <tr>
      <td>4</td>
      <td>1.102e-14</td>
      <td>1.081e-14</td>
    </tr>
  </tbody>
</table>
```


Check Poisson iterations
```bash
python3 -c "import convergence_study; convergence_study.plot_convergence_study_func()" ../Flower.jl/examples/convergence_Butler.yml poisson_00000032_00000001.h5 poisson_00000032_00000002.h5
```


```bash
python3 -c "import convergence_study; convergence_study.plot_convergence_study_func()" ../Flower.jl/examples/convergence_Butler.yml mesh_00000*
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="./assets/Butler_phi_1D.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential</figcaption>
</figure>
```

!!! todo "TODO plot" 
    x coord check 0...

!!! todo "TODO log plot "
 
plot iterations of Poisson
```bash
python3 -c "import convergence_study; convergence_study.plot_convergence_study_func()" ../Flower.jl/examples/convergence_Butler.yml poisson_00000032_00000001.h5 poisson_00000032_00000002.h5 poisson_00000032_00000003.h5 poisson_00000032_00000004.h5 poisson_00000032_00000005.h5
```

Plot schematics
```bash
python3 -c "import convergence_study; convergence_study.plot_schematics_func()" ../Flower.jl/examples/convergence_Butler_bubble_wall_Newton.yml mesh_00000*
```

```bash
python3 -c "import plot_flower; plot_flower.plot_all_fig_func()" ../Flower.jl/examples/convergence_Butler_bubble_wall_Newton.yml mesh_00000*
```
```bash
python3 -c "import plot_flower; plot_flower.plot_all_fig_func()" ../Flower.jl/examples/convergence_Butler_bubble_wall_Newton.yml mesh_00000128.h5 
```


```@raw html
<a href="https://github.com/flnt/Flower.jl/blob/electrolysis/examples/convergence_Butler.yml"> convergence_Butler.yml </a>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler/current_lines_32.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler/phi.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential</figcaption>
</figure>
```


```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler/poisson_iter_mesh32.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential: iterations for mesh 32</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler/poisson_iter_mesh64.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential: iterations</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler/poisson_iter_mesh128.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential: iterations</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler/poisson_iter_mesh256.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential: iterations</figcaption>
</figure>
```


## Butler-Volmer with bubble (no analytical solution)

```@raw html
Parameter file: <a href="https://github.com/flnt/Flower.jl/blob/electrolysis/examples/convergence_Butler_bubble.yml"> convergence_Butler_bubble.yml </a>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler_bubble/current_lines_32.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler_bubble/phi.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler_bubble/phi_through_bubble.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential through bubble</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler_bubble/phi_wall.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential at wall</figcaption>
</figure>
```

```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/examples/convergence_Butler.jl ../Flower.jl/examples/convergence_Butler_bubble.yml 
```

```bash
python3 -c "import plot_flower; plot_flower.plot_all_fig_func()" ../Flower.jl/test/butler_bubble.yml flower_00000001.h5
```


## Butler-Volmer with bubble at wall (no analytical solution)

```@raw html
<a href="https://github.com/flnt/Flower.jl/blob/electrolysis/examples/convergence_Butler_bubble_wall.yml"> convergence_Butler_bubble_wall.yml </a>
```


```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler_bubble_wall/current_lines_32.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler_bubble_wall/phi.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler_bubble_wall/phi_through_bubble.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential through bubble</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_Butler_bubble_wall/phi_wall.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential at wall</figcaption>
</figure>
```

```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/examples/convergence_Butler_bubble_wall.jl ../Flower.jl/examples/convergence_Butler_bubble_wall.yml 
```

```bash
python3 -c "import plot_flower; plot_flower.plot_all_fig_func()" ../Flower.jl/test/butler_bubble.yml flower_00000001.h5
```



## Diffusion

The following conditions have to be checked:
* check ``c_{\ce{H2}}  \geq c_{\ce{H2},0}``
* check ``c_{\ce{KOH}} \geq c_{\ce{KOH},0}``
* check ``0 < c_{\ce{H2O}} \leq c_{\ce{H2O},0}``

The discrete Laplacian operators are denoted L.

```@raw html
<a href="https://github.com/flnt/Flower.jl/blob/electrolysis/examples/convergence_diffusion.yml"> convergence_diffusion.yml </a>
```


!!! todo "Stopping criterion"

!!! danger "Bug or due to flow + production so conductivity increased hence the deformation of the current lines?"

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_diffusion/current_lines_32.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_diffusion/phi.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_diffusion/phi_through_bubble.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential through bubble</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_diffusion/phi_wall.svg" alt="Butler" title="Poisson">
    <figcaption>Electrical potential at wall</figcaption>
</figure>
```

!!! todo "Concentration"
    Not influenced by conductivity ? check conductivity

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_diffusion/concentration_H2_through_bubble.svg" alt="Butler" title="Poisson">
    <figcaption>Concentration through bubble</figcaption>
</figure>
```

```@raw html
<figure>
    <a name="Butler_slope"></a> 
    <img src="../../../convergence_diffusion/concentration_H2_wall.svg" alt="Butler" title="Poisson">
    <figcaption>Concentration at wall</figcaption>
</figure>
```

```@raw html
<!-- <body>
    <h1>Embed PDF Example</h1>
    <embed src="../../../convergence_diffusion/concentration_H2_wall.pdf" type="application/pdf" width="600" height="400">
</body>

<body>
    <h1>Object PDF Example</h1>
    <object data="../../../convergence_diffusion/concentration_H2_wall.pdf" type="application/pdf" width="600" height="400">
        <p>Your browser does not support PDFs. <a href="../../../convergence_diffusion/concentration_H2_wall.pdf">Download the PDF</a>.</p>
    </object>
</body> -->
```

```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/examples/convergence_Butler.jl ../Flower.jl/examples/convergence_diffusion.yml
```

```bash
python3 -c "import convergence_study; convergence_study.plot_convergence_study_func()" ../Flower.jl/examples/convergence_diffusion.yml mesh_0000*
```


```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/examples/convergence_Butler.jl ../Flower.jl/examples/convergence_diffusion.yml 
```

```bash
python3 -c "import convergence_study; convergence_study.plot_convergence_study_func()" ../Flower.jl/examples/convergence_diffusion.yml mesh_00000032.h5 
```

```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/examples/convergence_Butler.jl ../Flower.jl/examples/convergence_diffusion_constant_conductivity.yml
```


## Diffusion with a bubble


## Diffusion with a bubble at the wall

```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/examples/convergence_Butler.jl ../Flower.jl/examples/convergence_diffusion_constant_conductivity_bubble_wall.yml
```

## Diffusion and convection with a bubble at the wall

```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/examples/convergence_Butler.jl ../Flower.jl/examples/convergence_diffusion_constant_conductivity_bubble_wall.yml
```



## Radial flow (qualitative)

The velocity field is imposed with [`init_fields_multiple_levelsets!`](@ref) for "BC_uL" and "BC_vL". The gravity is deactivated.

```bash
python3 -c "import plot_flower; plot_flower.plot_all_fig_func()" ../Flower.jl/test/radial.yml flower_00000001.h5
```


## Poiseuille

!!! todo "output BC/results to html"
  https://github.com/JuliaPluto/HypertextLiteral.jl

### Flower

```julia
 BC_uL  left = Dirichlet 0.0 
 right = Dirichlet 0.0 
 bottom = Dirichlet 0.0 
 top = Neumann 0.0 
 LS 1 = Dirichlet 0.0 

 BC_vL  left = Dirichlet 0.0 
 right = Dirichlet 0.0 
 bottom = Dirichlet [4.1220703125000004e-5, 0.00011973632812500001, 0.00019301757812500002, 0.000261064453125, 0.000323876953125, 0.000381455078125, 0.00043379882812499994, 0.000480908203125, 0.000522783203125, 0.000559423828125, 0.0005908300781250002, 0.000617001953125, 0.0006379394531250002, 0.0006536425781250001, 0.0006641113281249999, 0.000669345703125, 0.0006693457031250001, 0.000664111328125, 0.0006536425781250001, 0.0006379394531250001, 0.000617001953125, 0.0005908300781249999, 0.0005594238281250001, 0.0005227832031250001, 0.0004809082031250002, 0.0004337988281250001, 0.00038145507812499985, 0.0003238769531250002, 0.000261064453125, 0.00019301757812499978, 0.000119736328125, 4.1220703125e-5] 
 top = Neumann 0.0 
 LS 1 = Dirichlet 0.0 

 BC_pL  left = Neumann 0.0 
 right = Neumann 0.0 
 bottom = Neumann 0.0 
 top = Dirichlet 0 

```

!!! danger "After one step"

```bash
[check_pressure_velocity_end]
min(u)  0.00e+00 min(v)  0.00e+00 min(p)  -0.00e+00
max(u)  0.00e+00 max(v)  6.69e-04 max(p)  3.34e-05
mean(u) 0.00e+00 mean(v) 2.22e-04 mean(p) 8.98e-06
min(grad_x)  -7.02e-02 min(grad_y)  -3.50e-01 
max(grad_x)  7.02e-02 max(grad_y)  -1.03e-09 
```



### Flower, top and bottom pressure imposed


```
[check_pressure_velocity_end]
min(u)  0.00e+00 min(v)  0.00e+00 min(p)  -4.52e-02
max(u)  0.00e+00 max(v)  6.69e-04 max(p)  -0.00e+00
mean(u) 0.00e+00 mean(v) 2.22e-04 mean(p) -1.20e-02
min(grad_x)  -2.76e-02 min(grad_y)  4.52e+02 
max(grad_x)  2.76e-02 max(grad_y)  4.52e+02 
```


### PmII, pressure imposed

```
[Before prediction]
min(u)  0.00e+00 min(v)  0.00e+00 min(p)  0.00e+00
max(u)  0.00e+00 max(v)  6.69e-04 max(p)  4.52e-02
mean(u) 0.00e+00 mean(v) 2.23e-04 mean(p) 1.20e-02
```


```@raw html
<table class="styled-table">
<thead>
<tr>
<th> Variable </th>
<th> left </th> <th> right </th> <th>bottom</th> <th>top </th> 
</tr>
</thead>
<tbody>
<tr> <td> u </td> <td> left = Dirichlet 0.0</td> <td> right = Dirichlet 0.0</td> <td>bottom = Dirichlet 0.0 </td> <td> top = Neumann 0.0 </td> </tr>
<tr> <td> v </td> <td> left = Dirichlet 0.0</td> <td> right = Dirichlet 0.0</td> <td>bottom = Neumann 0.0 </td> <td> top = Neumann 0.0 </td> </tr>
<tr> <td> p </td> <td> left = Neumann 0.0</td> <td> right = Neumann 0.0</td> <td>bottom = Dirichlet 0.045177296000000006 </td> <td> top = Dirichlet 0 </td> </tr>
</tbody>
</table>
```


!!! danger "With pressure BC for phi, a part of the increment of pressure"
    cf [`Brown et al. (2001)`](https://www.sciencedirect.com/science/article/pii/S0021999101967154)

    min(u)  0.00e+00 min(v)  0.00e+00 min(p)  -1.42e-06
    max(u)  0.00e+00 max(v)  6.69e-04 max(p)  4.52e-02
    mean(u) 0.00e+00 mean(v) 2.23e-04 mean(p) 1.33e-03
    min(grad_x)  -2.85e+04 min(grad_y)  -2.89e+04 
    max(grad_x)  1.04e-06 max(grad_y)  1.45e-02 


!!! info "With homogeneous Neumann for phi, a part of the increment of pressure"
    [check_pressure_velocity_end]
    min(u)  0.00e+00 min(v)  0.00e+00 min(p)  0.00e+00
    max(u)  0.00e+00 max(v)  6.69e-04 max(p)  4.52e-02
    mean(u) 0.00e+00 mean(v) 2.23e-04 mean(p) 1.20e-02
    min(grad_x)  -1.69e-13 min(grad_y)  -4.52e+02 
    max(grad_x)  1.69e-13 max(grad_y)  -4.52e+02 



### BC pressure imposed and BC for pressure increment

```@raw html
<table class="styled-table">
<thead>
<tr>
<th> Variable </th>
<th> left </th> <th> right </th> <th>bottom</th> <th>top </th> 
</tr>
</thead>
<tbody>
<tr> <td> u </td> <td>Dirichlet 0.0</td> <td>Dirichlet 0.0</td> <td>Dirichlet 0.0 </td> <td> Neumann 0.0 </td> </tr> 
 <tr> <td> v </td> <td>Dirichlet 0.0</td> <td>Dirichlet 0.0</td> <td>Neumann 0.0 </td> <td> Neumann 0.0 </td> </tr> 
<tr> <td> p </td> <td>Neumann 0.0</td> <td>Neumann 0.0</td> <td>Dirichlet 0.045177296000000006 </td> <td> Dirichlet 0 </td> </tr>
</tbody>
</table>
```

```@raw html
<figure>
    <a name="PoiseuillePmIIiposedpressureBCincrement"></a> 
    <img src="../../../convergence_Poiseuille_prediction_PmII_imposed_pressure_BC_increment/convergence_Poiseuille_errors.svg" alt="Poiseuille" title="Poiseuille">
    <figcaption>Errors on velocity</figcaption>
</figure>
```







!!! todo "compare"


Pressure-velocity coupling
* 1 : rh_v: add M * v
* prediction
* correction

!!! todo "BUG p"

!!! todo "BUG v"

We expect the following profile:

```bash
min(u)  0.00e+00 min(v)  0.00e+00 min(p)  0.00e+00
max(u)  0.00e+00 max(v)  6.69e-04 max(p)  4.45e-02
mean(u) 0.00e+00 mean(v) 6.38e-06 mean(p) 1.06e-02
```



```@docs
test_laplacian_pressure
```

We want to verify [`test_laplacian_pressure`](@ref)



Laplacian v, opC_p.iMy or else ? op_CvL.M 

```julia
M_pL = Diagonal(fzeros(grid))
iMx_pS = Diagonal(zeros((nx+1)*ny))
iMx_pL = Diagonal(zeros((nx+1)*ny))
iMy_pS = Diagonal(zeros(nx*(ny+1)))
iMy_pL = Diagonal(zeros(nx*(ny+1)))
```

```julia
M_vS = Diagonal(fzeros(grid_v)) #nx*(ny+1)
M_vL = Diagonal(fzeros(grid_v))
iMx_vS = Diagonal(zeros((grid_v.nx+1)*grid_v.ny))
iMx_vL = Diagonal(zeros((grid_v.nx+1)*grid_v.ny))
iMy_vS = Diagonal(zeros(grid_v.nx*(grid_v.ny+1)))
iMy_vL = Diagonal(zeros(grid_v.nx*(grid_v.ny+1)))
```

Based on:

```julia   
grav_y = g .* cos(β) .* opC_v.M * fones(grid_v)
```

```bash
python3 -c "import convergence_study; convergence_study.plot_convergence_study_func()" ../Flower.jl/examples/convergence_Poiseuille.yml mesh_00000*
```

```bash
python3 -c "import convergence_study; convergence_study.plot_errors_from_h5()" ../Flower.jl/examples/convergence_Poiseuille.yml convergence_Poiseuille.h5
```

```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/examples/convergence_Poiseuille.jl ../Flower.jl/examples/channel_Dirichlet_pressure.yml
```

```bash
python3 -c "import convergence_study; convergence_study.plot_convergence_study_func()" ../Flower.jl/examples/channel_Dirichlet_pressure.yml mesh_00000*
```


```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/examples/convergence_Poiseuille.jl ../Flower.jl/examples/channel_Dirichlet_pressure.yml
```

```bash
python3 -c "import convergence_study; convergence_study.plot_convergence_study_func()" ../Flower.jl/examples/channel_Dirichlet_pressure.yml mesh_00000*
```

```@raw html
<a href="https://github.com/flnt/Flower.jl/blob/electrolysis/test/channel_Dirichlet_pressure.yml"> See this test for the factor 4/3. </a>
```

To plot the results, you can use:
```bash
python3 -c "import plot_flower; plot_flower.plot_all_fig_func()" ../Flower.jl/test/channel_Dirichlet_pressure.yml flower_00000000.h5 flower_00000001.h5

python3 -c "import convergence_study; convergence_study.plot_errors_from_h5()" ../Flower.jl/test/channel_Dirichlet_pressure.yml convergence_Poiseuille.h5

```

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


Need to change opC_u.M


## Diffusion (old cf heat equation)
```math
\frac{\partial T}{\partial t} = LT
```

Here, ``L`` is the discrete linear operator. We couple the Cut Cell space discretization with a Crank-Nicolson time discretization, where ``τ`` is the time step, ``∆`` the uniform grid spacing and ``n`` the current iteration, resulting in the following discrete system  

```math
\frac{T^n-T^{n-1}}{\tau} = \frac{1}{2} \left[  LT^{n-1} +  LT^{n} \right]
```



which requires the solution of a linear system forming a pentadiagonal matrix. We validate the method in different stationary setups. 

!!!  todo "Not analytical reference?"
    A convergence study is carried out for these cases, where the reference solution is taken as the simulation with the highest number of points per dimension. 

In each case, the initial temperature field is set to zero and we impose a Dirichlet boundary condition at the interface. The ratio τ /∆2 = 0.5 is kept constant as we increase the number of points. 

Convergence study of the Cut Cell method coupled with a Crank-Nicolson scheme when solving the heat equation inside a stationary circle with a Dirichlet boundary condition TD = 1 imposed at the interface. The top gures show the position of the interface in red and the temperature eld at nal time tf = 0.03125 for N = 16, 32, 64, 128. The middle gure show the normalized error in temperature eld with respect to the reference solution taken for N = 256. The bottom gure shows the convergence rate of the method in mixed cells, full cells and in all cells.

### Inside circle (TODO)
1. A solid circle of radius ``R = 0.85``, initialized in a ``2 \times 2`` domain. The level set function is defined as  ``φ(x, y) = px2 + y2 − R``.  The Dirichlet TD = 1 boundary condition is imposed at the interface. We solve only for the phase inside of the circle until a nal time tf = 0.03125. The simulation is carried out for di erent resolutions N = 16, 32, 64, 128 corresponding to 4, 16, 64, 256 iterations, respectively. 

!!! todo "No analytical solution /Richardson extrapolation?"
    The reference solution is taken for N = 256.


 The results are summarized in Figure 3.8. As expected, the order of convergence of the error in full cells is close to 2 while the order of convergence in mixed cells is slightly less than 2. This drop in order in mixed cells is due to the assumption of a piece-wise linear interface approximation as well as the accumulation of errors of the bi-quadratic interpolation. The maximal errors are localized in cells where the wetted area is small (typically smaller than 5% of ∆2). Nevertheless, the global order of convergence in all cells is exactly 2.  

### Outside circle (TODO)
We initialize a solid circle of radius R = 0.75. The domain size are the same as well as the considered grid resolutions than in the previous. This time, we solve outside of the circle with a Dirichlet boundary condition TD = 1 at the interface and insulated boundary conditions at the domain boundary. The results, presented in Figure 3.9 are similar to case 1 with a slight drop in absolute error. This is due to the fact that there are less points per diameter than previously. This case validates the implementation of the Neumann boundary condition imposed at the domain boundaries.  

### Square (TODO)
3. In this third case, we initialize a square of area 1.6 × 1.6 in a 2 × 2 domain. The level set function is de ned as  ``φ(x, y) = max((x − 0.8), 9(x + 0.8), (y − 0.8), 9(y + 0.8)``.  We impose a Dirichlet boundary condition TD = 1 and solve inside of the square until the same nal time tf = 0.03125 with the same resolutions considered previously. In Figure 3.10, we can see that the maximal errors are located at the corners. The order of convergence for full cells is similar to the circle cases as well as for all cells. This case exhibits the robustness of the method when dealing with mesh aligned geometries as explained in Section 3.2.  

### Crystal (TODO)
4. Finally, in the last case, we consider a crystal in a 2 × 2 domain where the level set function is de ned as  φ(x, y) = px2 + y2 − R − 0.2 cos (6α) ,  where α is the angle of the interface with respect to the x axis and R = 0.7. 3.4. Validation on stationary geometries 45  At the interface, we impose the Gibbs-Thomson relation  TD = εκκ,  with κ the curvature and εκ = 0.01. The resulting temperature eld will now depend on the sign and amplitude of κ. We solve only for the phase inside of the circle until a nal time tf = 0.0078125. The simulation is performed for di erent resolutions N = 32, 64, 128, 256 corresponding to 4, 16, 64, 256 iterations respectively. The reference solution is taken for N = 512. In Figure 3.11, we can observe a drop in order of convergence for full cells with respect to the cases where TD was constant. This is explained by the accuracy of the curvature computation (Equation 1.8). The error is maximal in regions where the radius of curvature is large, where the interface is quasi aligned with the grid.  With these validation cases, we close the chapter on the Cut Cell method for di usive transport. In the next chapter, we describe the rest of the numerical steps of the two-phase Stefan problem.

" [`Fullana (2017)`](https://theses.hal.science/tel-04053531/)



## Advection

"The advection equation  ∂u  ∂t + (u · ∇) u = 0, (2.68)  is solved inside a cylinder of non-dimensionalized radius 0.5. The vector field u is initialized with an angular velocity ω0 = 1, representing the rotation of a rigid body, and slip boundary conditions at the wall, resulting in velocity components expressed as
" [`Rodriguez (2024)`](https://theses.fr/s384455)



## Growth of bubble away from the wall

!!! todo "Documentation divergence for phase change"
    can test by initializing c=0.16 at a planar bubble and c1>0.16 in bulk to verify that we have dc/dx*s


## Growth of bubble at the wall



## Test with dummy advection

```bash
julia +1.10.5 --project=../Flower.jl --threads=1 ../Flower.jl/examples/main_concise.jl levelset_Butler.yml
```

## Phase-change with diffusion
From [`Fleckenstein and bothe 2015`](https://www.sciencedirect.com/science/article/pii/S0021999115005306):
"
First, we compare a 1D simulation with FS2D to an analytically solvable test case. Similar validation cases have been previously used by, e.g., Welch and Wilson [40] for the validation of numerical methods for evaporation.

For this validation case, we consider mass transfer from a (pure) gas phase into a semi-infinite liquid phase. Initially, the interface between ``\Omega^G`` and ``\Omega^L`` is at position ``x = x_0``. The concentration of the single component inside ``\Omega^G`` remains constant throughout the simulation. Then, under the assumption that the transferred gas is highly diluted in the liquid phase, the equation for the species concentration ``c`` in the liquid phase can be written as:"

```math
\begin{aligned}
\partial_{t} c & = D \partial_{x}^{2} c, & & x > 0, t > 0 \\
c(0, x) & = 0, & & x > 0 \\
c(t, 0) & = c^G / H, & & t > 0
\end{aligned}
```

where the origin of the coordinate system moves with the interface position and H is the Henry coefficient. The time-dependent solution to this diffusion equation is given by [6]:

```math
c(t, x) = -\frac{c^G}{H} \operatorname{erf}\left(\frac{x}{2 \sqrt{D t}}\right) + \frac{c^G}{H}
```

With this species distribution in the liquid phase, we obtain:

```math
\frac{\partial c}{\partial \mathbf{n}_{x}} = -\frac{c^G}{H \sqrt{\pi D}} \frac{1}{\sqrt{t}}
```

the normal derivative of the transfer species at the interface. Taking this concentration derivative, we can finally compute the velocity of the interface as:

```math
V_{\Sigma} = -D \frac{c^G}{H \sqrt{\pi D}} \frac{1}{\sqrt{t}}
```

i.e., the traveled length ``l`` of the interface is:

```math
\ell = 2 \frac{c^G}{H} \sqrt{t D / \pi}
```


" To validate the numerical method, we use the position of the interface as given by Equation (89) and compare with the simulation results from FS3D. These simulations employ as diffusion coefficients ``D = 3.7 × 10^-7 m^2/s`` and ``D = 3.7 × 10^-8 m²/s`` and a Henry coefficient of ``H = 2``. As the analytical solution has been obtained for a semi-infinite domain, ample distance between the interface position and the right domain boundary (as in Fig. 4) has been chosen. 


!!! todo "The comparison of analytical and numerical solution 
  of the problem up to the final time ``t_{\text{fin}] = 10 s`` is presented in... "


```@raw html
<figure>
    <a name="phase_change"></a> 
    <img src="./assets/phase_change.svg" alt="Phase change" title="Phase change">
    <figcaption>"Phase change"</figcaption>
</figure>
```


"The minor, slowly growing deviation between numerical and analytical solution is due to the fact that during the motion of the interface, cells with very small phase fractions occur and in some cases, the mass flux into or from the small phase fraction is larger than the physically possible value, i.e., it would lead to a concentration either being negative or larger than the Henry law allows. In this case, the flux is limited to the physically maximum possible value which introduces a small systematic error. It is important to note that this problem is visible only in 1D computations, since the occurrence of several grid cells with interface and with only a very small fraction of one of the phases is very unlikely. We indeed did not observe this in any 3D computation."

```math
\llbracket \rho \mathbf{v} \otimes (\mathbf{v} - \mathbf{v}^{\Sigma}) \rrbracket \cdot \mathbf{n}_{\Sigma} = \tilde{m} \left( \mathbf{v} - \tilde{m} \left\llbracket \frac{1}{\rho} \right\rrbracket \mathbf{n}_{\Sigma}
```

---

**References:**

[40] S.W. Welch, J. Wilson, A volume of fluid based method for fluid flows with phase change, J. Comput. Phys. 160 (2) (2000) 662-682.

[6] J. Crunk, The Mathematics of Diffusion, Oxford University Press, 1975.
"


## TODO check H2O variations at electrode





## Variable coefficient

## Suggestion: use test case from electrostatics from Griffiths




## Plots
When the border values are plotted, they are shown at ``-dx/2``.
!!! todo 
Todo store centroids 