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



The following test cases are planned:
```@raw html
<div w3-include-html="./assets/test_table.html"></div> 
```


```@raw html
<table class="styled-table">

    <!-- <tr>
    
    <td> </td>
    <td> Steady </td>
    <td> Steady </td>
    <td> Steady </td>
    <td> Unsteady</td>
    <td> Unsteady</td>

    </tr> -->

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
            <td> <div> <img class="mid" src="./assets/Poisson_square.svg" width="75vh" > </div></td>
            <td><img src="./assets/Poisson_square_circle.svg" width="75vh" ></td>
            <!-- <td><img src="./assets/Poisson_square_circle_arc.svg" width="75vh" ></td> -->
            <td><a href="test.html#Poisson-equation-inside-a-square-with-circular-interface-at-wall"><img src="./assets/Poisson_square_circle_arc.svg" width="75vh" ></a></td>
            <td><img src="./assets/Poisson_square_circle_arrow.svg" width="75vh" ></td>
            <td><img src="./assets/Poisson_square_circle_arc_arrow.svg" width="75vh" ></td>
        </tr>
 
        <th>Poisson: $ \nabla \cdot (\kappa \nabla \phi)=0$ </th>
        <td>• • •</td>
        <td>  </td>
        <td> </td>
        <td> </td>
        <td> </td>
        </tr>

        <th>Diffusion</th>
        <td> </td>
        <td> </td>
        <td> </td>
        <td> </td>
        <td> </td>
        </tr>

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

## Manufactured solutions


See [`set_poisson`](@ref)

## Poisson equation inside a circle

See 2.4.4 in [`Rodriguez 2024`](https://theses.fr/s384455)

See [`(Rodriguez et al. 2024)`](https://link.springer.com/article/10.1007/s00707-024-04133-4)

``p (x, y) = cos (π^2 x y) sin (π^2 x y)``

### Axisymmetric case

``p(x, y) = x^2 + y^2``

### Poisson equation inside a square
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



See [`set_poisson_variable_coeff_SPD!`](@ref) in [Electrical potential](@ref) and [Poisson equation](@ref).

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

We study the relative errors in the discrete ``l_2`` and ``l_\infty`` norms in 2D:

```math
\begin{aligned}
   l_2&=\sqrt{\frac{\sum{S_i\left( p_i- p_{i}^e\right)}^2}{\sum{S_i\left( p_{i}^e\right)}^2}}
\end{aligned}
```

```math
\begin{aligned}
   l_\infty&=\frac{\max\lvert p_i- p_{i}^e\rvert}{\max( \lvert p_{i}^e \rvert )}
\end{aligned}
```


Recalling [Poisson equation](@ref)

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


In [`set_poisson`](@ref) system for `` + \nabla \cdot \nabla p = f`` is:

```math
\begin{cases}
+ \mathrm{div} (q^\omega, q^\gamma ) &= V f^\omega\\
I_a I p^\gamma + I_b \mathrm{div} (0, q^\omega) &= I g^\omega\\
q^\omega &= \mathrm{grad} ( p^\omega, p^\gamma ) \\
q^\gamma &= q^\omega
\end{cases}
```

```math
\begin{equation}
    \left [ \begin{array}{>{\centering\arraybackslash$} p{2.0cm} <{$} >{\centering\arraybackslash$} p{3.2cm} <{$}}
    -G ^ \top W ^ \dagger G & -G ^ \top W ^ \dagger H \\
    I _ b (-H ^ \top W ^ \dagger G) & I _ b (-H ^ \top W ^ \dagger H) + I _ a I _ \Gamma
    \end{array} \right ] \left [ \begin{array}{c}
    p ^ \omega \\
    p ^ \gamma
    \end{array} \right ] \simeq \left [ \begin{array}{c}
    V f ^ \omega \\
     I _ \Gamma g ^ \gamma
    \end{array} \right ],
\end{equation}
```

cf. :
```julia
A[end-nb+1:end,1:ni] = -b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
```

```julia
A[end-nb+1:end,1:ni] = b_b * (-1) * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
```
This corresponds to ``I _ b (-H ^ \top W ^ \dagger G)``


cf. :
```julia
A[end-nb+1:end,end-nb+1:end] = -pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .- χ_b * a1_b, 4.0)
```

```julia
A[end-nb+1:end,end-nb+1:end] = pad(b_b *(-1) * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .+ χ_b * a1_b, -4.0)
```
This corresponds to ``I _ b (-H ^ \top W ^ \dagger H) + I _ a I _ \Gamma``

```@docs
pad
```

```julia
A[sb,i*ni+1:(i+1)*ni] = -b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
```

##### At the interface
* Dirichlet __a1 = -1.0
               
* Neumann  __b = 1.0
* Robin __a1 = -1.0
        __a2 = 0.0
        __b = 1.0

In we have

In [`set_borders!`](@ref), we have:

Dirichlet: a1 = -1
Neumann: b =1
Robin: a1 = -1 b = 1
a0 : BC value 

!!! todo "Signs"
    * why a1 = -1
    * why -a0 


```julia
A[sb,sb] = -pad(
b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .- χ[iLS] * a1 .+
a2 * Diagonal(diag(fs_mat)), 4.0
)
```

```julia
A[sb,sb] = pad(
b * (-1) * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .+ χ[iLS] * a1 .-
a2 * Diagonal(diag(fs_mat)), -4.0
)
```


!!! todo
    ``a1 = -1``  

```julia
A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b)
```

Why ``+b``? and everywhere else ``-b`` ?



```julia
A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b)
# Boundary conditions for outer boundaries
A[end-nb+1:end,sb] = -b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
end

veci(rhs,grid,iLS+1) .= -χ[iLS] * vec(a0[iLS])
end

vecb(rhs,grid) .= -χ_b * vec(a0_b)
```


```julia
function set_poisson(
    bc_type, num, grid, a0, opC, opC_u, opC_v,
    A, L, bc_L, bc_L_b, BC,
    ls_advection)
    @unpack Bx, By, Hx, Hy, HxT, HyT, χ, M, iMx, iMy, Hx_b, Hy_b, HxT_b, HyT_b, iMx_b, iMy_b, iMx_bd, iMy_bd, χ_b = opC

    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny

    rhs = fnzeros(grid, num)

    a0_b = zeros(nb)
    _a1_b = zeros(nb)
    _b_b = zeros(nb)
    # Dirichlet: a1 = -1
    # Neumann: b =1
    # Robin: a1 = -1 b = 1
    # a0 : BC value 
    for iLS in 1:num.nLS
        set_borders!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    end
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    if ls_advection
        # Poisson equation
        A[1:ni,1:ni] = pad(L, -4.0)
        A[1:ni,end-nb+1:end] = bc_L_b

        # Boundary conditions for outer boundaries
        A[end-nb+1:end,1:ni] = -b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
        A[end-nb+1:end,end-nb+1:end] = -pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .- χ_b * a1_b, 4.0)
    end

    for iLS in 1:num.nLS
        if ls_advection
            if is_dirichlet(bc_type[iLS])
                __a1 = -1.0
                __a2 = 0.0
                __b = 0.0
            elseif is_neumann(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_robin(bc_type[iLS])
                __a1 = -1.0
                __a2 = 0.0
                __b = 1.0
            elseif is_fs(bc_type[iLS])
                __a1 = 0.0
                __a2 = 1.0
                __b = 0.0
            elseif is_wall_no_slip(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_navier(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_navier_cl(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            else
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            end
    
            _a1 = ones(grid) .* __a1
            a1 = Diagonal(vec(_a1))
            _a2 = ones(grid) .* __a2
            a2 = Diagonal(vec(_a2))
            _b = ones(grid) .* __b
            b = Diagonal(vec(_b))

            fs_mat = HxT[iLS] * Hx[iLS] .+ HyT[iLS] * Hy[iLS]

            sb = iLS*ni+1:(iLS+1)*ni
            
            # Poisson equation
            A[1:ni,sb] = bc_L[iLS]
            # Boundary conditions for inner boundaries
            A[sb,1:ni] = -b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
            # Contribution to Neumann BC from other boundaries
            for i in 1:num.nLS
                if i != iLS
                    A[sb,i*ni+1:(i+1)*ni] = -b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
                end
            end
            A[sb,sb] = -pad(
                b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .- χ[iLS] * a1 .+
                a2 * Diagonal(diag(fs_mat)), 4.0
            )
            A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b)
            # Boundary conditions for outer boundaries
            A[end-nb+1:end,sb] = -b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
        end

        veci(rhs,grid,iLS+1) .= -χ[iLS] * vec(a0[iLS])
    end

    vecb(rhs,grid) .= -χ_b * vec(a0_b)
    
    return rhs
end

```


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
<table border="1" class="dataframe">
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


## Diffusion

"

With the discrete Laplacian operators (one for each phase) now constructed, we can solve the heat equations on both domains  

```math
\frac{\partial T}{\partial t} = LT
```

on both domains. Here, ``L`` is the discrete linear operator. We couple the Cut Cell space discretization with a Crank-Nicolson time discretization, where ``τ`` is the time step, ``∆`` the uniform grid spacing and ``n`` the current iteration, resulting in the following discrete system  

```math
\frac{T^n-T^{n-1}}{\tau} = \frac{1}{2} \left[  LT^{n-1} +  LT^{n} \right]
```

!!! todo "Crank-Nicolson"
    why is there ``\frac{1}{\Delta ^2}`` in CN equation in Tomas's thesis?

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



## Suggestion: with a simple boundary condition similar to Butler-Volmer
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
    At ``x = L``, ``\Phi = 0`` so we have ``\Phi(x) = a(\frac{x}{L}-1)`` 
    
!!! todo ""
    with ``a>0`` because 

!!! todo
    At ``0^+``: ``\Phi(0^+) = -a`` and ``\frac{\partial\Phi}{\partial x }(0^+) = \frac{a}{L}``
    ``\frac{a}{L} = 2 i_0 \mathrm{sinh}(\frac{\alpha F}{RT}(-0.6-(-a)) )``
## Suggestion: use test case from electrostatics from Griffiths
