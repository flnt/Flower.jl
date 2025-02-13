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

# Documentation

This document contains extracts from the theses and articles listed in [References](@ref).

## Abbreviations
* BC: boundary conditions
* LS: levelset (distance to a fluid-fluid interface or fluid-solid interface i.e. a wall)
* p-grid: scalar grid, see [staggered grid](https://www.cfd-online.com/Wiki/Staggered_grid)
* u-grid: staggered grid (with regard to p-grid, in x-axis), used to discretize the x-component of the velocity 
* v-grid: staggered grid (with regard to p-grid, in y-axis), used to discretize the y-component of the velocity 
* ghost cells: cells defined outside of the computational domain so that the same stencils may be used near the wall (see [`allocate_ghost_matrices_2`](@ref), [`init_ghost_neumann_2`](@ref), [`IIOE_normal_indices_2!`](@ref); for MPI parallelization, ghost cells may also refer to cells outside the computational domain of the current processor, not implemented here)


## Algorithm, main function

```@docs
run_forward!
```
!!! todo "change order scalar transport electrical potential"


# Algorithm
* Initialize simulation parameters
* Allocations
* Initialisation of bulk and interfacial values
* Set matrices/operators
* Solve heat equation
* Poisson equation (electrical potential)
* Reevaluate the current i_current for boundary conditions 
* Scalar transport: update boundary conditions from electrical current
* Scalar transport: solve
* Phase change
* Advection of the levelset
* Pressure-velocity coupling (Navier-Stokes)



!!! info "Viewing the main parts of the algorithm in VScode"
    Comments starting with #region and ending with #endregion
    You can visualise the regions correponding to this algorithm more easily with tools like "Region marker" in VScode.

!!! danger "update LS"
    See [Updating the operator from Levelset](@ref)
    
## Parameters for Flower.jl

```@docs
Numerical
```

## Definitions

### Convention and orientation of the bubble interface
The signed distance is positive in the liquid and negative in the bubble. 

The normal defined by ``LS.\alpha`` is oriented towards the liquid, so we have to take the opposite to define the outward normal for manufactured solutions for example (so that it points towards the interior of the bubble when solving in the liquid). 


### Multiple levelsets

In Flower, multiple levelsets may be defined. The last levelset is the combination of all levelsets and is used to define the cell and its centroid.


!!! todo "TODO"
    Logic center refers to indices in the matrices.


"

In this work, a levelset function ``\phi`` \cite{Sethian1999} is defined on the computational domain ``\Omega`` to map the locus of one of its iso-levels (``\left \{ \left . x \in \Omega \right | \phi \left ( x, t \right ) = \phi _ 0 \right \}``) to an interface ``\Gamma \left ( t \right )`` that separates two non-overlapping domains, ``\Omega _ 1 \left ( t \right )`` and ``\Omega _ 2 \left ( t \right )``, each occupied by a different phase. The value ``\phi`` is defined as the signed distance to the interface,
```math
\begin{equation}
	\phi \left ( x, t \right ) = \left \{ \begin{aligned}
		-d \left ( x, \Gamma \left ( t \right ) \right ),& \ x \in \Omega _ 1 \left ( t \right ) \\
		0,& \ x \in \Gamma \left ( t \right ) \\
		d \left ( x, \Gamma \left ( t \right ) \right ),& \ x \in \Omega _ 2  \left ( t \right )
	\end{aligned} \right .,
\end{equation}
```
where ``d \left ( x, \Gamma \left ( t \right ) \right )`` denotes the minimal distance between the point ``x`` and the interface ``\Gamma \left ( t \right )``,

"

```@raw html
<figure>
    <a name="levelset_doc"></a> 
    <img src="./assets/levelset_doc.svg" alt="Levelset" title="Levelset defining an interface $\Gamma$ separating two phases $\Omega _ 1$ and $\Omega _ 2$">
    <figcaption>Levelset defining an interface $\Gamma$ separating two phases $\Omega _ 1$ and $\Omega _ 2$ </figcaption>
</figure>
```

### Geometric moments

The mesh is assumed to be rectilinear with ``n _ x`` points along ``x`` and ``n _ y`` along ``y``. For wet and partially wet cells numbered ``i`` along ``x`` and ``j`` along ``y``, the coordinates of the fluid center of mass are denoted ``x ^ \omega _ {i, j} \quad \mathrm{and} \quad y ^ \omega _ {i, j}`` and represent the components of two vectors of size ``n _ x n _ y`` (denoted ``x ^ \omega`` and ``y ^ \omega``). Likewise, the vectors ``x ^ \gamma`` and ``y ^ \gamma``, of size ``n _ x n _ y``, store the coordinates of the center of mass of the boundary of partially wet cells (referred to as boundary cells).

In the following, discrete values can be cell-, face- or node-centered, therefore a numbering convention must be adopted for the indices ``\left ( i, j \right )``. Fig.~\ref{fig:numbering} showcases the convention used in this work for both cell and face-centered quantities (the ``\left ( i, j \right )`` node-centered quantities is located at the bottom-left of the cell).

```@raw html
<figure>
    <a name="numbering_doc"></a> 
    <img src="./assets/numbering_doc.svg" alt="Numbering convention" title="Numbering convention">
    <figcaption> Numbering convention </figcaption>
</figure>
```

"
Cut-cell methods are firmly grounded in the Finite Volume Method, which defines the primary discrete variables as cell-wise averages over mesh elements. The design of the Finite Volume operators is then based on the application of the Divergence theorem. For example, given a scalar field ``p``,  this theorem states that in a Cartesian coordinate system, the ``x`` component of the gradient ``q \equiv \nabla p`` averaged over a cell ``\Omega`` may be computed as
\begin{equation}
\Omega q _ x = \int _ \Omega \frac{\partial p}{\partial x} V = \oint _ {\partial \Omega} p e _ x \cdot  S
\label{eq:Stokes}
\end{equation}
where ``S`` denotes the outward-pointing surface element and ``e_x`` the unit vector along the ``x`` direction.

For the sake of presentation, the case displayed in Fig.~\ref{fig:moments} is considered where $\Omega \equiv \mathcal{V} _ {i,j}$ consists of the intersection of a phase domain and a computational cell (a right hexahedron). The contour $\partial \Omega$ then consists of the union of the four planar faces $\mathcal{A} _ {x, i, j}$, $\mathcal{A} _ {x, i+1, j}$, $\mathcal{A} _ {y, i, j}$ and $\mathcal{A} _ {y, i, j+1}$ as well as the boundary surface $\Gamma _ {i, j}$. A piece-wise linear approximation of $\Gamma _ {i, j}$, denoted $\widetilde{\Gamma} _ {i, j}$ and of unit normal $\left ( n _ {x, i, j}, n _ {y, i, j} \right )$, can be defined by applying Eq.~\eqref{eq:Stokes} to $\mathcal{V} _ {i, j}$ with $p = 1$, yielding
```math
\int _ {\widetilde{\Omega}} \frac{\partial 1}{\partial x} \mathrm{d} V = \mathcal{A} _ {x, i+1, j} - \mathcal{A} _ {x, i, j} + n _ {x, i, j} \widetilde{\Gamma} _ {i, j} = 0
```
which highlights the existence of a fundamental relation
\begin{equation}
    \mathcal{A} _ {x, i+1, j}  - \mathcal{A} _ {x, i, j}  = -n _ {x, i, j} \widetilde{\Gamma} _ {i, j}
    \label{eq:approximate}
\end{equation}
sometimes referred to as a Geometric Surface Conservation Law. The same can be done in the $y$-direction in order to obtain a relation between $\mathcal{A} _ {y, i, j+1}, \mathcal{A} _ {y, i, j}, n _ {y, i, j}$ and $\widetilde{\Gamma} _ {i, j}$. Note that to simplify the exposition, the notations ($\mathcal{A} _ {x, i, j}$, $\mathcal{V} _ {i, j}$...) are used to denote both a region in space and its measure (length or area, in 2D).
"


### Domains
``\Gamma`` is the interface and ``\omega`` is the wet volume of the cell without the interface: `` \omega = \Omega \setminus \Gamma``. 

### Centroids

In the code, the cell centroid for the p-grid (scalar grid, here noted gp) corresponding to the first levelset (``LS[1]``) is ``x^\omega = (x_{\text{centroid}},y_{\text{centroid}})``:
```julia
x_centroid = gp.x .+ getproperty.(gp.LS[1].geoL.centroid, :x) .* gp.dx
y_centroid = gp.y .+ getproperty.(gp.LS[1].geoL.centroid, :y) .* gp.dy
```

For the v-grid:

```julia
x_centroid = gv.x .+ getproperty.(gv.LS[1].geoL.centroid, :x) .* gv.dx
y_centroid = gv.y .+ getproperty.(gv.LS[1].geoL.centroid, :y) .* gv.dy
```

The interface centroid ``x^\gamma`` is defined as "the mid point of the segment crossing the cell which will be used in the computation of the Stefan condition." [`(Fullana 2017)`](https://theses.hal.science/tel-04053531/).

```julia
x_bc = gp.x .+ getproperty.(gp.LS[1].mid_point, :x) .* gp.dx
y_bc = gp.y .+ getproperty.(gp.LS[1].mid_point, :y) .* gp.dy
```

!!! info
    ``(gp.x, gp.y)`` is the cell center, not the centroid. In the current post-processings the cell center is used if not said otherwise.

!!! info
    The divergence ``\mathrm{div}(0,\mathrm{grad})`` is a bulk field, it is located at the cell centroids. 



### H 

The boundary field is always located at the borders of the p-grid. So for scalars, the distance is always half a cell in both directions because the cell centroid is in the cell center. For the u-grid, since it's shifted half a cell in x, then from the wet cell centroid to the boundary there's a quarter of a cell at the left and right borders. And the same for the v-grid, it's displaced half a cell in y, so from the wet cell centroid to the boundary there's a quarter of a cell of distance.


!!! todo "TODO"
    Without interface, with a constant mesh spacing in every direction, non-dimensionalized by the mesh spacing H should be equal to:

| Border  | p-grid   | u-grid   | v-grid   |
|:--------|:--------:|:--------:|:--------:|
| left   | 0.5    | 0.25   | 0.5    |
| bottom | 0.5    | 0.5    | 0.25   |
| right  | 0.5    | 0.25   | 0.5    |
| top    | 0.5    | 0.5    | 0.25   |


### Capacities
The face capacities for the liquid phase of levelset n° iLS for ``grid`` are stored in grid.LS[iLS].geoL.cap[II,1:4]:

```julia
grid.LS[iLS].geoL.cap[II,1] # left, also known as A1 in the code
``` 
```julia
grid.LS[iLS].geoL.cap[II,2] # bottom, also known as A2 in the code
```

```julia
grid.LS[iLS].geoL.cap[II,3] # right, also known as A3 in the code
```

```julia
grid.LS[iLS].geoL.cap[II,4] # top, also known as A4 in the code
```

The cell-centered volume is:
```julia 
grid.LS[iLS].geoL.cap[II,5] # cell-centered volume
``` 
The cell-centered heights are:

```julia 
grid.LS[iLS].geoL.cap[II,6] # cell-centered, also known as B1
``` 

```julia 
grid.LS[iLS].geoL.cap[II,7] # B2
``` 
Others:
```julia 
grid.LS[iLS].geoL.cap[II,8] # W1
``` 
See [`set_cutcell_matrices!`](@ref)
```julia 
grid.LS[iLS].geoL.cap[II,9] # W2
``` 

```julia 
grid.LS[iLS].geoL.cap[II,10] # W3
``` 

```julia 
grid.LS[iLS].geoL.cap[II,10] # W4
``` 


!!! info 
    In the code, B1 B2 are placed at the cell centroid, not the cell center.



``iMx`` is the inverse of a volume. The border cells for u and v have half the normal cell size because the borders of the domain pass through the cell centers, not the cell faces like they do for p.






## Example for 1D (equations)

This section explains the equations and their implementation for a simple geometry with one levelset describing a square bubble and a wall at the limits of the computational domain. 


### Operators in a cell of the mesh
Let us consider a square bubble in contact with a wall. For simplicity, the faces of a regular cartesian p-grid are displayed in black. 
The wall (in red) and bubble (in green) interfaces are shifted so that two boundary conditions can be imposed simultaneously ( see [Simultaneous boundary conditions](@ref)).

```@raw html
<figure>
    <img src="./assets/mesh_square_bubble.svg" alt="Square bubble in contact with wall" title="Square bubble in contact with wall">
    <figcaption>Square bubble in contact with wall</figcaption>
</figure>
```

Each boundary condition is solved in a separate cell.

```@raw html
<figure>
    <img src="./assets/mesh_square_bubble_2.svg" alt="Square bubble in contact with wall" title="Square bubble in contact with wall">
    <figcaption>Square bubble in contact with wall</figcaption>
</figure>
```


```@raw html
The scalar control volume $\mathcal{V}$, $\mathcal{B}$ capacities (height of wetted part in x and y directions) are associated with the centroid of the current phase (e.g. liquid): $x^\omega$ <span class="hover_img"> <a href="#centroidcontrolvolume"> (see figure) <span><img src="./assets/Vt_mesh_square_bubble_wall_1.svg" alt="image" height="800" /></span></a></span>.
```

The centroid is required even if it does not intervene in the discretization. 
The [`LibGEOS.jl`](https://github.com/JuliaGeo/LibGEOS.jl) library, based on [`GEOS`](https://libgeos.org/), is used for the geometrical computations, the VOFI library could also be used: [`(bna et al., 2016)`](https://doi:10.1016/j.cpc.2015.10.026).


The scalar control volume ``\mathcal{W}``, ``\mathcal{A}`` capacities (height of wetted part in x and y directions) are associated with the center of the faces of the cells. 

The first kind capacities are ``x^\omega`` (cell), V (cell), ``A^x`` (face), the second kind capacities are ``B^x`` (cell) and ``\mathcal{W}^x`` (face).

The interface of the left wall is shifted at ``o^-`` and the bubble interface is also shifted so that two different boundary conditions (BC) may be imposed simultaneously.


When a volume is null, the associated capacity ``\mathcal{B}`` must be null as well. 


#### At wall

The scalar control volume of cell ``(i=1,j)`` and its associated centroid are represented below:
```@raw html
<figure>
    <a name="centroidcontrolvolume"></a> 
    <img src="./assets/V_mesh_square_bubble_wall_1.svg" alt="Centroid control volumes at i = 1" title="Centroid control volumes at i = 1">
    <figcaption>Centroid control volumes at i = 1</figcaption>
</figure>
```

The control volume ``\mathcal{W}_{x,i=1}`` associated with the face is represented below:

```@raw html
<figure>
    <a name="facecontrolvolume"></a> 
    <img src="./assets/Wx_mesh_square_bubble_wall_2.svg" alt="Face control volume " title="Face control volume">
    <figcaption>Face control volume</figcaption>
</figure>
```



The values of the capacities are detailed here:

```@raw html
<figure>
    <a name="capacities"></a> 
    <img src="./assets/Vt_mesh_square_bubble_wall_3.svg" alt="Centroid control volumes at i=0 and i = 1" title="Centroid control volumes at i=0 and i = 1">
    <figcaption>Centroid control volumes at i=0 and i = 1</figcaption>
</figure>
```





```@raw html
\begin{array}{cccccccc}
	\hline
	Variable & Face ~-\frac 12 & Bulk ~ i = 0 & \color{flred}{0^-} & Face ~\frac 12 & Bulk ~ i = 1 & Face ~\frac 32 &  Bulk ~ i = 2 \\
    \hline
	V &  & 0 & \color{flred}{|} & & h_xh_y & & h_xh_y\\
	A_x & 0 &  &\color{flred}{|} &h_y &  & h_y & \\
	B_x &  & 0 &\color{flred}{|} & & h_y & & h_y \\
	W_x & 0 &  &\color{flred}{|} &\frac{h_xh_y}{2} &  & h_xh_y &\\
    \hline
\end{array}
```


The location of the centroid of ``\mathcal{W}_x`` is not used in the current discretization.


#### Gradient

##### x-component

At cell ``(i,j)``, with ``\Gamma`` the interface and ``c`` the scalar, the x-component of the cell-averaged gradient is:

```math
\begin{aligned}
{\mathrm{grad}(c^ω, c^{γ})}_{i +\frac 12} &= \mathcal{W}_{i +\frac 12}^\dagger (B_{i+1}c_{i+1}^\omega-B_{i}c_{i}^\omega \\
&+ (A_{i+\frac 12} - B_{i+1})c_{i+1}^\gamma + (B_{i}-A_{i+\frac 12})c_{i}^\gamma)
\end{aligned}
```

```@raw html
In the case of 
<span class="hover_img"> <a href="#capacities">this figure <span><img src="./assets/Vt_mesh_square_bubble_wall_3.svg" alt="image" height="800" /></span></a></span>, the above expression is the same as with finite differences:
```

```math
\begin{aligned}
{\mathrm{grad}(c^ω, c^{γ})}_{0 +\frac 12}  &= \mathcal{W}_{\frac 12}^\dagger (B_{1}c_{1}^\omega-B_{0}c_{0}^\omega \\
&+ (A_{\frac 12} - B_{1})c_{1}^\gamma + (B_{0}-A_{\frac 12})c_{0}^\gamma) \\
&= \frac{1}{\frac{h_x}{2}h_y} (h_yc_{1}^\omega- 0c_{0}^\omega \\
&+ (h_y - h_y)c_{1}^\gamma + (0-h_y)c_{0}^\gamma) \\
&=\frac{c_1^\omega-c_0^\gamma}{\frac{h_x}{2}}
\end{aligned}
```

Away from the wall and without interface, the formula is:

```math
\begin{aligned}
{\mathrm{grad}(c^ω, c^{γ})}_{i +\frac 12} =\frac{c_{i+1}^\omega-c_i^\omega}{h_x}
\end{aligned}
```

##### y-component


The y-component is given below:

```math
\begin{aligned}
{\mathrm{grad}(c^ω, c^{γ})}_{i +\frac 12} &= \mathcal{W}_{i +\frac 12}^\dagger (B_{i+1}^yc_{j+1}^\omega-B_{i}^yc_{j}^\omega \\
&+ (A_{j+\frac 12} - B_{j+1})c_{i+1}^\gamma + (B_{j}-A_{j+\frac 12})c_{j}^\gamma)
\end{aligned}
```

Since ``\Gamma`` is shifted upwards in the y-direction ``A_{\frac 12}=0``:
```math
\begin{aligned}
{\mathrm{grad}(c^ω, c^{γ})}_{0 +\frac 12}  &= \mathcal{W}_{\frac 12}^\dagger (B_{1}c_{1}^\omega-B_{0}c_{0}^\omega \\
&+ (A_{\frac 12} - B_{1})c_{1}^\gamma + (B_{0}-A_{\frac 12})c_{0}^\gamma) \\
&= \frac{1}{\frac{h_x}{2}h_y} (h_x c_{1}^\omega- 0c_{0}^\omega \\
&+ (0 - h_x)c_{1}^\gamma + (0-0)c_{0}^\gamma) \\
&=\frac{c_i^\omega-c_0^\gamma}{\frac{h_y}{2}}
\end{aligned}
```


##### Right wall
For the right wall, with ``nx`` the number of cells in the x direction:

!!! todo "TODO"
    finish
```math
\begin{aligned}
{\mathrm{grad}(c^ω, c^{γ})}_{nx +\frac 12}  &= \mathcal{W}_{\frac 12}^\dagger (B_{nx+1}c_{nx+1}^\omega-B_{nx}c_{nx}^\omega \\
&+ (A_{nx+\frac 12} - B_{nx+1})c_{1}^\gamma + (B_{nx}-A_{nx+\frac 12})c_{nx}^\gamma) \\
&= \frac{1}{\frac{h_x}{2}h_y} (0 c_{nx+1}^\omega- h_yc_{nx}^\omega \\
&+ (h_y - 0)c_{nx+1}^\gamma + (h_y-h_y)c_{nx}^\gamma) \\
&=\frac{-c_1^\omega+c_{nx+1}^\gamma}{\frac{h_x}{2}}
\end{aligned}
```

##### Interface below

!!! todo "TODO"
    finish and check

Since ``\Gamma`` is shifted downwards in the y-direction ``A_{\frac 12}=0``:
```math
\begin{aligned}
{\mathrm{grad}(c^ω, c^{γ})}_{}  &= \mathcal{W}_{\frac 12}^\dagger (B_{}c_{1}^\omega-B_{0}c_{0}^\omega \\
&+ (A_{\frac 12} - B_{1})c_{1}^\gamma + (B_{0}-A_{\frac 12})c_{0}^\gamma) \\
&= \frac{1}{\frac{h_x}{2}h_y} (0 c_{1}^\omega- h_x c_{0}^\omega \\
&+ (0 - 0)c_{1}^\gamma + (h_x-0)c_{0}^\gamma) \\
&=\frac{-c_i^\omega+c_0^\gamma}{\frac{h_y}{2}}
\end{aligned}
```


#### Divergence

At cell ``(i,j)``, with ``\Gamma`` the interface

```@raw html
<a name="divergencetag"></a> 
```


```math
\begin{equation}
\label{divergence}
\begin{aligned}
{\mathrm{div}(q^ω, q^{γ})}_i &= \color{florange}{A^x_{i+\frac{1}{2}} q_{i+\frac{1}{2}}^\omega-A^x_{i-\frac{1}{2}} q_{i-\frac{1}{2}}^\omega} \\
&+\color{flblue}{(B_i^x -A^x_{i+\frac{1}{2}}) q_{i+\frac{1}{2}}^\gamma +(A^x_{i-\frac{1}{2}}-B_i^x ) q_{i-\frac{1}{2}}^\gamma } \\
&+\color{flgreen}{A^y_{j+\frac{1}{2}} q_{j+\frac{1}{2}}^\omega -A^y_{j-\frac{1}{2}} q_{j-\frac{1}{2}}^\omega } \\
&+\color{flblue}{(B_j^y -A^y_{j+\frac{1}{2}}) q_{j+\frac{1}{2}}^\gamma +(A^y_{j-\frac{1}{2}}-B_j^y ) q_{j-\frac{1}{2}}^\gamma } \\
\end{aligned}
\end{equation}
```


```math
\begin{aligned}
\int_\Omega \nabla \cdot q dV &= \oint q \cdot n dS \\
& = \color{florange}{\oint_{i+\frac 12} q \cdot (+e_x) dS} + \color{florange}{\oint_{i-\frac 12} q \cdot (-e_x) dS} \\
& + \color{flgreen}{\oint_{j+\frac 12} q \cdot (+e_y) dS} + \color{flgreen}{\oint_{j-\frac 12} q \cdot (-e_y) dS} \\
& + \color{flblue}{\oint_\Gamma q \cdot n dS} 
\end{aligned}
```



```@raw html
In the case of <span class="hover_img">
     <a href="#capacities">this figure <span><img src="./assets/Vt_mesh_square_bubble_wall_3.svg" alt="image" height="800" /></span></a>
</span>, the x-contribution from the wall (from the above expression) is:
```

```@raw html
<a name="tagdivleftwall"></a> 
```
```math
\color{flblue}{(0-h_y) (\frac{c_i^\omega-c_0^\gamma}{\frac{h_x}{2}}) + (0-0) q_{i-\frac{1}{2}}^\gamma }
```
Since 

```@raw html
<a name="tagdivbubbleinterface"></a> 
```
The y-contribution for the green bubble is (by denoting j=1 the index of the cell cut by the interface, in the phase under consideration, i.e. above the bubble here):
 ```math
\color{flblue}{(h_x-h_x) q_{\frac 32}^y + (0-h_x) q_{\frac 12}^y=-h_x (\frac{c_1^\omega-c_1^\gamma}{\frac{h_y}{2}}) }
```
The two former equations give an example of discretization of the boundary conditions at the left wall and at the interface.



##### Right wall

At nx+1

```math
\color{flblue}{(0-0) (\frac{c_i^\omega-c_0^\gamma}{\frac{h_x}{2}}) + (h_y-0) q_{i-\frac{1}{2}}^\gamma }
```

```math
\color{flblue}{ h_y q_{i-\frac{1}{2}}^\gamma \frac{-c_1^\omega+c_{nx+1}^\gamma}{\frac{h_x}{2}}}
```

Since 

```@raw html
<a name="tagdivbubbleinterface"></a> 
```
The y-contribution for the green bubble (under the bubble) is:
 ```math
\color{flblue}{(h_x-0) q_{i+\frac 12}^y + (h_x-h_x) q_{i-\frac 12}^y=h_x (\frac{-c_1^\omega+c_1^\gamma}{\frac{h_y}{2}}) }
```

If ``h_x=h_y``, we have the following coefficients on the rows of the interfacial values: +2 for interfacial, -2 for bulk

!!! danger "Check sign"


##### Coefficients in a simple configuration
```@raw html
<a name="Coefficients-in-a-simple-configuration"></a> 
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

##### Example away from interfaces

Away from the wall, without an interface, with a regular mesh of constant spacings ``h_x`` and ``h_y`` is:


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


##### Approximations inside a cell 

The following hypothesis is made, ``q^\omega=q^\gamma=q``, which means that the cell-averaged gradient (corresponding to the control volume ..., associated with the centroid ...), is associated at the middle of the part of interface, as illustrated by the black arrows below.


```@raw html
<figure>
    <img src="./assets/Wx_mesh_square_bubble_wall_4.svg" alt="Face control volume " title="Face control volume">
    <figcaption>Face control volume</figcaption>
</figure>
```

## Poisson equation

!!! todo "Old implementation"
    In [`set_poisson`](@ref) the system is `` + \nabla \cdot \nabla p = f`` :

In [`solve_poisson`](@ref) system is `` + \nabla \cdot \nabla p = f`` :




```math
\begin{cases}
+ \mathrm{div} (q^\omega, q^\gamma ) &= V f^\omega\\
I_a I p^\gamma + I_b \mathrm{div} (0, q^\omega) &= I g^\omega\\
q^\omega &= \mathrm{grad} ( p^\omega, p^\gamma ) \\
q^\gamma &= q^\omega
\end{cases}
```


Second line: rechecking

!!! danger "Question"
    since div in bulk and BC, and BC part of div is of opposite sign, why not opposite signs in matrix ?

full divergence

```@raw html
<a name="tagfulldivergence"></a> 
```
```math
\begin{equation}
\mathrm{div}(q^ω, q^{γ}  ) = − (G^⊤ + H^{Γ1,⊤} + H^{Γ2,⊤}) q^ω + H^{Γ1,⊤}q^{γ1} + H^{Γ2,⊤}q^{γ2}
\end{equation}
```

```@raw html
<a name="grad_concise"></a> 
```
```math
\mathrm{grad}(q^ω, q^{γ}  ) = W ^ \dagger \left ( G p ^ \omega + H p ^ \gamma \right )
```

```math
\begin{equation}
    \left [ \begin{array}{>{\centering\arraybackslash$} p{2.0cm} <{$} >{\centering\arraybackslash$} p{3.2cm} <{$}}
    -G ^ \top W ^ \dagger G & -G ^ \top W ^ \dagger H \\
    I _ b (H ^ \top W ^ \dagger G) & I _ b (H ^ \top W ^ \dagger H) + I _ a I _ \Gamma
    \end{array} \right ] \left [ \begin{array}{c}
    p ^ \omega \\
    p ^ \gamma
    \end{array} \right ] \simeq \left [ \begin{array}{c}
    V f ^ \omega \\
     I _ \Gamma g ^ \gamma
    \end{array} \right ],
\end{equation}
```

### Interface length

The interface length is defined as:

``\chi = \sqrt{{(A_{\frac 12}^x-B_0^x)}^2+{(B_0^x-A_{-\frac 12}^x)}^2 + {(A_{\frac 12}^y-B_0^y)}^2+{(B_0^y-A_{-\frac 12}^y)}^2}``

In [`Rodriguez 2024`](https://theses.fr/s384455), it is defined as:

``Y_\Gamma = \mathrm{diag}(|H^⊤1|)``

With ``G=\begin{bmatrix} D^{x}_- B^x \\ D^{y}_- B^y \end{bmatrix}`` and ``H=\begin{bmatrix}A^x D^x_- -D^x_- B^x_-\\A^y D^y_- -D^y_- B^y_- \end{bmatrix}``



In [`set_cutcell_matrices!`](@ref) and [`scalar_transport!`](@ref), the interface length ``\chi`` is thus computed:

```julia
χx = (grid.LS[iLS].geoL.dcap[:,:,3] .- grid.LS[iLS].geoL.dcap[:,:,1]) .^ 2
χy = (grid.LS[iLS].geoL.dcap[:,:,4] .- grid.LS[iLS].geoL.dcap[:,:,2]) .^ 2
op.χ[iLS].diag .= sqrt.(vec(χx .+ χy))
```

### System

#### For the bulk: 1:ni
With the [`pad`](@ref) function:

```@docs
pad
```

```julia
A[1:ni,1:ni] = pad(L, -4.0)
A[1:ni,end-nb+1:end] = bc_L_b
```
The function [`laplacian`](@ref) returns

```
L = BxT * iMx * Bx .+ ByT * iMy * By
bc_L[iLS]= BxT * iMx * Hx[iLS] .+ ByT * iMy * Hy[iLS]
bc_L_b = (BxT * iMx_b * Hx_b .+ ByT * iMy_b * Hy_b)
```
 
!!! info "Correpondence article-implementation"
    cf. [`divergence_B!`](@ref)
    ni: bulk so BxT →-G^T , Bx→G in the article 
    Hx → H
    Hx_b →H

So ``A[1:ni,1:ni] = pad(L, -4.0)`` corresponds to ``\mathrm{div}``: ``-G ^ \top W ^ \dagger G ``

```julia
A[1:ni,end-nb+1:end] = bc_L_b
``` 
corresponds to `` -G ^ \top W ^ \dagger H `` (``Hx_b`` is for the wall)

#### For the wall: A[end-nb+1:end,:]
```julia
A[end-nb+1:end,1:ni] = b_b * (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By)
```
corresponds to 

``
I _ b (H ^ \top W ^ \dagger G) 
``



in cf. [`set_border_matrices!`](@ref):
cf. [`bc_matrix_borders!`](@ref) : - sign to the left (cf the outward normal), +sigh to the right
and with

```julia 
mat_assign_T!(opC_p.HxT_b, sparse(opC_p.Hx_b'))
```
``Hx_b`` is the transpose of ``HxT_b``



cf. :
```julia
A[end-nb+1:end,end-nb+1:end] = pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .+ χ_b * a1_b, -4.0)
```
corresponds to
```math
I _ b (H ^ \top W ^ \dagger H) + I _ a I _ \Gamma
```

!!! todo "-4 or 4"
    -4 (convention, value not important ?)





```julia
A[1:ni,sb] = bc_L[iLS]
```
corresponds `` -G ^ \top W ^ \dagger H `` (``Hx[iLS]`` is for the wall in [`laplacian`](@ref))



#### For the interface (other than wall)
```julia
# Boundary conditions for inner boundaries
A[sb,1:ni] = b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
# Contribution to Neumann BC from other boundaries
for i in 1:num.nLS
    if i != iLS
        A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
    end
end
A[sb,sb] = pad(
    b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .+ χ[iLS] * a1 .-
    a2 * Diagonal(diag(fs_mat)), -4.0
)
#TODO was +b... here in set_poisson 
A[sb,end-nb+1:end] = b * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b) #why was it different in set_poisson ?
# Boundary conditions for outer boundaries
A[end-nb+1:end,sb] = b_b * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
end

veci(rhs,grid,iLS+1) .= +χ[iLS] * vec(a0[iLS]) #was - in set_poisson
end

vecb(rhs,grid) .= +χ_b * vec(a0_b) #was - in set_poisson
    
```

```julia
A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
```

!!! todo "Sign free-surface, check"


In [`solve_poisson`](@ref):

```julia
function solve_poisson(
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
    for iLS in 1:num.nLS
        set_borders_poisson!(grid, grid.LS[iLS].cl, grid.LS[iLS].u, a0_b, _a1_b, _b_b, BC, num.n_ext_cl)
    end
    a1_b = Diagonal(vec(_a1_b))
    b_b = Diagonal(vec(_b_b))

    if ls_advection
        # Poisson equation
        A[1:ni,1:ni] = pad(L, -4.0)
        A[1:ni,end-nb+1:end] = bc_L_b

        # Boundary conditions for outer boundaries
        A[end-nb+1:end,1:ni] = b_b * ( (HxT_b * iMx_b' * Bx .+ HyT_b * iMy_b' * By) )
        A[end-nb+1:end,end-nb+1:end] = pad(b_b * (HxT_b * iMx_bd * Hx_b .+ HyT_b * iMy_bd * Hy_b) .+ χ_b * a1_b,- 4.0)
    end

    for iLS in 1:num.nLS
        if ls_advection
            if is_dirichlet(bc_type[iLS])
                __a1 = 1.0
                __a2 = 0.0
                __b = 0.0
            elseif is_neumann(bc_type[iLS])
                __a1 = 0.0
                __a2 = 0.0
                __b = 1.0
            elseif is_robin(bc_type[iLS])
                __a1 = 1.0
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
            A[sb,1:ni] = b * (-1) * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)
            # Contribution to Neumann BC from other boundaries
            for i in 1:num.nLS
                if i != iLS
                    A[sb,i*ni+1:(i+1)*ni] = b * (-1) * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i])
                end
            end
            A[sb,sb] = pad(
                b * (-1) *(HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .+ χ[iLS] * a1 .-
                a2 * Diagonal(diag(fs_mat)), -4.0
            )
            #TODO was +b... here in set_poisson 
            A[sb,end-nb+1:end] = b * (-1) * (HxT[iLS] * iMx_b * Hx_b .+ HyT[iLS] * iMy_b * Hy_b) #why was it different in set_poisson ?
            # Boundary conditions for outer boundaries
            A[end-nb+1:end,sb] = b_b * (-1) * (HxT_b * iMx_b' * Hx[iLS] .+ HyT_b * iMy_b' * Hy[iLS])
        end

        veci(rhs,grid,iLS+1) .= +χ[iLS] * vec(a0[iLS]) #was - in set_poisson
    end

    vecb(rhs,grid) .= +χ_b * vec(a0_b) #was - in set_poisson
    
    return rhs
end

```








##### Simultaneous boundary conditions

##### Corners

!!! todo "TODO"
    corners



## Poisson equation with variable coefficient 

The gradient is not collocated with the scalar, so they don't have the same shapes for each direction. Five diagonal matrices with coeffD inside are required: 

```julia
mat_coeffDx = Diagonal(vec(coeffDx_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Bx
mat_coeffDy = Diagonal(vec(coeffDy_bulk)) # coeffDx_bulk is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies By
mat_coeffDx_i = Diagonal(vec(coeffDx_interface)) # coeffDx_interface is a 2d matrix with shape (grid_u.ny, grid_u.nx), multiplies Hx
mat_coeffDy_i = Diagonal(vec(coeffDy_interface)) # coeffDx_interface is a 2d matrix with shape (grid_v.ny, grid_v.nx), multiplies Hy
mat_coeffDx_b = Diagonal(vec(coeffD_borders)) # coeffDx_interface is a 1d vector with shape (2grid.ny + 2grid.nx), multiplies Hx_b and Hy_b
```


\begin{equation}
    \grad\left( p^\omega,p^\gamma \right)=W^\dagger \left( G p^\omega + H p^\gamma \right)
\end{equation}

With variable coefficient:

\begin{equation}
    \kappa \grad\left( p^\omega,p^\gamma \right)=\kappa^\omega W^\dagger \left( G p^\omega + H p^\gamma \right)
\end{equation}



\begin{equation}
    \mathrm{div}\left( \mathbf{q}^\omega,\mathbf{q}^\gamma \right)=-\left( G^T + H^T \right) \mathbf{q}^\omega + H^T \mathbf{q}^\gamma 
\end{equation}




```math
\begin{cases}
-\mathrm{div}\left( \mathbf{q}^\omega, \mathbf{q}^\gamma \right) =& Vf^\omega \\
I_a I_\Gamma p^\gamma - I_b \mathrm{div}\left( 0,\kappa^\gamma \mathbf{q}^\gamma \right) =& I_\Gamma g^\gamma \\
\mathbf{q}^\omega =& \kappa \grad\left( p^\omega,p^\gamma \right)\\
\mathbf{q}^\gamma =& \mathbf{q}^\omega
\end{cases}
```

```math
\begin{cases}
-\mathrm{div}\left( \mathbf{q}^\omega,\mathbf{q}^\gamma \right) =& Vf^\omega \\
I_a I_\Gamma p^\gamma - I_b \mathrm{div}\left( 0,\kappa^\gamma \mathbf{q}^\gamma \right) =& I_\Gamma g^\gamma \\
\mathbf{q}^\omega =& W^\dagger \left( \kappa^\omega G p^\omega + \kappa^\gamma H p^\gamma \right)\\
\mathbf{q}^\gamma =& \mathbf{q}^\omega
\end{cases}
```



```math
\begin{bmatrix}
G^T W^\dagger \kappa^\omega G & G^T W^\dagger \kappa^\gamma H \\
I_b H^T W^\dagger \kappa^\omega G & I_b H^T W^\dagger \kappa^\gamma H + I_a I_\Gamma  
\end{bmatrix}
\begin{bmatrix}
p^\omega \\
p^\gamma  
\end{bmatrix} \simeq
\begin{bmatrix}
V f^\omega \\
I_\Gamma g^\gamma   
\end{bmatrix}
```math


\begin{equation}
    \mathrm{div}(UA)=U(divA)+A⋅\grad(U)
\end{equation}

\begin{equation}
    \mathrm{div}(\kappa \mathbf{q} )=\kappa (\mathrm{div}(\mathbf{q}))+ \mathbf{q}⋅\grad(\kappa)
\end{equation}


% div(kappa grad(A))=U(divgrad(A))+grad(A)⋅gradU

    


```julia
ni = grid.nx * grid.ny
nb = 2 * grid.nx + 2 * grid.ny
nt = (num.nLS + 1) * ni + nb
Aphi_eleL = spzeros(nt, nt)

phL.phi_ele .= reshape(veci(phL.phi_eleD,grid,1), grid)

A[end-nb+1:end,1:ni] = -b_b * (HxT_b * iMx_b' * coeff[end-nb+1:end,1:ni] * Bx .+ HyT_b * iMy_b' * coeff[end-nb+1:end,1:ni] * By)

@inline zeros(g::Grid) = zeros(g.ny, g.nx)

@inline fnzeros(g::Grid, n::Numerical) = zeros((n.nLS + 1) * g.ny * g.nx + 2 * g.nx + 2 * g.ny)

set_border_matrices!
```


```julia
# Implicit part of heat equation
        A[1:ni,1:ni] = pad_crank_nicolson(M .- 0.5 .* τ .* diffusion_coeff_scal .* LT, grid, τ)
        A[1:ni,ni+1:2*ni] = - 0.5 .* τ .* diffusion_coeff_scal .* LD
        A[1:ni,end-nb+1:end] = - 0.5 .* τ .* diffusion_coeff_scal .* LD_b

        # Interior BC
        A[ni+1:2*ni,1:ni] = b * (HxT[1] * iMx * Bx .+ HyT[1] * iMy * By)
        A[ni+1:2*ni,ni+1:2*ni] = pad(b * (HxT[1] * iMx * Hx[1] .+ HyT[1] * iMy * Hy[1]) .- χ[1] * a1)
        A[ni+1:2*ni,end-nb+1:end] = b * (HxT[1] * op.iMx_b * op.Hx_b .+ HyT[1] * op.iMy_b * op.Hy_b)

        # Border BCs
        A[end-nb+1:end,1:ni] = b_b * (op.HxT_b * op.iMx_b' * Bx .+ op.HyT_b * op.iMy_b' * By)
        A[end-nb+1:end,ni+1:2*ni] = b_b * (op.HxT_b * op.iMx_b' * Hx[1] .+ op.HyT_b * op.iMy_b' * op.Hy[1])
        A[end-nb+1:end,end-nb+1:end] = pad(b_b * (op.HxT_b * op.iMx_bd * op.Hx_b .+ op.HyT_b * op.iMy_bd * op.Hy_b) .- op.χ_b * a1_b, 4.0)

        # Explicit part of heat equation
        B[1:ni,1:ni] = M .+ 0.5 .* τ .* diffusion_coeff_scal .* LT .- τ .* CT
        B[1:ni,ni+1:2*ni] = 0.5 .* τ .* diffusion_coeff_scal .* LD
        B[1:ni,end-nb+1:end] = 0.5 .* τ .* diffusion_coeff_scal .* LD_b
```

## Poisson equation (old implementation)

From [`(Rodriguez et al. 2024)`](https://link.springer.com/article/10.1007/s00707-024-04133-4)

We consider the Poisson equation ``-\nabla \cdot \nabla p = f~\text{in}~\Omega`` and the Robin boundary condition ``\frac{\partial p}{\partial n} = q = g``  where a and b are two scalar coefficients


is written:

`` \color{flblue}{\oint_\Gamma q \cdot n dS} = g \chi``
```@raw html
<a name="tagRobin"></a> 
```
Robin boundary condition `` ap + b \frac{\partial p}{\partial n} =  g``


The hypothesis ``q^\omega=q^\gamma=q`` can be made for system resolution, which simplifies the above expression:

```math
\begin{aligned}
\int_\Omega \nabla \cdot q dV = B_i^x (q_{i+ \frac 12}^x-q_{i - \frac 12}^x) + B_i^y (q_{j+ \frac 12}^y-q_{i - \frac 12}^y)
\end{aligned}
```


"The scalar value along the immersed boundary ( ``p^\gamma`` ) is not explicitly known, and has to be inferred from the values of the normal component of the gradient along the immersed  boundary and ``g``. The discretization of this equation requires the definition of the boundary vectors ``g^\gamma`` , ``a^\gamma`` and ``b^\gamma`` for boundary cells, whose components are set to,

``g^\gamma = g(x_{i,j}^\gamma, y_{i,j}^\gamma)``
``a^\gamma = a(x_{i,j}^\gamma, y_{i,j}^\gamma)``
``b^\gamma = b(x_{i,j}^\gamma, y_{i,j}^\gamma)``

where ``a^\gamma`` and ``b^\gamma`` define the diagonal matrices ``I_a`` and ``I_b``,  ``I_a = \mathrm{diag} (a^\gamma)`` and ``I_b = \mathrm{diag} (b^\gamma)``.  Likewise, the components of the discrete load ``f^\omega`` are defined as  ``f^\omega_{i,j} = f  (  x_{i,j}^\gamma, y_{i,j}^\gamma  )``  in fluid cells and zero elsewhere.
"


```@raw html
See this equation for the
<a href="documentation.html#divergencetag">divergence</a> and:
```

"
In order to discretize the Poisson and Robin 
```@raw html
See this equation for the
<a href="test.html#tagPoisson">Poisson</a> and <a href="test.html#tagRobin">Robin</a> equations
```
, the term with the derivative in the normal direction ``(\partial_n p)`` is set as the boundary contribution to the divergence. The resulting system is given by:


```math
\begin{cases}
− \mathrm{div} (q^\omega, q^\gamma ) &= V f^\omega\\
I_a I p^\gamma − I_b \mathrm{div} (0, q^\omega) &= I g^\omega\\
q^\omega &= \mathrm{grad} ( p^\omega, p^\gamma ) \\
q^\gamma &= q^\omega
\end{cases}
```

where ``I`` measures the boundary within a given cell.
"

"
The intermediate variables ``q^\omega`` and ``q^\gamma`` can be eliminated from the above linear system, which can then be expressed in terms of the previously defined matrices as
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


```@raw html
by eliminating all intermediate variables and using Eqs.
<a href="documentation.html#grad_concise">(gradient)</a> and <a href="documentation.html#div_concise">(divergence)</a>.
```
"The matrix on the left-hand side is symmetric if the coefficients b are equal to one and it is straightforward to obtain a symmetric matrix if they are not. Appendix C presents the resulting expression for a given cell of the one-dimensional Poisson's equation, for the sake of clarity."



From [`(Rodriguez et al. 2024)`](https://link.springer.com/article/10.1007/s00707-024-04133-4)
"
$I^\Gamma$ is needed in order to have a dimensionally consistent equation for the 
```@raw html
See this equation for the
<a href="test.html#tagRobin">Robin</a> equations.
```

The intermediate variables ``q^\omega`` and ``q^\gamma`` can be eliminated from the above linear system, which can then be expressed in terms of the previously defined matrices as 


```math
\begin{bmatrix}
G W^\dagger G & G W^\dagger H\\
I_b H W ^\dagger G & I_b H W ^\dagger  H + I_a I_\Gamma 
\end{bmatrix}
\begin{bmatrix}
p^\omega\\
p^\gamma
\end{bmatrix} = 
\begin{bmatrix}
V f^\omega\\
I g^\gamma
\end{bmatrix}
```

"



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

```@docs
set_poisson
```



cf. [`(Rodriguez et al. 2024)`](https://link.springer.com/article/10.1007/s00707-024-04133-4) for a 1D expression in the i-th cell, x component 


```math
\begin{aligned}
-\mathcal{B}_{x,i} [&\mathcal{W}_{x,i+1} (\mathcal{B}_{x,i+1} p^\omega_{i+1} - \mathcal{B}_{x,i} p^\omega _i ) - \mathcal{W}_{x,i} (\mathcal{B}_{x,i} p^\omega_i - \mathcal{B}_{x,i-1} p^\omega_{i-1} )] \\
-\mathcal{B}_{x,i} \{&\mathcal{W}_{x,i+1} [(\mathcal{B}_{x,i+1}  - \mathcal{A}_{x,i+1} )p^\gamma_{i+1} - (\mathcal{A}_{x,i+1} - \mathcal{B}_{x,i} ) p^\gamma_i ] \\
&-\mathcal{W}_{x,i} [(\mathcal{B}_{x,i} - A_{x,i} ) p^\gamma_i + (A_{x,i} - \mathcal{B}_{x,i-1}) p^\gamma_{i-1} ] \} \\
= V_i f^\omega_i&
\end{aligned}
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


### Boundary conditions

#### Neumann boundary condition
The interface length is used to impose boundary conditions, for instance, a Neumann boundary condition ``\frac{\partial p}{\partial n} = q = g`` is written:

`` \color{flblue}{\oint_\Gamma q \cdot n dS} = g \chi``

#### Robin boundary condition 

`` ap + b \frac{\partial p}{\partial n} =  g``

#### Assumption
The hypothesis ``q^\omega=q^\gamma=q`` can be made for system resolution, which simplifies the above expression:

```math
\begin{aligned}
\int_\Omega \nabla \cdot q dV = B_i^x (q_{i+ \frac 12}^x-q_{i - \frac 12}^x) + B_i^y (q_{j+ \frac 12}^y-q_{i - \frac 12}^y)
\end{aligned}
```


### Other formalisms 

The aforementioned operators may be rewritten in a more concise form.

```@raw html
<a name="grad_concise"></a> 
```
```math
\mathrm{grad}(q^ω, q^{γ}  ) = W ^ \dagger \left ( G p ^ \omega + H p ^ \gamma \right )
```


```math
\mathrm{div}(q^ω, q^{γ}  ) = − G^⊤ + H^{Γ1,⊤} + H^{Γ2,⊤} q^ω + H^{Γ1,⊤}q^{γ1} + H^{Γ2,⊤}q^{γ2}
```
where  
```math
−G^⊤ − H^{Γ1,⊤} − H^{Γ2,⊤} = D_x+A_x^\Omega D_y+A_y^\Omega
```
and  
```math
H^{Γi,⊤} = B_x^{Γi} D_x+ − D_x+A_x^{Γi} B_y^{Γi} D_y+ − D_y+A_y^{Γi}
```

"
It should be noted that if the vector field represents a gradient, which is a higher order quantity, no distinction is made between its values in the fluid and at the boundary, i.e. $q^\gamma$ is simply set to $q^ \omega$, simplifying the divergence operator to

```@raw html
<a name="div_concise"></a> 
```
```math
\mathrm{div}(q^ω, q^{γ}  ) = -G ^T q^ \omega
```
" [`Rodriguez et al. 2024`](https://link.springer.com/article/10.1007/s00707-024-04133-4)

cf fig 2.3 in [`Rodriguez 2024`](https://theses.fr/s384455)

Divergence of a face-centered vector field

```@raw html
<figure>
    <a name="geometric_moments"></a> 
    <img src="./assets/moments_doc.svg" alt="Geometric moments" title="Geometric moments">
    <figcaption>"Geometric moments. $\mathcal{A} _ x$ in red dashed lines, $\mathcal{A} _ y$ in blue dashed lines, $\mathcal{B} _ x$ in green dashed lines. $\mathcal{B} _ y$ are not shown."</figcaption>
</figure>
```

cf *A levelset based cut-cell method for two-phase flows. Part 2: Free-surface flows and dynamic contact angle treatment*:

"Just like the gradient operator, the (volume-integrated) divergence operator is multilinear and takes one more argument than the number of interfaces:  

```math
\mathrm{div}(q^ω, q^{γ1} , q^{γ2} ) = − G^⊤ + H^{Γ1,⊤} + H^{Γ2,⊤} q^ω + H^{Γ1,⊤}q^{γ1} + H^{Γ2,⊤}q^{γ2}
```
where  
```math
−G^⊤ − H^{Γ1,⊤} − H^{Γ2,⊤} = D_x+A_x^\Omega D_y+A_y^\Omega
```
and  
```math
H^{Γi,⊤} = B_x^{Γi} D_x+ − D_x+A_x^{Γi} B_y^{Γi} D_y+ − D_y+A_y^{Γi}
```

The reader may find details in [`Rodriguez 2024`](https://theses.fr/s384455) and [`Fullana (2017)`](https://theses.hal.science/tel-04053531/)
for Morinishi's notation for proofs.




### Capacities

See section 2.2.1 Geometric moments in "Numerical methods for modeling and optimization of interfacial flows"
```julia
opC_p.Hx_b: left: -dy (cell height along y), bottom: 0, right: +dy, top: 0
opC_p.Hy_b: left: 0 , bottom: -dx, right: 0, top: +dx
```


## Geometry

Example:

```julia
x = LinRange(0, L0, n+1)
y = LinRange(0, L0, n+1)
```
Or
```julia
x = collect(LinRange(0, L0, n+1))
y = collect(LinRange(0, L0, n+1))
```

That sets the location of the cell faces, not the cell center. The borders of the domain will be at 0 and at L0.


### Along x-axis, without interface, constant mesh spacing

The first cell centroid (when there is no two-fluid interface and the mesh spacing is constant) is at:

```julia
x[1]+dx[1]/2
```

i.e.

```julia
xmin+dx[1]/2
```
The cell centroids are: ``x[1]+dx[1]/2``, ``x[1]+dx[1]*3/2``,..., ``x[end]-dx[1]3/2``,``x[end]-dx[1]/2``

The meshes are initialized by [`init_meshes`](@ref) which in turn calls [`Mesh`](@ref)

```@docs
init_meshes
```

```@docs
Mesh
```

The borders of the domain are defined by x and y, now stored ax ``x\_node`` and ``y\_node`` in the Mesh structure, you can also compute the coordinates of the borders of the domain with: 
```julia
x_bc_left = gp.x[:,1] .- gp.dx[:,1] ./ 2.0

y_bc_bottom = gp.y[1,:] .- gp.dy[1,:] ./ 2.0

y_bc_top = gp.y[end,:] .+ gp.dy[end,:] ./ 2.0

x_bc_right = gp.x[:,end] .+ gp.dx[:,end] ./ 2.0
```

## Storage of bulk and interfacial variables
Variables are stored in the following order: bulk, interfacial (based on levelsets) and interfacial (borders, in the order left bottom right top, cf ``grid.ind.b_left[1], grid.ind.b_bottom[1], grid.ind.b_right[1] and grid.ind.b_top[1]``).


 The system size is ``nt \times nt``:

```julia
ni = gp.nx * gp.ny #size of bulk/interfacial (LS) data
nb = 2 * gp.nx + 2 * gp.ny #size of border 
nt = (num.nLS + 1) * ni + nb 
A = spzeros(nt, nt)
```



Example at ``(j=1,i=1)``, index of corner values
```julia
i_corner_vecb_left = (num.nLS + 1) * ni + 1
i_corner_vecb_bottom = (num.nLS + 1) * ni + gp.ny + 1
```

You can use the following functions to access bulk and interfacial fields:

```@docs
veci
```

```julia
vecb(a, g::G) where {G<:Grid} = @view a[end-2*g.ny-2*g.nx+1:end]
```

```@docs
vecb
```

```julia
vecb_L(a,g::G) where {G<:Grid} = @view vecb(a, g)[1:g.ny]
```

```@docs
vecb_L
```

```julia
vecb_B(a,g::G) where {G<:Grid} = @view vecb(a, g)[g.ny+1:g.ny+g.nx]
```

```@docs
vecb_B
```

```julia
vecb_R(a,g::G) where {G<:Grid} = @view vecb(a, g)[g.ny+g.nx+1:2*g.ny+g.nx]
```

```@docs
vecb_R
```

```julia
vecb_T(a,g::G) where {G<:Grid} = @view vecb(a, g)[2*g.ny+g.nx+1:2*g.ny+2*g.nx]
```

```@docs
vecb_T
```



### Indices
In Flower.jl, the dimensions are stored in this order (y,x) i.e. to access the ``data`` of cell ``(i,j)``, one may use ``data[j,i]``.

```@docs
Indices
```

### Mesh and how to access data at (i,j)
For a 1D vector, grid p
```julia
II = CartesianIndex(jplot, iplot) #(id_y, id_x)
pII = lexicographic(II, grid.ny)
```

!!! todo "TODO CHECK"
    For a 1D vector, grid v
    ```julia
    II = CartesianIndex(jplot, iplot) #(id_y, id_x)
    pII = lexicographic(II, grid.ny +1)
    ```

diag[pII] is faster than [pII,pII]

### Accessing iMx (inverse of weight)

The function lexicographic can be used to access iMx (inverse of weight) in a cell.

```julia
II = CartesianIndex(j,i)
tmp = lexicographic(II,gp.ny)
```

II is a two-dimensional index, and in ``iMx`` 1D index is required.
It is recommended to use ``iMx.diag[id]``, it's faster than accessing ``iMx[id,id]``.

!!! todo "TODO
    Do I need to add 1 to n for the border terms ? ``pII = lexicographic(δy⁺(II), grid.ny)``          

In [`set_cutcell_matrices!`](@ref), ``grid.ny`` already has the good size no matter what grid you're using

But for all the grids, when you set ``iMy``, you have to do

```julia
pII = lexicographic(δy⁺(II), grid.ny+1) #for the top border
pII = lexicographic(II, grid.ny+1) #for the rest
```

Because ``iMy`` has one extra layer on y (true for all the grids).


!!! todo "TODO"
    v mesh left wall : in the x-direction the distances are h and h/2





## Initialisation

When a Neumann Boundary Condition (BC) is imposed at the interface, an extrapolation can be used to get an approximation of the value at the interface. 

```@docs
get_height!
init_fields_multiple_levelsets!
```

get height for u grid: solid then liquid

H,LS.mid_point, dx, dy and geo.centroid are of size (ny, nx+1) 


get height for v grid: solid then liquid

H,LS.mid_point, dx, dy and geo.centroid are of size (ny+1, nx) 


get height for p grid: solid then liquid

H, LS.mid_point, dx, dy and geo.centroid are of size (ny, nx) 

See [Definitions](@ref)

!!! todo "TODO"
    finish doc on centroids
    
    ```julia
    x_centroid = gp.x .+ getproperty.(gp.LS[1].geoL.centroid, :x) .* gp.dx
    y_centroid = gp.y .+ getproperty.(gp.LS[1].geoL.centroid, :y) .* gp.dy

    x_bc = gp.x .+ getproperty.(gp.LS[1].mid_point, :x) .* gp.dx
    y_bc = gp.y .+ getproperty.(gp.LS[1].mid_point, :y) .* gp.dy

    #Initialize bulk value
    vec1(phL.TD,gp) .= vec(ftest.(x_centroid,y_centroid))
    #Initialize interfacial value
    vec2(phL.TD,gp) .= vec(ftest.(x_bc,y_bc))

    #Initialize border value

    vecb_L(phL.TD, gp) .= ftest.(gp.x[:,1] .- gp.x[:,1] ./ 2.0, gp.y[:,1])
    vecb_B(phL.TD, gp) .= ftest.(gp.x[1,:], gp.y[1,:] .- gp.y[1,:] ./ 2.0)
    vecb_R(phL.TD, gp) .= ftest.(gp.x[:,end] .+ gp.x[:,1] ./ 2.0, gp.y[:,1])
    vecb_T(phL.TD, gp) .= ftest.(gp.x[1,:], gp.y[end,:] .+ gp.y[1,:] ./ 2.0)
    ```



## Electrical potential
### Poisson equation for electrical potential

Variable coefficients are interpolated. The Laplacian matrices `L, bc_L_b and bc_L` are defined in the functions [`laplacian`](@ref) and [`laplacian_bc`](@ref). ``coeffD`` is a coefficient involving the electrical conductivity. Then, you have to put the coefficient ``coeffD`` before the matrices ``Bx, By, Hx, Hy, Hx_b`` and ``Hy_b``. The thing is that the gradient is not collocated with the scalar, so they don't have the same shapes for each direction. Five diagonal matrices with coeffD inside are required.


In [`set_borders!`](@ref), we have:
* Dirichlet: a1 = -1
* Neumann: b =1
* Robin: a1 = -1 b = 1
* a0 : BC value 


```@docs
interpolate_scalar!(grid, grid_u, grid_v, u, uu, uv)
```

```@docs
aux_interpolate_scalar!(II_0, II, u, x, y, dx, dy, u_faces)
static_stencil(a, II::CartesianIndex)
faces_scalar(itp, II_0, II, x, y, dx[II], dy[II])
```

```math
\phi(\bar{x}, \bar{y}) = \begin{bmatrix} \bar{x}^2 & \bar{x} & 1 \end{bmatrix}
\begin{bmatrix}
m_{0,0} & m_{0,1} & m_{0,2} \\
m_{1,0} & m_{1,1} & m_{1,2} \\
m_{2,0} & m_{2,1} & m_{2,2}
\end{bmatrix}
\begin{bmatrix} \bar{y}^2 \\ \bar{y} \\ 1 \end{bmatrix} = X^T MY
```
This matrix is computed in [`B_BT`](@ref).
```@docs
B_BT(II_0, x, y)
```
where $\bar{x} = \frac{x - x_{i,j}}{x_{i+1,j} - x_{i-1,j}}$ and $\bar{y} = \frac{y - y_{i,j}}{y_{i,j+1} - y_{i,j-1}}$ are normalized Cartesian coordinates, and $M$ is a matrix with the unknown interpolation coefficients. The values of the level-set at the cells of the $3 \times 3$ stencil can be used to set the following system of equations

```math
\Phi = AMB.
```

where

```math
\Phi = \begin{bmatrix}
\phi_{i-1,j-1} & \phi_{i-1,j} & \phi_{i-1,j+1} \\
\phi_{i,j-1} & \phi_{i,j} & \phi_{i,j+1} \\
\phi_{i+1,j-1} & \phi_{i+1,j} & \phi_{i+1,j+1}
\end{bmatrix},
```

```math
A = \begin{bmatrix}
(\bar{x}_{i-1,j})^2 & \bar{x}_{i-1,j} & 1 \\
0 & 0 & 1 \\
(\bar{x}_{i+1,j})^2 & \bar{x}_{i+1,j} & 1
\end{bmatrix},
```

```math
B = \begin{bmatrix}
(\bar{y}_{i,j-1})^2 & 0 & (\bar{y}_{i,j+1})^2 \\
\bar{y}_{i,j-1} & 0 & \bar{y}_{i,j+1} \\
1 & 1 & 1
\end{bmatrix}.
```

The interpolation matrix $M$ can be computed after inverting matrices $A$ and $B$

```math
M = A^{-1}\Phi B^{-1}.
```

Once the interpolation matrix is obtained, the following equation can be used to interpolate the value of the level-set at a normalized position $(\bar{x}, \bar{y})$:
```math
\phi(\bar{x}, \bar{y}) = \begin{bmatrix} \bar{x}^2 & \bar{x} & 1 \end{bmatrix}
\begin{bmatrix}
m_{0,0} & m_{0,1} & m_{0,2} \\
m_{1,0} & m_{1,1} & m_{1,2} \\
m_{2,0} & m_{2,1} & m_{2,2}
\end{bmatrix}
\begin{bmatrix} \bar{y}^2 \\ \bar{y} \\ 1 \end{bmatrix} = X^T MY
```




"
The level-set has to be interpolated to compute the required geometric moments. 
To that end, a biquadratic interpolation defined on the ```3 \times 3``` stencil centered at a cell ```(i, j)``` will be used. 
The interpolant is defined as
"

```@docs
biquadratic(m, x, y)
```


#### Example 
To determine the values of the corners, we perform a bi-quadratic interpolation on the 3 x 3 stencil centered on the cell of interest. 
For this exercise, we assume a constant spacing of the Cartesian grid in all dimensions. We want to calculate the

The interpolation matrix \( M \) is given by:

\[
\phi(x, y) = \sum_{i=0}^{2} \sum_{j=0}^{2} m_{ij} x^i y^j
\]

This can be expressed as:

\[
\phi(x, y) = \begin{bmatrix} x^2 & x & 1 \end{bmatrix}
\begin{bmatrix}
m_{2,2} & m_{2,1} & m_{2,0} \\
m_{1,2} & m_{1,1} & m_{1,0} \\
m_{0,2} & m_{0,1} & m_{0,0}
\end{bmatrix}
\begin{bmatrix}
y^2 \\
y \\
1
\end{bmatrix}
= X M Y^T
\]

With \(\phi(x, y)\) being the set of known data points of the level set function, yielding:

\[
\Phi = G M G^T
\]

Where:

\[
\Phi = \begin{bmatrix}
\phi_{-1,-1} & \phi_{-1,0} & \phi_{-1,1} \\
\phi_{0,-1} & \phi_{0,0} & \phi_{0,1} \\
\phi_{1,-1} & \phi_{1,0} & \phi_{1,1}
\end{bmatrix}
= \begin{bmatrix}
(-1)^2 & -1 & 1 \\
0^2 & 0 & 1 \\
1^2 & 1 & 1
\end{bmatrix}
G
\begin{bmatrix}
m_{2,2} & m_{2,1} & m_{2,0} \\
m_{1,2} & m_{1,1} & m_{1,0} \\
m_{0,2} & m_{0,1} & m_{0,0}
\end{bmatrix}
M
\begin{bmatrix}
(-1)^2 & 0 & (1)^2 \\
-1 & 0 & 1 \\
1 & 1 & 1
\end{bmatrix}
G^T
\]



!!! todo "Interpolation for Poisson with variable coefficient"
    does not take border and interfacial values into account

 

```@raw html
For the wall term on the left part of the domain, corresponding face control volume <a href="#facecontrolvolume"> (see figure), we should use: $\kappa = \frac{\kappa^\gamma+\kappa^\omega_{1,j}}{2}$.<span><img src="./assets/Wx_mesh_square_bubble_wall_2.svg" alt="image" height="800" /></span></a></span>
```

```@docs
interpolate_scalar_to_staggered_u_v_grids_at_border!
```

```@docs
solve_poisson_variable_coeff!
```
Laplacian: bulk 
```julia
L = BxT * iMx * Bx
```
L is of size nx*ny
tmp_x ((nx+1)*ny , nx*ny)
BxT (nx*ny, (nx+1)*ny )
Bx ((nx+1)*ny, nx*ny)
iMx ((nx+1)*ny, (nx+1)*ny )
iMx_b (nx+1)*ny,2*nx +2*ny
Hx_b  2*nx +2*ny,2*nx +2*ny


The main loop for the resolution of the Poisson equation is in:
```@docs 
solve_poisson_loop!
```

### Electrical current
```@docs
compute_grad_phi_ele!
```

## Scalar transport
```@docs
scalar_transport!
```

## Phase change
```@docs
integrate_mass_flux_over_interface
update_free_surface_velocity_electrolysis!
```

## Electrical conductivity

```@docs
update_electrical_conductivity!(num,grid,elec_cond,elec_condD)
```

## Butler-Volmer equation

```@docs
update_BC_electrical_potential!(num,grid,BC_phi_ele,elec_cond,elec_condD,i_butler)
```

```@docs
update_electrical_current_from_Butler_Volmer!(num,grid,heat,phi_eleD,i_butler;T=nothing)
```


```@docs
butler_volmer_concentration(alpha_a,alpha_c,c_H2,c0_H2,c_H2O,c0_H2O,c_KOH,c0_KOH,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0)
butler_volmer_no_concentration(alpha_a,alpha_c,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0)
butler_volmer_no_concentration!(alpha_a,alpha_c,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0,i_current)
butler_volmer_no_concentration_concentration_Neumann(alpha_a,alpha_c,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0,diffusion_coeff,inv_stoechiometric_coeff)
```

## Heat equation

```@docs
set_heat!
```

## Navier-Stokes


!!! info "Operators"
    "I do this:
    ```julia
    opC_u.AxT * opC_u.Rx
    ```
    so that the same staggered volume is used for all the terms. Otherwise, there were some velocity cells that were not seeing any pressure gradient"



```@docs
FE_set_momentum
```
The resolution is coupled if a Navier BC is activated.
```@docs
FE_set_momentum_coupled
```


```@docs
pressure_projection!
```


In "An important issue is that the boundary condition for u§ must be consistent with Eq. (10), although at the time the boundary conditions are applied the function phi is not yet known and hence must be approximated"

```@docs
set_Forward_Euler!
```

```@docs
set_Crank_Nicolson!
```


## Velocity-pressure coupling




Either pressure gradient in prediction or remove component in velocity boundary conditions


With the parameter ``num.prediction``, the following methods can be called:

* 0: original pressure-velocity coupling in Flower, the pressure gradient is not in the prediction and the bouondary conditions are not modified accordingly
* 1: pressure gradient in prediction
*

"
In the present work, the method referred to as projection method II (PmII) by [Brown et al. 2001](https://www.sciencedirect.com/science/article/pii/S0021999101967154), which ensures a second order discretization of the equations, is employed. The first step consists in updating an intermediate velocity field $\mathbf{u}^*$ using the momentum equation (referred to as prediction step), where the diffusive transport term is solved using a Crank--Nicolson scheme and the convective transport term using a second-order Adams--Bashforth integrator,
"

```math
\frac{\mathbf{u}^* - \mathbf{u}^n}{\tau} + [\mathbf{u} \cdot \nabla \mathbf{u}]^{n+1/2} = -\nabla p^{n-1/2} + \frac{\mu}{2} \Delta (\mathbf{u}^n + \mathbf{u}^*).
```

Here, $\tau$ refers to the time step and the superscript $n$ the time iteration. In general, the velocity field $\mathbf{u}^*$ is not divergence-free, and the next step enforces this condition by first solving the Poisson's equation

```math
\tau \Delta \psi^{n+1} = \nabla \cdot \mathbf{u}^*,
```

for the intermediate pressure field $\psi$, prior to correcting the intermediate velocity field according to

```math
\mathbf{u}^{n+1} = \mathbf{u}^* - \tau \nabla \psi^{n+1}.
```

Finally, the pressure is updated for the velocity prediction at the next time-step:

```math
p^{n+1/2} = p^{n-1/2} + \psi^{n+1} - \frac{\tau \mu}{2} \Delta \psi^{n+1}.
```


"

The discrete cut-cell operators are used to discretize the pressure-projection method. The prediction step reads

```math
\frac{u_{\alpha}^{\omega_{i, *}} - u_{\alpha}^{\omega_{i, n}}}{\tau} + \frac{3}{2} \text{conv}_{\alpha} \left( \boldsymbol{u}^{\omega_{i, n}}, \boldsymbol{u}^{\gamma, n} \right) - \frac{1}{2} \text{conv}_{\alpha} \left( \boldsymbol{u}^{\omega_{i, n-1}}, \boldsymbol{u}^{\gamma, n-1} \right) = -V_{\alpha} \left[ \text{grad} \left( p^{\omega_{i, n-1/2}}, p^{\gamma, n-1/2} \right) \right]_{\alpha} + \frac{1}{2 \text{Re}} \left[ \text{div}_{\alpha} \left( \text{grad}_{\alpha} \left( u_{\alpha}^{\omega_{i, *}}, u_{\alpha}^{\gamma, *} \right) \right) + \text{div}_{\alpha} \left( \text{grad}_{\alpha} \left( u_{\alpha}^{\omega_{i, n}}, u_{\alpha}^{\gamma, n} \right) \right) \right] \quad \forall \ \alpha \in \{x, y\}
```


The boundary conditions applicable to $\boldsymbol{u}^{*}$ are those of the velocity field at the next time step $\boldsymbol{u}^{n+1}$. The discrete Poisson’s equation results in

```math
\text{div} \left( \text{grad} \left( \psi^{\omega_{i, n+1}}, \psi^{\gamma, n+1} \right) \right) = \frac{1}{\tau} \text{div} \left( \boldsymbol{u}^{\omega_{i, *}}, \boldsymbol{u}^{\gamma, *} \right)
```



The velocity field is ultimately corrected as

```math
u_{\alpha}^{\omega_{i, n+1}} = u_{\alpha}^{\omega_{i, *}} - \tau V_{\alpha} \left[ \text{grad} \left( \psi^{\omega_{i, n+1}}, \psi^{\gamma, n+1} \right) \right]_{\alpha} \quad \forall \ \alpha \in \{x, y\}
```



and the pressure is finally updated as

```math
p^{\omega_{i, n+1/2}} = p^{\omega_{i, n-1/2}} + \psi^{n+1} - \frac{\tau}{2 \text{Re}} \text{div} \left( \text{grad} \left( \psi^{\omega_{i, n+1}}, \psi^{\gamma, n+1} \right) \right)
```


where the last term ensures the second order accuracy of the pressure field.

As explained in Sec. 1.1.3, in order to obtain an energy preserving discretization of the incompressible Navier–Stokes equations, the pressure gradient must be equal to the negative transpose of the velocity divergence. This means that in the term $\left[ \text{grad} \left( p^{\omega_{i, n-1/2}}, p^{\gamma, n-1/2} \right) \right]_{\alpha}$, instead of using the definition given by Eq. (2.37), the transpose of the operators used in Eq. (2.39) must be employed. Regarding the boundary conditions for the pressure, to ensure the skew-symmetricity of the gradient and the divergence, opposite boundary conditions must be used. This means that whenever a Dirichlet boundary conditions is used in the velocity, a Neumann boundary condition is employed in the pressure and vice-versa.


"


Either pressure gradient in prediction or remove component in velocity boundary conditions

!!! todo "TODO"
    BC document and test 

"For moving boundaries  the gradient of pressure was removed from the prediction because the simulations weren't very stable."

!!! info ""
    To print information about the matrix in the velocity-pressure coupling:  ``gp.ny +1`` in  lexicographic for v  and ``gp.ny``  in  lexicographic for u 

All the operators in the code are volume integrated, 

so if you want to obtain the actual gradient or divergence you always have to divide by the volume 

(In the finite volume method you have the volume that appears dividing)

For the viscous term, you cannot divide by ``opC_v.iMx_bd``, the viscous term is collocated with the velocity component, and you're dividing by a staggered volume. You have to divide by the collocated volume. You have to divide by ``opC_p.iMy`` .


The parameter 
```julia
num.prediction
```
 is used to choose the prediction method.

cf. [Brown et al. 2001](https://www.sciencedirect.com/science/article/pii/S0021999101967154)

!!! todo "Adams"

```@docs
set_convection!
vector_convection!
```




"
For moving boundaries  the gradient of pressure was removed from the prediction because the simulations weren't very stable."

To print information about the matrix in the velocity-pressure coupling:  
gp.ny +1 in  lexicographic for v  and gp.ny  in  lexicographic for u 


All the operators in the code are volume integrated, so if you want to obtain the actual gradient or divergence you always have to divide by the volume (In the finite volume method you have the volume that appears dividing)

For the viscous term, you cannot divide by ``opC_v.iMx_bd``, the viscous term is collocated with the velocity component, and you're dividing by a staggered volume. You have to divide by the collocated volume. You have to divide by ``opC_p.iMy`` .

p not v



## Convection
!!! todo "multiple LS" 


```math
\begin{equation}
    \text{conv}_\alpha \left( u^\omega, u^\gamma, v^\omega, v^\gamma \right)= C_x \left( \mathbf{u}^\omega \right) v_x^\omega - \frac{1}{2} K_x\left( \right)
\end{equation}
```

```math
\begin{equation}
    \text{conv}_\alpha \left( \mathbf{u}^\omega, \mathbf{u}^\gamma, T^\omega, T^\gamma \right)= \sum_\alpha \left\{ \frac{\delta A_\alpha u_\alpha^\omega \overline{T^\omega}^\alpha }{\delta \xi_\alpha} + \left[ \frac{\delta\left( \overline{B_\alpha}^\alpha-A_\alpha\right) \mathbf{u}_\alpha^\gamma}{\delta \xi_\alpha} -\overline{\frac{\delta B_\alpha \mathbf{u}_\alpha^\gamma}{\delta \xi_\alpha}}^\alpha\right]\frac{T^\omega + T^\gamma}{2}\right\}
\end{equation}

```


In [`set_convection!`](@ref) which uses  [`vector_convection!`](@ref)
```julia
Du_y[1,:] .= vecb_B(uD,grid_u) + dt* grd_x[1,:]
```

Extract from [`(Rodriguez et al. 2024)`](https://link.springer.com/article/10.1007/s00707-024-04133-4):

"
This section presents the proposed discretization of the incompressible Navier-Stokes equations for an isotropic Newtonian fluid

```math
\begin{cases}
\frac{\partial u}{\partial t} + \left ( u \cdot \nabla \right ) u & = -\dfrac{\nabla p}{\rho} + \dfrac{1}{\mathrm{Re}} \nabla \cdot \nabla u \\
    \nabla \cdot u & = 0
\label{eq:ns}
\end{cases}
```

where $u$ and $p$ respectively denote the fluid's velocity and pressure fields, $\rho$ its constant density and $\mathrm{Re}$ denotes the Reynolds number. This equation can be rewritten in semi-discrete form using the linear operators previously defined, yielding the momentum equation in direction $\alpha$

\begin{equation} \label{eq:nsSemi}
V _ \alpha \frac{\mathrm{d} u _ \alpha ^ \omega}{\mathrm{d} t} + \mathrm{conv} _ \alpha \left ( u ^ \omega, u ^ \gamma \right ) = - \dfrac{1}{\rho} \left [ \mathrm{grad} \left ( p ^ \omega, p ^ \gamma \right ) \right ] _ \alpha + \dfrac{1}{\mathrm{Re}} \mathrm{div} _ \alpha \left ( \mathrm{grad} _ \alpha \left ( u _ \alpha ^ \omega, u _ \alpha ^ \gamma \right ) \right ),
\end{equation}

and the divergence-free condition
\begin{equation} \label{eq:cont}
    \mathrm{div} \left ( u ^ \omega, u ^ \gamma \right ) = 0,
\end{equation}
where $u ^ \omega$ and $u ^ \gamma$ are two face-centered vector fields and $p ^ \omega$ and $p ^ \gamma$ are two cell-centered scalar fields. [...] The method referred to as projection method II (PmII) by [Brown et al. 2001](https://www.sciencedirect.com/science/article/pii/S0021999101967154) is employed. It ensures a second order discretization of the equations.

In this projection method, the convective term is discretized using the explicit second-order Adams-Bashforth scheme:


!!! todo "boundary conditions"

 and the viscous term is discretized using the implicit Crank-Nicolson scheme. The first step of the method consists of obtaining an intermediate velocity field $u ^ \star$ by solving:
```math
\begin{multline} \label{eq:pmPred}
V _ \alpha \frac{u _ \alpha ^ {\omega, \star} - u _ \alpha ^ {\omega, n}}{\tau} + \frac{3}{2} \mathrm{conv} _ \alpha \left ( u ^ {\omega, n}, u ^ {\gamma, n} \right ) \\ 
- \frac{1}{2} \mathrm{conv} _ \alpha \left ( u ^ {\omega, n - 1}, u ^ {\gamma, n - 1} \right ) = - \dfrac{1}{\rho} \left [ \mathrm{grad} \left ( p ^ {\omega, n - 1 / 2}, p ^ {\gamma, n - 1 / 2} \right ) \right ] _ \alpha \\
+ \dfrac{1}{2 \mathrm{Re}} \left [ \mathrm{div} _ \alpha \left ( \mathrm{grad} _ \alpha \left ( u _ \alpha ^ {\omega, \star}, u _ \alpha ^ {\gamma, \star} \right ) \right ) + \mathrm{div} _ \alpha \left ( \mathrm{grad} _ \alpha \left ( u _ \alpha ^ {\omega, n}, u _ \alpha ^ {\gamma, n} \right ) \right ) \right ]\ \forall\ \alpha \in \{ x, y \},
\end{multline}
```

where $\tau$ denotes the time step and the superscript $n$ the iteration number. The boundary conditions applicable to $u ^ \star$ (the predicted velocity field) are those of the velocity field at the next time step. 


The Dirichlet boundary conditions applied for example to impose the no-slip condition or the inflow conditions read
\begin{equation}
    u _ \alpha ^ {\gamma, \star} = u _ {\alpha, \mathrm{ref}}, \quad \forall \alpha \in \{ x, y \},
\end{equation}
where $u _ {\alpha, \mathrm{ref}}$ is the value of the velocity to be imposed at the boundary in the $\alpha$-direction. The homogeneous Neumann boundary condition applied to impose outflow conditions is given by
\begin{equation}
    -\mathrm{div} _ \alpha \left ( 0, \, \mathrm{grad} _ \alpha \left ( u _ \alpha ^ {\omega, \star}, u _ \alpha ^ {\gamma, \star} \right ) \right ) = 0, \quad \forall \alpha \in \{ x, y \},
\end{equation}

In the projection step, the velocity field is updated by projecting $u ^ \star$ using the intermediate pressure field $\psi ^ {n + 1}$, which is obtained by solving the following Poisson equation
\begin{equation} \label{eq:pmPoisson}
\mathrm{div} \left ( \mathrm{grad}  \left ( \psi ^ {\omega, n + 1}, \psi ^ {\gamma, n + 1} \right ) \right ) = \dfrac{\rho}{\tau}\mathrm{div} \left ( u ^ {\omega, \star}, u ^ {\gamma, \star} \right ).
\end{equation}
Whenever no-slip or inflow conditions are to be imposed, a homogeneous Neumann boundary condition is applied to the intermediate pressure field following
\begin{equation}
    -\mathrm{div} \left ( 0, \, \mathrm{grad} \left ( \psi ^ {\omega, n + 1}, \psi ^ {\gamma, n + 1} \right ) \right ) = 0.
\end{equation}
With outflow boundary conditions, a Dirichlet boundary condition is imposed following
\begin{equation}
    \psi ^ {\gamma, n + 1} = \psi _ {\mathrm{ref}},
\end{equation}

The velocity field is ultimately corrected as
\begin{equation}
u_\alpha ^ {\omega, n + 1} = u ^ {\omega, \star} _\alpha - \dfrac{\tau}{\rho} \left [ \mathrm{grad} \left ( \psi ^ {\omega, n + 1}, \psi ^ {\gamma, n + 1} \right ) \right ] _ \alpha\ \forall\ \alpha \in \{ x, y \}.
\label{eq:pm_projection}
\end{equation}
and the pressure is finally updated as
\begin{equation}
p ^ {\omega, n + 1 / 2} = p ^ {\omega, n - 1 / 2} + \psi ^ {n + 1} - \frac{\tau}{2 \mathrm{Re}} \mathrm{div} \left ( \mathrm{grad} \left ( \psi ^ {\omega, n + 1}, \psi ^ {\gamma, n + 1} \right ) \right ),
\label{eq:pm_update}
\end{equation}
where the last term ensures the second order accuracy of the pressure field.

"


## Boundary conditions

### Boundary conditions for border
The idea there is to have at the corners the contributions from both borders. 
bcTx is only used for derivatives in x and bcTy for derivatives in y, so it doesn't matter that
bcTx is not filled at the top and bottom and the same for bcTy. 
Otherwise we would have to select one of the contributions.

```@docs
set_borders!
```

## Interface definition

Liquid cells have an isovalue of 0, solid cells fave an isovalue of 15.0 and mixed cells have an other value (cf marching squares).

### Convention, orientation of the surface
The signed distance is positive in the liquid and negative in the bubble. The normal is oriented towards the liquid.


```@docs
update_all_ls_data
```

```@docs
marching_squares!
get_cells_indices
```

## Interface transport (with Levelset)


```@docs
IIOE_normal!
IIOE_normal_indices_2!
indices_extension
inflow_outflow
field_extension!
S2IIOE!
```


## Operators

### Gradient of a scalar field

```@docs
compute_grad_T_x_T_y_array!
compute_grad_T_x_array!
compute_grad_T_y_array!
```

```@docs
compute_grad_p!
normal_gradient
```


### Laplacian

```@docs
set_matrices!
```

```@docs
laplacian
laplacian_bc

```





There are several ways to present the method: G and H, Morinishi's notations and a geometric formulation: [`(Rodriguez et al. 2024)`](https://link.springer.com/article/10.1007/s00707-024-04133-4).
```@docs
diamond
```


### Functions for divergence

```math
-G^T= [ B_xD_x^+ B_yD_y^+ ]
```

with 

```math
\begin{equation}
( D ^ - _ x ) ^ \top = -D ^ + _ x \quad \mathrm{and} \quad ( S ^ - _ x ) ^ \top = S ^ + _ x,
\label{eq:multitrans}    
\end{equation}
```

```math
\Delta ^ + = \left [ \begin{array}{rrrrr}
  -1   &    1   &        &        &         \\
       &   -1   &    1   &        &         \\
       &        & \ddots & \ddots &         \\
       &        &        &   -1   &     1   \\
       &        &        &        &     0
\end{array} \right ]
```

```@docs
divergence_B!
```

Then, for the gradient, ``B_x`` corresponds to ``-BxT'`` i.e. ``-(-G^T)^T=G``




### Updating the operator from Levelset
!!! todo "Order of updates" 
    At every iteration, [`update_all_ls_data`](@ref) is called twice, once inside run.jl and another one 
    (if there's advection of the levelset) inside [`set_heat!`](@ref). 
    The difference between both is a flag as last argument, inside run.jl is implicitly defined as true 
    and inside [`set_heat!`](@ref) is false. If you're calling your version of [`set_heat!`](@ref) several times,
    then you're calling the version with the flag set to false, but for the convective term it has to be set to true.

The flag=true, the capacities are set for the convection, the flag=false they are set for the other operators


In Flower.jl:


```@docs
DiscreteOperators
```


In scalar_transport, this corresponds to: 

```julia
# Interior BC
A[sb,1:ni] = b * (HxT[iLS] * iMx * Bx .+ HyT[iLS] * iMy * By)

# Contribution to Neumann BC from other boundaries
for i in 1:num.nLS
    if i != iLS
        A[sb,i*ni+1:(i+1)*ni] = b * (HxT[iLS] * iMx * Hx[i] .+ HyT[iLS] * iMy * Hy[i]) #-b in Poisson
    end
end

A[sb,sb] = pad(b * (HxT[iLS] * iMx * Hx[iLS] .+ HyT[iLS] * iMy * Hy[iLS]) .- op.χ[iLS] * a1)

#TODO no need for fs_mat

A[sb,end-nb+1:end] = b * (HxT[iLS] * op.iMx_b * op.Hx_b .+ HyT[iLS] * op.iMy_b * op.Hy_b)
```

### Interpolation
```@docs
Acpm
interpolate_grid_liquid!
interpolated_temperature
static_stencil
```

### Interpolating to display solid and liquid

```julia
@views fwd.trans_scal[1,:,:,iscal] .= phL.trans_scal[:,:,iscal].*LS[end].geoL.cap[:,:,5] .+ phS.trans_scal[:,:,iscal].*LS[end].geoS.cap[:,:,5]
```


### Convection

## Capacities

The mass matrix: 
M: face-centered mass matrices, diagonal with coefficients ``V_\alpha`` (the volume of the staggered control volumes)
is stored in the struct Operators and defined in:  
```julia   
M.diag .= vec(geo[end].dcap[:,:,5]) in [`set_matrices!`](@ref)
```
Example: 
```julia   
M_vL = Diagonal(fzeros(grid_v))
```   
Example:
V(1,1)  ⋅           ⋅

⋅       V(...,...)  ⋅

⋅       ⋅           V

New capacities from ``x^2``



```@docs
bc_matrix_borders!
capacities
mat_assign_T!
mass_matrix_borders!
periodic_bcs_borders!
set_boundary_indicator!
```

!!! todo "Rx Ry"

```@docs
periodic_bcs_R!
```


B vs H

```math
\chi
```
```julia
for iLS in 1:num.nLS
    χx = (geo.dcap[:,:,3] .- geo.dcap[:,:,1]) .^ 2
    χy = (geo.dcap[:,:,4] .- geo.dcap[:,:,2]) .^ 2
    χ[iLS].diag .= sqrt.(vec(χx .+ χy))
end
```

## Small cells

Small cells are cut with ``num.\epsilon`` and ``num.\epsilon_{wall}``. This is used to prevent abnormally high values

The mid points are not exactly on the interface when the epsilon parameter is not null:

![](./assets/interface_eps0.png)
*``num.\epsilon = 0 ``*

![](./assets/interface_eps0v05.png)
*``num.\epsilon = 5% ``*


## Solver

For LS: gmres
```julia
grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u) .+ rhs_LS), grid)
```

Elsewhere: Julia's solver: /

* Diagonal scaling does not help (direct resolution)
* MUMPS similar


Changing
* diag + eps
* eps small cell

with diagonal scaling (not useful for direct resolution?): norm2(Ax-b)/norm2(b)=21.5550725612907207e-12

## System

``A[end-nb+testn,end-nb+1:end]``  is the contribution to the gradient in the normal direction from the border, and
  ``A[end-nb+testn,ni+1:2*ni]`` from the interface


If the interface is nearly perpendicular to the border, there's no much contribution from the interface to the gradient in the x-direction.





## Cutcell

```@docs
ilp2cap
```

```@docs
set_cutcell_matrices!
set_other_cutcell_matrices!
set_border_matrices!
```

### Levelset

!!! todo "TODO"
    The levelset is not extended for the evaluation of capacities. At the borders, the interpolation is shifted one cell. You can see this in [marching_squares!](@ref)  (the ``II\_0`` variable).  

#### Ghost cells
```@docs
allocate_ghost_matrices_2
init_ghost_neumann_2
```

## Temporal scheme

!!! todo "Time-step
    change coefficients for extrapolation in case of adaptive time-step


### Crank-Nicolson

### Forward-Euler

### Initialization

#### Interface position for storage

"Moreover, we define the interface centroid as the mid point of the segment crossing the cell which will be used in the computation of the Stefan condition." \citet{fullanaSimulationOptimizationComplex2022}

Store segments in 2D to save up space

```julia 
#Initialize liquid phase

x_centroid = gp.x .+ getproperty.(gp.LS[1].geoL.centroid, :x) .* gp.dx
y_centroid = gp.y .+ getproperty.(gp.LS[1].geoL.centroid, :y) .* gp.dy

x_bc = gp.x .+ getproperty.(gp.LS[1].mid_point, :x) .* gp.dx
y_bc = gp.y .+ getproperty.(gp.LS[1].mid_point, :y) .* gp.dy

#Initialize bulk value
vec1(phL.TD,gp) .= vec(ftest.(x_centroid,y_centroid))
#Initialize interfacial value
vec2(phL.TD,gp) .= vec(ftest.(x_bc,y_bc))
```

## Convection

!!! todo "Advection equation of a vector-valued field  2.5.1 Discrete operators for staggered quantities" [`Rodriguez (2024)`](https://theses.fr/s384455)



## Small cells

Small cells introduce errors, the parameter ``\epsilon`` is used to delete small cells. Merging the cells could help with this problem, though it would involve modifying the capacities and the structure of the cut-cell operators.


Need concentration near wall for conductivity for Poisson equation                      
take the value of the concentration in the bulk field instead 

```julia
BC_phi_ele.left.val .= -i_butler./elec_cond[:,1]
```

for the conductivity used in the Poisson equation, instead of doing :

```julia
BC_phi_ele.left.val .= -i_butler./vecb_L(elec_condD, grid)
```


The conductivity is computed from the scalar. Taking the value in the bulk field is recommended as long as cell merging is not implemented. Due to small cells, we may have slivers/small cells at the left wall, then the divergence term is small,
which produces a higher concentration in front of the contact line.


## Solver

\section{Solvers}

In Flower for LS: gmres
```julia
grid.LS[iLS].u .= reshape(gmres(grid.LS[iLS].A, grid.LS[iLS].B * vec(grid.LS[iLS].u) .+ rhs_LS), grid)
```

Elsewhere: Julia's solver: /
* Diagonal scaling does not help (direct resolution)
* MUMPS similar

Changing
* diag + eps
* eps small cell

with diagonal scaling (not useful for direct resolution?): norm2(Ax-b)/norm2(b)=21.5550725612907207e-12


## Possible improvements 
* Restart (TODO restart with PDI) cf hello_access.jl
```julia
if sim.restart == 1:
    cf hello_access.jl
end
```
* Interface length
* Remove ``replace!(N, NaN=>0.0)``
* imposed constant velocity field : after 100 it : error after scalar transport 3.88e-14
* pressure BC during prediction 
* small cells: linear algebra, time step ...
* Epsilon to prevent NaN, with eps : coefficients not exact 
* sometimes Julia reports out of memory 
* Compilation time: upcoming development: \url{https://jbytecode.github.io/juliac/}
* Store segments in 2D to save up space
* Memory, allocations
* check conservation of variables after n iterations
* stencils
* parallelisation
*  check BC wall pressure (wall no slip : Neumann it seems)
* List of variables which can be removed:
    * emptied
    * cap (already dcap)
    * ...

List all variables
* bulk+ interface \rightarrow *2
* phL.u phL.uD duplicate for bulk 
* u, v, T, LS, moments (heights,...), phi elec, i current, scalars ...
* grids
* liquid/solid
* successive substitution reuse LU decomposition



### Epsilon/handling small or cancelled cells

* ``num.\epsilon`` is uded for small cells
* weights for cancelled cells

```@docs
inv_weight_clip
```


### Sources of errors
* [`inv_weight_clip`](@ref)
* eps added to diagonal (``A[i,i] += 1e-10``)
* epsilon for small cells (to prevent NaN) (``num.\epsilon``), the coefficients are no longer exact in simple cases (like Laplacian: 1 1 -4 1 1)
* Approximation of gradient inside a cell (``q^\omega \approx q^\gamma``)
* approximation of interface





## IO/post-processing

```@docs
convert_interfacial_D_to_segments
```

The output files are generated with [`PDI`](https://github.com/pdidev). The on-line PDI documentation is available at [`PDI documentation`](https://pdi.dev). 

## TODO
Need details on dimensions: construction not explained (sparsity)

## References

* [`Fullana (2017)`](https://theses.hal.science/tel-04053531/)
* [`Rodriguez et al. 2024`](https://link.springer.com/article/10.1007/s00707-024-04133-4)
* [`Rodriguez (2024)`](https://theses.fr/s384455)
* [`Brown et al. (2001)`](https://www.sciencedirect.com/science/article/pii/S0021999101967154)