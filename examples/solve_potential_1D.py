# Script to compare Newton and successive substitutions (also known as Picard iteration)for 1D Poisson equation with constant conductivity
# Butler-Volmer on the left boundary condition and zero potential on the right
# Examples can be found in Langtangen's "Computational Partial Diﬀerential Equations"
# Note : did not multiply by the conductivity in the matrix (constant and electroneutrality)
# Possible improvements : choose U during iterations, linearization


import numpy as np
from scipy.linalg import solve
from scipy.optimize import fsolve
import pandas as pd
from termcolor import colored

# Define a formatting function for scientific notation with three significant digits
def sci_notation(x):
    return f"{x:.3e}"

# For analytical solution
def f(x, L, i0, alpha, Faraday, Ru, temperature0, elec_cond):
    return 2*i0/elec_cond*np.sinh((alpha*Faraday)/(Ru*temperature0)*(-0.6-(-x*L))) + x

# Function to compute F(U)
def compute_residual(U,elec_cond,fact,h,N):
    F = np.zeros(N+1)
    for i in range(1, N):
        F[i] = -(U[i+1] - 2*U[i] + U[i-1]) / h**2

    F[0] = (U[1] - U[0]) / (h) + 1.0/elec_cond * 2*np.sinh(fact * (-0.6 - U[0]))
    F[N] = U[N]  # Right boundary condition
    return F

# Function to compute the Jacobian J(U)
def compute_jacobian(U,elec_cond,fact,h,N):
    J = np.zeros((N+1, N+1))
    for i in range(1, N):
        J[i, i-1] = -1 / h**2
        J[i, i] = 2 / h**2
        J[i, i+1] = -1 / h**2

    J[0, 0] = -1 / (h) - 1.0/elec_cond * fact * 2*np.cosh(fact * (-0.6 - U[0]))
    J[0, 1] = 1 / (h)
    J[N, N] = 1  # Right boundary condition
    return J

# Newton-Raphson iteration
def newton(U_first_guess,elec_cond,fact,solved_slope,h,N,tol,max_iter):
    U = U_first_guess.copy()
    print('U first guess',U)

    columns = ['k', 'phi_wall', 'diff']
    df = pd.DataFrame(columns=columns)
    df['k'].astype("Int64")

    for k in range(max_iter):
        F_U = compute_residual(U,elec_cond,fact,h,N)
        J_U = compute_jacobian(U,elec_cond,fact,h,N)
        delta_U = solve(J_U, -F_U)
        U += delta_U
        slope = (U[-1]-U[0])/(xmax-xmin)
        diff = abs(slope-solved_slope)/solved_slope
        print("iter {} phi at wall {:.3e} relative error {:.3e} ".format(k+1,U[0],diff))

        # Create the list of values
        data_list = [k, U[0], diff]

        # Append the list to the DataFrame
        df.loc[len(df)] = data_list

        if np.linalg.norm(delta_U) < tol:
            break

    return U,df
   

def assemble(h,N):
    A = np.zeros((N+1, N+1))
    for i in range(1, N):
        A[i, i-1] = -1 / h**2
        A[i, i] = 2 / h**2
        A[i, i+1] = -1 / h**2

    A[0, 0] = -1 / (h) 
    A[0, 1] = 1 / (h)
    A[N, N] = 1  # Right boundary condition
    return A


# Successive substitution
def successive_substitutions(U_first_guess,elec_cond,fact,solved_slope,h,N,tol,max_iter,relaxation=1.0):
    
    U = U_first_guess.copy()
    U_previous = U.copy()

    columns = ['k', 'phi_wall', 'diff']
    df = pd.DataFrame(columns=columns)
    df['k'].astype("Int64")

    print('U first guess',U)
    rhs = np.zeros(N+1)

    A=assemble(h,N)

    for k in range(max_iter):
        rhs[0] = -1.0/elec_cond * 2*np.sinh(fact * (-0.6 - U[0]))
        rhs[N] = 0.0
        U = solve(A, rhs)

        slope = (U[-1]-U[0])/L
        diff = abs(slope-solved_slope)/solved_slope

        print("iter {} phi {:3e} relative error {:3e}".format(k+1,U[0],diff))

        # Create the list of values
        data_list = [k, U[0], diff]
        
        # print('U_previous',U_previous[0],(U_previous[0]+solved_slope*L)/(solved_slope*L))
        # print('U',U[0],(U[0]+solved_slope*L)/(solved_slope*L))

        U_new = relaxation * U + (1.0-relaxation) * U_previous 
        U_previous = U.copy()
        U = U_new.copy()

        # print('U_previous',U_previous[0],(U_previous[0]+solved_slope*L)/(solved_slope*L))
        # print('U',U[0],(U[0]+solved_slope*L)/(solved_slope*L))

        # Append the list to the DataFrame
        df.loc[len(df)] = data_list


        # if np.linalg.norm(delta_U) < tol:
        #     break

    return U,df


# Parameters
xmin = 0.0
xmax = 1e-4

#Parameters for Newton/successive substitution
tol = 1e-8
max_iter = 10

#Parameters for fsolve in Python to get an "analytical" solution (by assuming linearity of phi)
xtol=1e-14
maxfev = 1000


print_latex = True
print_latex = False


# Mesh
N = 32 
L= xmax - xmin
h = L / N 

mesh_x_coord = [h*float(i) for i in range(N+1)]
print('mesh_x_coord',mesh_x_coord)


Ru = 8.314
Faraday = 9.64853321233100184e4 #C⋅mol−1
temperature0=353.0
alpha = 0.5
i0 = 1.0
c0_KOH=6700
DKOH=3.2e-9

fact = alpha * Faraday / (Ru * temperature0)

elec_cond = 2*Faraday**2*c0_KOH*DKOH/(Ru*temperature0)
print("elec_cond",elec_cond) # elec_cond = 136.0168282397691


#"analytical" solution if we assume linearity of potential with fsolve 
solved_slope = fsolve(f, 0.0, args=(L, i0, alpha, Faraday, Ru, temperature0, elec_cond),xtol=xtol,maxfev=maxfev)[0]
print('with fsolve, slope:',solved_slope,' residual ',f(solved_slope,L, i0, alpha, Faraday, Ru, temperature0, elec_cond),'phi at electrode ',-solved_slope*L)

print('i',-elec_cond * solved_slope)

# Initial guess
U = np.zeros(N+1)

#Change U[0] if you want to change the initial guess for the potential at the electrode
# For example U[0] = -1e-2


U_first_guess = np.zeros(N+1)
print(colored('\nNewton', 'cyan'))
U_newton,df_newton = newton(U_first_guess,elec_cond,fact,solved_slope,h,N,tol,max_iter)

# Output the solution
# print("Solution U:", U_newton)

slope_newton = (U_newton[-1]-U_newton[0])/L


# print(df_newton)

# diff = [(U_newton[i+1]-U_newton[i])/h for i in range(N)]

# print('x component of gradient',diff)

if print_latex:
    # Convert the DataFrame to a LaTeX table with scientific notation and three significant digits
    latex_table = df_newton.to_latex(index=False,formatters={
    'k': lambda x: f'{x:.0f}',
    'phi_wall': lambda x: f'{x:.3e}',
    'diff': lambda x: f'{x:.3e}'
    })

    # Print the LaTeX table
    print(latex_table)

   
    print(df_newton.to_html(index=False,
        formatters={
        'k': lambda x: f'{x:.0f}',
        'phi_wall': lambda x: f'{x:.3e}',
        'diff': lambda x: f'{x:.3e}'
        }))
    # formatters={"name": str.upper},
    
print(colored('\nSuccessive substitutions', 'cyan'))
U_successive,df_successive = successive_substitutions(U_first_guess,elec_cond,fact,solved_slope,h,N,tol,max_iter)

# Output the solution
# print("Solution U:", U_successive)

slope_successive = (U_successive[-1]-U_successive[0])/(xmax-xmin)

# print(df_successive)

if print_latex:
    # Convert the DataFrame to a LaTeX table with scientific notation and three significant digits
    latex_table = df_successive.to_latex(index=False,formatters={
    'k': lambda x: f'{x:.0f}',
    'phi_wall': lambda x: f'{x:.3e}',
    'diff': lambda x: f'{x:.3e}'
    })
    
    # Print the LaTeX table
    print(latex_table)

    print(df_successive.to_html(index=False,
    formatters={
    'k': lambda x: f'{x:.0f}',
    'phi_wall': lambda x: f'{x:.3e}',
    'diff': lambda x: f'{x:.3e}'
    }))

# diff = [(U_successive[i+1]-U_successive[i])/h for i in range(N)]

# print('x component of gradient',diff)

print(colored('\nSuccessive substitutions with relaxation', 'cyan'))
U_successive_relax,df_successive_relax = successive_substitutions(U_first_guess,elec_cond,fact,solved_slope,h,N,tol,max_iter,relaxation=0.5)

# Output the solution
# print("Solution U:", U_successive_relax)

slope_successive_relax = (U_successive_relax[-1]-U_successive_relax[0])/(xmax-xmin)

# print(df_successive)

if print_latex:
    # df_successive_relax['k'].astype("Int64")

    # Convert the DataFrame to a LaTeX table with scientific notation and three significant digits
    latex_table = df_successive_relax.to_latex(index=False,formatters={
    'k': lambda x: f'{x:.0f}',
    'phi_wall': lambda x: f'{x:.3e}',
    'diff': lambda x: f'{x:.3e}'
    })

    # Print the LaTeX table
    print(latex_table)

    print(df_successive_relax.to_html(index=False,
    formatters={
    'k': lambda x: f'{x:.0f}',
    'phi_wall': lambda x: f'{x:.3e}',
    'diff': lambda x: f'{x:.3e}'
    }))

# diff = [(U_successive_relax[i+1]-U_successive_relax[i])/h for i in range(N)]

# print('x component of gradient',diff)
print('\n')

print('fsolve                                  , slope:{:.3e} residual {:.3e} phi at electrode {:.3e}'.format(
    solved_slope,f(solved_slope,L, i0, alpha, Faraday, Ru, temperature0, elec_cond),-solved_slope*L))

print('Newton                                  , slope:{:.3e} residual {:.3e} phi at electrode {:.3e}'.format(
    slope_newton,f(slope_newton,L, i0, alpha, Faraday, Ru, temperature0, elec_cond),-slope_newton*L))

print('slope_successive                        , slope:{:.3e} residual {:.3e} phi at electrode {:.3e}'.format(
    slope_successive,f(slope_successive,L, i0, alpha, Faraday, Ru, temperature0, elec_cond),-slope_successive*L))

print('successive substitutions, relaxation 0.5, slope:{:.3e} residual {:.3e} phi at electrode {:.3e}'.format(
    slope_successive_relax,f(slope_successive_relax,L, i0, alpha, Faraday, Ru, temperature0, elec_cond),-slope_successive_relax*L))

# print('with fsolve, slope:',solved_slope,' residual ',f(solved_slope,L, i0, alpha, Faraday, Ru, temperature0, elec_cond),'phi at electrode ',-solved_slope*L)

# print('with newton, slope:',slope_newton,' residual ',f(slope_newton,L, i0, alpha, Faraday, Ru, temperature0, elec_cond),'phi at electrode ',-slope_newton*L)

# print('with successive substitutions, slope:',slope_successive,' residual ',f(slope_successive,L, i0, alpha, Faraday, Ru, temperature0, elec_cond),'phi at electrode ',-slope_successive*L)

# print('with successive substitutions, relaxation 0.5, slope:',slope_successive_relax,' residual ',f(slope_successive_relax,L, i0, alpha, Faraday, Ru, temperature0, elec_cond),'phi at electrode ',-slope_successive_relax*L)
