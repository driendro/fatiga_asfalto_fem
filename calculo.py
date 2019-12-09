#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import geometria as g
import fenics as fnc
import numpy as np
import os
import matplotlib.pyplot as plt

# Dimensiones en cm
l, h_prob, b = 30.0, 5.0, 10.0	# largo, alto, ancho del asfalto
e_geo = 0.4 					# Espesor del geosintetico
h_base = 2.0  					# alto base de hormigon
e = 0.5 						# espesor de la ranura
d = 1.0 						# ancho de la zonad de detalle

resmax=0.8						# tamaño de los elementos de los extremos
resmin=0.5						# tamaño de los elementos del centro

archivo = 'calc/modelo'			# PATH/AL/ARCHIVO
archivovtu = 'calc/vtu/modelo'	# PATH/A/LOS/ARCHIVOS/VTU

# No tocar esta lista
parametros=[l, h_prob, b, e_geo, h_base, e, d, resmax, resmin, archivo]

# Mostrar mallado
malla=False
# Mostrar geometria
geometria=False

# Guardamos los parametros del mallado en un archivo llamado "parametros"
# pero antes de guardar los nuevos parametros, copiamos el archivo original
# y lo guardamos como parametros.old
s= 'cp calc/parametros calc/parametros.old'
os.system(s)
a = open('calc/parametros', 'w')
s = parametros
s.append(malla)
s.append(geometria) 
a.write(str(s))
a.close()

# Comparamos ambos archivos para verificar si existen cambios en los parametros
# en caso afirmativo, recalculamos la malla

a = open('calc/parametros', 'r')
file1= a.read()
a.close()
a = open('calc/parametros.old', 'r')
file2=a.read()
a.close()

# comparamos la configuracion actual con la anterior, si son iguales no recalcula
if file1==file2:
	print('parametros anteriores = %s' %(file2))
	print('parametros actuales   = %s' %(file1))
	print('Al ser iguales no se regenerla la malla')
else:
	print('parametros anteriores = %s' %(file2))
	print('parametros actuales   = %s' %(file1))
	print('Se regenerla la malla')
	# llamo a la funcion para regenerar el calculo del mashado completo
	g.probeta(parametros, malla, geometria)

# Obtenemos la malla creada en el archivo geometria.py
mesh = fnc.Mesh("%s.xml" %(archivo)) # mesh
subdominio = fnc.MeshFunction('size_t', mesh, "%s_physical_region.xml" %(archivo)) # Volumenes
bordes = fnc.MeshFunction('size_t', mesh, "%s_facet_region.xml" %(archivo)) # Bordes

#Propidades mecanicas de los materiales
class E(fnc.UserExpression):
	def __init__(self, subdominio, E, **kwargs):
		super().__init__(**kwargs)
		self.subdominio = subdominio
		self.E = E

	def eval_cell(self, values, x, cell):
		# Fijo
		if self.subdominio[cell.index] == 3:
			values[0]=self.E[0]
		# Movil
		elif self.subdominio[cell.index] == 4:
			values[0]=self.E[0]
		# Geosintetico
		elif self.subdominio[cell.index] == 5:
			values[0]=self.E[1]
		# Asafaltos
		else:
			values[0]=self.E[2]

class Nu(fnc.UserExpression):
	def __init__(self, subdominio, nu, **kwargs):
		super().__init__(**kwargs)
		self.subdominio = subdominio
		self.nu = nu

	def eval_cell(self, values, x, cell):
		# Fijo
		if self.subdominio[cell.index] == 3:
			values[0]=self.nu[0]
		# Movil
		elif self.subdominio[cell.index] == 4:
			values[0]=self.nu[0]
		# Geosintetico
		elif self.subdominio[cell.index] == 5:
			values[0]=self.nu[1]
		# Asafaltos
		else:
			values[0]=self.nu[2]

class Rho(fnc.UserExpression):
	def __init__(self, subdominio, rho, **kwargs):
		super().__init__(**kwargs)
		self.subdominio = subdominio
		self.rho = rho

	def eval_cell(self, values, x, cell):
		# Fijo
		if self.subdominio[cell.index] == 3:
			values[0]=self.rho[0]
		# Movil
		elif self.subdominio[cell.index] == 4:
			values[0]=self.rho[0]
		# Geosintetico
		elif self.subdominio[cell.index] == 5:
			values[0]=self.rho[1]
		# Asafaltos
		else:
			values[0]=self.rho[2]

# Parámetros de tiempo
T      = 2.0
Nsteps = 10
dt     = fnc.Constant(T/Nsteps)

# Carga aplicada
p0 = 1000.
#cutoff_Tc = T/4
#p = fnc.Expression(("t <= tc ? p0*t/tc : p0", "0", "0"), t=0, tc=cutoff_Tc, p0=p0, degree=0)
cutoff_Tc = T
p = fnc.Expression(("t <= tc ? p0*t/tc : p0", "0", "0"), t=0, tc=cutoff_Tc, p0=p0, degree=0)

# Defino las propiedades de los materiales
# Homigon
E_h= 1.0 
nu_h = 0.2
rho_h = 0.0
# Geosintetico
E_g= 1.0 
nu_g = 0.2
rho_g = 0.0
# Asfalto
E_a= 1.0 
nu_a = 0.2
rho_a = 0.0

# metemos todo en listas
listE	=[E_h,		E_g,	E_a]
listNu	=[nu_h,		nu_g,	nu_a]
listRho	=[rho_h,	rho_g,	rho_a]

e  	= E(subdominio,		listE,		degree=1)
nu 	= Nu(subdominio,	listNu,		degree=1)
rho = Rho(subdominio,	listRho,	degree=1)

mu = e/(2.0*(1.0+nu))
lmbda = e*nu/((1.0+nu)*(1.0-2.0*nu))
lmbda = 2*mu*lmbda/(lmbda+2*mu)

#Coeficientes de amortiguamiento de Rayleigh
eta_m = fnc.Constant(0.01)
eta_k = fnc.Constant(0.05)

#------------------------
#Discretización en tiempo

#Modelos alfa generalizados
alpha_m = fnc.Constant(0.2)
alpha_f = fnc.Constant(0.4)
gamma   = fnc.Constant(0.5+alpha_f-alpha_m)
beta    = fnc.Constant((gamma+0.5)**2/4.)


# Se define igual que siempre
V = fnc.VectorFunctionSpace(mesh, "CG", 1)
# Acá es distinto, definimos un tensor!
Vsig = fnc.TensorFunctionSpace(mesh, "DG", 0)

# Funciones Test y trial
du = fnc.TrialFunction(V)
w =  fnc.TestFunction(V)
# Desplazamiento actual Desconocido!
u =  fnc.Function(V, name="Desplazamiento")
#Campos del paso anterior (Desplazamiento, velocidad, aceleración)
u_old = fnc.Function(V)
v_old = fnc.Function(V)
a_old = fnc.Function(V)

#Define la medida de la integral de borde (contorno)
dss = fnc.ds(subdomain_data=bordes)

fnc.dx=fnc.dx(metadata={'quadrature_degree': 3})

#Condición de borde del costado izquierdo
zero = fnc.Constant((0.0, 0.0, 0.0))
bc1 = fnc.DirichletBC(V, zero, bordes, 1)
#bc2 = fnc.DirichletBC(V, zero, bordes, 9)
#bc3 = fnc.DirichletBC(V, zero, bordes, 10)
bc  = [bc1, bc1, bc1]

def eps(v):
    return fnc.sym(fnc.grad(v))

# Tensor de tensiones
def sigma(r):
    return 2.0*mu*fnc.sym(fnc.grad(r)) + lmbda*fnc.tr(fnc.sym(fnc.grad(r)))*fnc.Identity(len(r))

# Masa
def m(u,w):
    return rho*fnc.inner(u, w)*fnc.dx

#Rigidez
def k(u, w):
    return fnc.inner(sigma(u), fnc.sym(fnc.grad(w)))*fnc.dx

#Amortiguamiento de Rayleigh
def c(u, w):
    return eta_m*m(u, w) + eta_k*k(u, w)

#Trabajo de fuerzas externas
def Wext(w):
    return fnc.dot(w, p)*dss(2)# + fnc.dot(w, p)*dss(8) + fnc.dot(w, p)*dss(1)

#Actualiza la aceleración
# a = 1/(2*beta)*((u - u0 - v0*dt)/(0.5*dt*dt) - (1-2*beta)*a0)
def update_a(u, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        beta_ = beta
    else:
        dt_ = float(dt)
        beta_ = float(beta)
    return (u-u_old-dt_*v_old)/beta_/dt_**2 - (1-2*beta_)/2/beta_*a_old

#Acutaliza la velocidad
# v = dt * ((1-gamma)*a0 + gamma*a) + v0
def update_v(a, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        gamma_ = gamma
    else:
        dt_ = float(dt)
        gamma_ = float(gamma)
    return v_old + dt_*((1-gamma_)*a_old + gamma_*a)

def update_fields(u, u_old, v_old, a_old):
    """Actualiza los campos al final de cada paso.""" 

    # Get vectors (references)
    u_vec, u0_vec  = u.vector(), u_old.vector()
    v0_vec, a0_vec = v_old.vector(), a_old.vector() 

    # use update functions using vector arguments
    a_vec = update_a(u_vec, u0_vec, v0_vec, a0_vec, ufl=False)
    v_vec = update_v(a_vec, u0_vec, v0_vec, a0_vec, ufl=False)

    # Update (u_old <- u)
    v_old.vector()[:], a_old.vector()[:] = v_vec, a_vec
    u_old.vector()[:] = u.vector()

#Tiempos intermedio entre tn y tn+1 utilizando modelo alfa generalizado
def avg(x_old, x_new, alpha):
	return alpha*x_old + (1-alpha)*x_new

#Formulación variacional
a_new = update_a(du, u_old, v_old, a_old, ufl=True)
v_new = update_v(a_new, u_old, v_old, a_old, ufl=True)
res = m(avg(a_old, a_new, alpha_m), w) + c(avg(v_old, v_new, alpha_f), w) + k(avg(u_old, du, alpha_f), w) - Wext(w)
a_form = fnc.lhs(res)
L_form = fnc.rhs(res)


#Solvers (veremos luego)
K, res = fnc.assemble_system(a_form, L_form, bc)
solver = fnc.LUSolver(K, 'default')#"mumps")
solver.parameters["symmetric"] = True

# We now initiate the time stepping loop. We will keep track of the beam vertical tip
# displacement over time as well as the different parts of the system total energy. We
# will also compute the stress field and save it, along with the displacement field, in
# a ``XDMFFile``.  
# The option `flush_ouput` enables to open the result file before the loop is finished,
# the ``function_share_mesh`` option tells that only one mesh is used for all functions
# of a given time step (displacement and stress) while the ``rewrite_function_mesh`` enforces
# that the same mesh is used for all time steps. These two options enables writing the mesh
# information only once instead of :math:`2N_{steps}` times::

# Time-stepping
time = np.linspace(0, T, Nsteps+1)
u_tip = np.zeros((Nsteps+1,))
energies = np.zeros((Nsteps+1, 4))
E_damp = 0
E_ext = 0
sig = fnc.Function(Vsig, name="sigma")
xdmf_file = fnc.XDMFFile("%s.xdmf" %(archivo))
xdmf_file.parameters["flush_output"] = True
xdmf_file.parameters["functions_share_mesh"] = True
xdmf_file.parameters["rewrite_function_mesh"] = False
ofile = fnc.File("%s.pvd" %(archivovtu))


# The time loop is now started, the loading is first evaluated at :math:`t=t_{n+1-\alpha_f}`. The
# corresponding system right-hand side is then assembled and the system is solved. The different
# fields are then updated with the newly computed quantities. Finally, some post-processing is
# performed: stresses are computed and written to the result file and the tip displacement and
# the different energies are recorded:

def local_project(v, V, u=None):
    """Element-wise projection using LocalSolver"""
    dv = fnc.TrialFunction(V)
    v_ = fnc.TestFunction(V)
    a_proj = fnc.inner(dv, v_)*fnc.dx
    b_proj = fnc.inner(v, v_)*fnc.dx
    solver = fnc.LocalSolver(a_proj, b_proj)
    solver.factorize()
    if u is None:
        u = fnc.Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return

for (i, dt) in enumerate(np.diff(time)):

    t = time[i+1]
    print("Time: ", t)

    # Forces are evaluated at t_{n+1-alpha_f}=t_{n+1}-alpha_f*dt
    p.t = t-float(alpha_f*dt)

    # Solve for new displacement
    res = fnc.assemble(L_form)
    bc[0].apply(res)
    bc[1].apply(res)
    bc[2].apply(res)
    solver.solve(K, u.vector(), res)

    # Update old fields with new quantities
    update_fields(u, u_old, v_old, a_old)

    # Save solution to XDMF format
    xdmf_file.write(u, t)

    # Compute stresses and save to file
    local_project(sigma(u), Vsig, sig)
    xdmf_file.write(sig, t)
    
    ofile << sig, t

    p.t = t
    #Calculo de energias
    if fnc.MPI.comm_world.size == 1:
        u_tip[i+1] = u(1., 0.05, 0.)[1]
    E_elas = fnc.assemble(0.5*k(u_old, u_old))
    E_kin = fnc.assemble(0.5*m(v_old, v_old))
    E_damp += dt*fnc.assemble(c(v_old, v_old))
    # E_ext += assemble(Wext(u-u_old))
    E_tot = E_elas+E_kin+E_damp #-E_ext
    energies[i+1, :] = np.array([E_elas, E_kin, E_damp, E_tot])
