import geometria as g
import fenics as fnc
import os

# Dimensiones
l, h_prob, b = 30.0, 5.0, 10.0  # largo, alto, ancho del asfalto
e_geo = 0.4 # Espesor del geosintetico
h_base = 2.0  # alto base de hormigon
e = 0.5 # espesor de la ranura
d = 1 # ancho de la zonad de detalle

resmax=1.0 # tamaño de los elementos de los extremos
resmin=0.5 # tamaño de los elementos del centro

archivo = 'calc/modelo' # PATH/AL/ARCHIVO

parametros=[l, h_prob, b, e_geo, h_base, e, d, resmax, resmin, archivo]

# Mostrar mallado
malla=False
# Mostrar geometria
geometria=True

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

if file1==file2:
	print('parametros anteriores = %s' %(file2))
	print('parametros actuales   = %s' %(file1))
	print('Al ser iguales no se regenerla la malla')
else:
	print('parametros anteriores = %s' %(file2))
	print('parametros actuales   = %s' %(file1))
	print('Se regenerla la malla')
	# llamo a la funcion para regenerar el calculo del mashado completo
	g.parametrosProbeta(parametros, malla, geometria)

# Obtenemos la malla creada en el archivo geometria.py
mesh = fnc.Mesh("%s.xml" %(archivo)) # mesh
subdominio = fnc.MeshFunction('size_t', mesh, "%s_physical_region.xml" %(archivo)) # Volumenes
bordes = fnc.MeshFunction('size_t', mesh, "%s_facet_region.xml" %(archivo)) # Bordes

# Creamos la Funcion Espacio
V=fnc.FunctionSpace(mesh,'P', 1)

# Defino las Condiciones de borde
bhf= fnc.DirichletBC(V, fnc.Constant(5.0), bordes, 1)
bhm= fnc.DirichletBC(V, fnc.Constant(0.0), bordes, 2)
bc=[bhf, bhm]

# Defino los Volumenes

class PropiedadesMaterial(fnc.UserExpression):  # Ahora anda con python3, hay que usar UserExpression
    def __init__(self, subdomains, E, nu, **kwargs):
        super().__init__(**kwargs)  # Ahora anda con python3
        self.subdomains = subdomains
        self.E = E
        self.nu = nu

    def eval_cell(self, values, cell):
        # Fijo
        if self.subdomains[cell.index] == 3:
            values[0] = self.E[0]
        # Movil
        elif self.subdomains[cell.idex] == 4:
            values[1] = self.E[0]
        # Geosintetico
        elif self.subdomains[cell.idex] == 5:
            values[2] = self.E[1]
        # Asafaltos
        else:
            values[3] = self.E[2]