import meshio
import pygmsh
import numpy as np
import os

# Dimensiones
# Asfalto
l,h_prob,b = 30.0, 5.0, 10.0  # largo, alto, ancho
# Geosintetico
e_geo=0.1 # Espesor del geosintetico
# Base de hormigon o asfalto
h_base = 2.0  # alto
# Ranura
e=0.5 # espesor

elem=0.5 # tamaño de los elementos

nombre_archivo = 'calc/modelo' # PATH/AL/ARCHIVO

# Dimensiones de los elementos
geo = pygmsh.built_in.Geometry()

b=[0.0,b,0.0]

# Bases
base_movil = geo.add_polygon([
	[0.0,    0.0,0.0   ],
	[(l-e)/2,0.0,0.0   ],
	[(l-e)/2,0.0,h_base],
	[0.0,    0.0,h_base],
	],
	lcar=elem)
	
base_fija = geo.add_polygon([
	[(l+e)/2,0.0,0.0   ],
	[l,      0.0,0.0   ],
	[l,      0.0,h_base],
	[(l+e)/2,0.0,h_base],
	],
	lcar=elem)

# Geosintetico
geosintetico = geo.add_polygon([
	[0.0,0.0,h_base      ],
	[l,  0.0,h_base      ],
	[l,  0.0,h_base+e_geo],
	[0.0,0.0,h_base+e_geo],
	],
	lcar=elem)

# Probeta de asfalto
probeta = geo.add_polygon([
	[0.0,0.0,h_base+e_geo       ],
	[l,  0.0,h_base+e_geo       ],
	[l,  0.0,h_base+e_geo+h_prob],
	[0.0,0.0,h_base+e_geo+h_prob],
	],
    lcar=elem)

# Extrusion de la geometria
geo.extrude(base_movil,  b)
geo.extrude(base_fija,   b)
geo.extrude(geosintetico,b)
geo.extrude(probeta,     b)

# Escribimos el archivo .geo
a=open('%s.geo' %(nombre_archivo), 'w')
a.write(geo.get_code())
a.close()

# Mashado
mesh = pygmsh.generate_mesh(geo)

meshio.write('%s.xml' %(nombre_archivo), mesh)
meshio.write('%s.vtk' %(nombre_archivo), mesh)

# Mostrar el resultado en gmsh
com = 'gmsh %s.geo' %(nombre_archivo) 
os.system(com)