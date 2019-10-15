#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import meshio
import pygmsh as pg
import numpy as np
import os

# Dimensiones
l,h_prob,b = 30.0, 5.0, 10.0  # largo, alto, ancho del asfalto
e_geo=0.1 # Espesor del geosintetico
h_base = 2.0  # alto base de hormigon
e=0.5 # espesor de la ranura



resmax=1.0 # tamaño de los elementos
resmin=0.1 

nombre_archivo = 'calc/modelo' # PATH/AL/ARCHIVO

# Dimensiones de los elementos
geo = pg.built_in.Geometry()

# Coordenadas de la probeta
x0=0
x1=(l-e)/2
x2=(l+e)/2
x3=l

y0=0
y1=b

z0=0
z1=h_base
z2=z1+e_geo
z3=z2+h_prob

# Puntos de la Geometria
# Fondo
p1= geo.add_point([x0,y0,z0], resmax)
p2= geo.add_point([x0,y1,z0], resmax)
p3= geo.add_point([x1,y0,z0], resmax)
p4= geo.add_point([x1,y1,z0], resmax)
p5= geo.add_point([x2,y0,z0], resmax)
p6= geo.add_point([x2,y1,z0], resmax)
p7= geo.add_point([x3,y0,z0], resmax)
p8= geo.add_point([x3,y1,z0], resmax)
# Geosintetico y hormigon
p101= geo.add_point([x0,y0,z1], resmax)
p102= geo.add_point([x0,y1,z1], resmax)
p103= geo.add_point([x1,y0,z1], resmax)
p104= geo.add_point([x1,y1,z1], resmax)
p105= geo.add_point([x2,y0,z1], resmax)
p106= geo.add_point([x2,y1,z1], resmax)
p107= geo.add_point([x3,y0,z1], resmax)
p108= geo.add_point([x3,y1,z1], resmax)
# Geosintetico y Asfalto
p201= geo.add_point([x0,y0,z2], resmax)
p202= geo.add_point([x0,y1,z2], resmax)
p207= geo.add_point([x3,y0,z2], resmax)
p208= geo.add_point([x3,y1,z2], resmax)
# Asfalto superior
p301= geo.add_point([x0,y0,z3], resmax)
p302= geo.add_point([x0,y1,z3], resmax)
p307= geo.add_point([x3,y0,z3], resmax)
p308= geo.add_point([x3,y1,z3], resmax)

# Lineas que conforman la geometria
# Lineas Horizontales de la cara frontal
l101= geo.add_line(p1,p3)
l102= geo.add_line(p5,p7)
l103= geo.add_line(p101,p103)
l104= geo.add_line(p105,p107)
l105= geo.add_line(p101,p107)
l106= geo.add_line(p201,p207)
l107= geo.add_line(p301,p307)

# Lineas Verticales de la Cara Frontal
l111= geo.add_line(p1,p101)
l112= geo.add_line(p101,p201)
l113= geo.add_line(p201,p301)
l121= geo.add_line(p3,p103)
l131= geo.add_line(p5,p105)
l141= geo.add_line(p7,p107)
l142= geo.add_line(p107,p207)
l143= geo.add_line(p207,p307)

# Lineas Horizontales de la cara trasera
l201= geo.add_line(p2,p4)
l202= geo.add_line(p6,p8)
l203= geo.add_line(p102,p104)
l204= geo.add_line(p106,p108)
l205= geo.add_line(p102,p108)
l206= geo.add_line(p202,p208)
l207= geo.add_line(p302,p308)

# Lineas Verticales de la cara trasera
l211= geo.add_line(p2,p102)
l212= geo.add_line(p102,p202)
l213= geo.add_line(p202,p302)
l221= geo.add_line(p4,p104)
l231= geo.add_line(p6,p106)
l241= geo.add_line(p8,p108)
l242= geo.add_line(p108,p208)
l243= geo.add_line(p208,p308)

# Lineas Transeverales
l11= geo.add_line(p1,p2)
l12= geo.add_line(p101,p102)
l13= geo.add_line(p201,p202)
l14= geo.add_line(p301,p302)
l21= geo.add_line(p3,p4)
l22= geo.add_line(p103,p104)
l31= geo.add_line(p5,p6)
l32= geo.add_line(p105,p106)
l41= geo.add_line(p7,p8)
l42= geo.add_line(p107,p108)
l43= geo.add_line(p207,p208)
l44= geo.add_line(p307,p308)

# Superficies de la geometria y Volumenes
# Base de Hormigon Mobil
# Inferior
bhmin= geo.add_line_loop([l101,l21,-l201,-l11])
shmin= geo.add_surface(bhmin)
# Superior
bhms= geo.add_line_loop([l103,l22,-l203,-l12])
shms= geo.add_surface(bhms)
# Frontal
bhmf= geo.add_line_loop([l101,l121,-l103,-l111])
shmf= geo.add_surface(bhmf)
# Trasera
bhmt= geo.add_line_loop([l201,l221,-l203,-l211])
shmt= geo.add_surface(bhmt)
# Derecha
bhmd= geo.add_line_loop([-l21,l121,l22,-l221])
shmd= geo.add_surface(bhmd)
# Izquierda
bhmiz= geo.add_line_loop([-l11,l111,l12,-l211])
shmiz= geo.add_surface(bhmiz)
# Volumen de la Base movil
sbhm=geo.add_surface_loop([shmin, shms, shmf, shmt, shmd, shmiz])
vbhm=geo.add_volume(sbhm)

# Base de Hormigon Fija
# Inferior
bhfin= geo.add_line_loop([l102,l41,-l202,-l31])
shfin= geo.add_surface(bhfin)
# Superior
bhfs= geo.add_line_loop([l104,l42,-l204,-l32])
shfs= geo.add_surface(bhfs)
# Frontal
bhff= geo.add_line_loop([l102,l141,-l104,-l131])
shff= geo.add_surface(bhff)
# Trasera
bhft= geo.add_line_loop([l202,l241,-l204,-l231])
shft= geo.add_surface(bhft)
# Derecha
bhfd= geo.add_line_loop([-l41,l141,l42,-l241])
shfd= geo.add_surface(bhfd)
# Izquierda
bhfiz= geo.add_line_loop([-l31,l131,l32,-l231])
shfiz= geo.add_surface(bhfiz)
# Volumen de la Base Fija
sbhf=geo.add_surface_loop([shfin, shfs, shff, shft, shfd, shfiz])
vbhf=geo.add_volume(sbhf)

# Geosintetico
# Inferior
bgsin= geo.add_line_loop([-l12,l105,l42,-l205])
sgsin= geo.add_surface(bgsin)
# Superior
bgss= geo.add_line_loop([-l13,l106,l43,-l206])
sgss= geo.add_surface(bgss)
# Frontal
bgsf= geo.add_line_loop([l105,l142,-l106,-l112])
sgsf= geo.add_surface(bgsf)
# Trasera
bgst= geo.add_line_loop([l205,l242,-l206,-l212])
sgst= geo.add_surface(bgst)
# Derecha
bgsd= geo.add_line_loop([-l42,l142,l43,-l242])
sgsd= geo.add_surface(bgsd)
# Izquierda
bgsiz= geo.add_line_loop([-l12,l112,l13,-l212])
sgsiz= geo.add_surface(bgsiz)
# Volumen del Geosintetico
sgs=geo.add_surface_loop([sgsin, sgss, sgsf, sgst, sgsd, sgsiz])
vgs=geo.add_volume(sgs)

# Probeta Asfalto
# Inferior
bpain= geo.add_line_loop([l106,l43,-l206,-l13])
spain= geo.add_surface(bpain)
# Superior
bpas= geo.add_line_loop([l107,l44,-l207,-l14])
spas= geo.add_surface(bpas)
# Frontal
bpaf= geo.add_line_loop([l106,l143,-l107,-l113])
spaf= geo.add_surface(bpaf)
# Trasera
bpat= geo.add_line_loop([l206,l243,-l207,-l213])
spat= geo.add_surface(bpat)
# Derecha
bpad= geo.add_line_loop([-l43,l143,l44,-l243])
spad= geo.add_surface(bpad)
# Izquierda
bpaiz= geo.add_line_loop([-l13,l113,l14,-l213])
spaiz= geo.add_surface(bpaiz)
# Volumen del Probeta Asfalto
spa=geo.add_surface_loop([spain, spas, spaf, spat, spad, spaiz])
vpa=geo.add_volume(spa)

# Superficies Fisicas
geo.add_physical_surface(shmin, label="supBaseMovil")
geo.add_physical_surface(shfin, label="supBaseFija")

# Volumenes Fisicos
geo.add_physical_volume(vbhm, label="volBaseMovil")
geo.add_physical_volume(vbhf, label="volBaseFija")
geo.add_physical_volume(vgs, label="volGeosintetico")
geo.add_physical_volume(vpa, label="volProbetaAsfalto")

# Escribimos el archivo .geo
a=open('%s.geo' %(nombre_archivo), 'w')
a.write(geo.get_code())
a.close()

# Realizar el mashado y el archivo xml
s= 'gmsh -3  %s.geo -format msh2' %(nombre_archivo)
os.system(s)
s= 'dolfin-convert %s.msh %s.xml' %(nombre_archivo, nombre_archivo)
os.system(s)

# Mostrar el resultado en gmsh
s = 'gmsh %s.msh' %(nombre_archivo) 
os.system(s)