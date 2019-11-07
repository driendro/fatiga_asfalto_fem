#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import meshio
import pygmsh as pg
import os

def probeta(x, malla=False, geometria=False, *args, **kwargs):
    
	l, h_prob, b = x[0], x[1], x[2]  # largo, alto, ancho del asfalto
	e_geo = x[3] # Espesor del geosintetico
	h_base = x[4]  # alto base de hormigon
	e = x[5] # espesor de la ranura
	d = x[6] # ancho de la zonad de detalle

	resmax=x[7] # tamaño de los elementos de los extremos
	resmin=x[8] # tamaño de los elementos del centro

	nombre_archivo = x[9] # PATH/AL/ARCHIVO

	# Dimensiones de los elementos
	geo = pg.built_in.Geometry()

	# Coordenadas de la probeta
	x0=0
	x1=(l-e)/2
	x2=(l-e-d)/2
	x3=(l+e+d)/2
	x4=(l+e)/2
	x5=l

	y0=0
	y1=b

	z0=0
	z1=h_base
	z2=z1+e_geo
	z3=z2+h_prob

	# Puntos de la Geometria
	# Fondo
	p1=  geo.add_point([x0,y0,z0], resmax)
	p2=  geo.add_point([x0,y1,z0], resmax)
	p3=  geo.add_point([x1,y0,z0], resmax)
	p4=  geo.add_point([x1,y1,z0], resmax)
	p9=  geo.add_point([x4,y0,z0], resmax)
	p10= geo.add_point([x4,y1,z0], resmax)
	p11= geo.add_point([x5,y0,z0], resmax)
	p12= geo.add_point([x5,y1,z0], resmax)
	# Geosintetico y hormigon
	p101= geo.add_point([x0,y0,z1], resmax)
	p102= geo.add_point([x0,y1,z1], resmax)
	p103= geo.add_point([x1,y0,z1], resmax)
	p104= geo.add_point([x1,y1,z1], resmax)
	p105= geo.add_point([x2,y0,z1], resmin)
	p106= geo.add_point([x2,y1,z1], resmin)
	p107= geo.add_point([x3,y0,z1], resmin)
	p108= geo.add_point([x3,y1,z1], resmin)
	p109= geo.add_point([x4,y0,z1], resmax)
	p110= geo.add_point([x4,y1,z1], resmax)
	p111= geo.add_point([x5,y0,z1], resmax)
	p112= geo.add_point([x5,y1,z1], resmax)
	# Geosintetico y Asfalto
	p201= geo.add_point([x0,y0,z2], resmax)
	p202= geo.add_point([x0,y1,z2], resmax)
	p205= geo.add_point([x2,y0,z2], resmin)
	p206= geo.add_point([x2,y1,z2], resmin)
	p207= geo.add_point([x3,y0,z2], resmin)
	p208= geo.add_point([x3,y1,z2], resmin)
	p211= geo.add_point([x5,y0,z2], resmax)
	p212= geo.add_point([x5,y1,z2], resmax)
	# Asfalto superior
	p301= geo.add_point([x0,y0,z3], resmax)
	p302= geo.add_point([x0,y1,z3], resmax)
	p305= geo.add_point([x2,y0,z3], resmin)
	p306= geo.add_point([x2,y1,z3], resmin)
	p307= geo.add_point([x3,y0,z3], resmin)
	p308= geo.add_point([x3,y1,z3], resmin)
	p311= geo.add_point([x5,y0,z3], resmax)
	p312= geo.add_point([x5,y1,z3], resmax)

	# Lineas que conforman la geometria
	# Lineas Horizontales de la cara frontal
	l101= geo.add_line(p1,p3)
	l102= geo.add_line(p9,p11)
	l103= geo.add_line(p101,p103)
	l104= geo.add_line(p109,p111)
	l105= geo.add_line(p101,p105)
	l106= geo.add_line(p105,p107)
	l107= geo.add_line(p107,p111)
	l108= geo.add_line(p201,p205)
	l109= geo.add_line(p205,p207)
	l110= geo.add_line(p207,p211)
	l111= geo.add_line(p301,p305)
	l112= geo.add_line(p305,p307)
	l113= geo.add_line(p307,p311)

	# Lineas Verticales de la Cara Frontal
	lv111= geo.add_line(p1,p101)
	lv112= geo.add_line(p101,p201)
	lv113= geo.add_line(p201,p301)
	lv122= geo.add_line(p105,p205)
	lv123= geo.add_line(p205,p305)
	lv131= geo.add_line(p3,p103)
	lv141= geo.add_line(p9,p109)
	lv152= geo.add_line(p107,p207)
	lv153= geo.add_line(p207,p307)
	lv161= geo.add_line(p11,p111)
	lv162= geo.add_line(p111,p211)
	lv163= geo.add_line(p211,p311)

	# Lineas Horizontales de la cara trasera
	l201= geo.add_line(p2,p4)
	l202= geo.add_line(p10,p12)
	l203= geo.add_line(p102,p104)
	l204= geo.add_line(p110,p112)
	l205= geo.add_line(p102,p106)
	l206= geo.add_line(p106,p108)
	l207= geo.add_line(p108,p112)
	l208= geo.add_line(p202,p206)
	l209= geo.add_line(p206,p208)
	l210= geo.add_line(p208,p212)
	l211= geo.add_line(p302,p306)
	l212= geo.add_line(p306,p308)
	l213= geo.add_line(p308,p312)

	# Lineas Verticales de la cara trasera
	lv211= geo.add_line(p2,p102)
	lv212= geo.add_line(p102,p202)
	lv213= geo.add_line(p202,p302)
	lv222= geo.add_line(p106,p206)
	lv223= geo.add_line(p206,p306)
	lv231= geo.add_line(p4,p104)
	lv241= geo.add_line(p10,p110)
	lv252= geo.add_line(p108,p208)
	lv253= geo.add_line(p208,p308)
	lv261= geo.add_line(p12,p112)
	lv262= geo.add_line(p112,p212)
	lv263= geo.add_line(p212,p312)

	# Lineas Transeverales
	l11= geo.add_line(p1,p2)
	l12= geo.add_line(p101,p102)
	l13= geo.add_line(p201,p202)
	l14= geo.add_line(p301,p302)
	l22= geo.add_line(p105,p106)
	l23= geo.add_line(p205,p206)
	l24= geo.add_line(p305,p306)
	l31= geo.add_line(p3,p4)
	l32= geo.add_line(p103,p104)
	l41= geo.add_line(p9,p10)
	l42= geo.add_line(p109,p110)
	l52= geo.add_line(p107,p108)
	l53= geo.add_line(p207,p208)
	l54= geo.add_line(p307,p308)
	l61= geo.add_line(p11,p12)
	l62= geo.add_line(p111,p112)
	l63= geo.add_line(p211,p212)
	l64= geo.add_line(p311,p312)

	# Superficies de la geometria y Volumenes
	# Base de Hormigon Mobil
	# Inferior
	bhmin= geo.add_line_loop([l101,l31,-l201,-l11])
	shmin= geo.add_surface(bhmin)
	# Superior
	bhms= geo.add_line_loop([l103,l32,-l203,-l12])
	shms= geo.add_surface(bhms)
	# Frontal
	bhmf= geo.add_line_loop([l101,lv131,-l103,-lv111])
	shmf= geo.add_surface(bhmf)
	# Trasera
	bhmt= geo.add_line_loop([l201,lv231,-l203,-lv211])
	shmt= geo.add_surface(bhmt)
	# Derecha
	bhmd= geo.add_line_loop([l11,lv211,-l12,-lv111])
	shmd= geo.add_surface(bhmd)
	# Izquierda
	bhmiz= geo.add_line_loop([l31,lv231,-l32,-lv131])
	shmiz= geo.add_surface(bhmiz)
	# Volumen de la Base movil
	sbhm=geo.add_surface_loop([shmin, shms, shmf, shmt, shmd, shmiz])
	vbhm=geo.add_volume(sbhm)

	# Base de Hormigon Fija
	# Inferior
	bhfin= geo.add_line_loop([l102,l61,-l202,-l41])
	shfin= geo.add_surface(bhfin)
	# Superior
	bhfs= geo.add_line_loop([l104,l62,-l204,-l42])
	shfs= geo.add_surface(bhfs)
	# Frontal
	bhff= geo.add_line_loop([l102,lv161,-l104,-lv141])
	shff= geo.add_surface(bhff)
	# Trasera
	bhft= geo.add_line_loop([l202,lv261,-l204,-lv241])
	shft= geo.add_surface(bhft)
	# Derecha
	bhfd= geo.add_line_loop([l61,lv261,-l62,-lv161])
	shfd= geo.add_surface(bhfd)
	# Izquierda
	bhfiz= geo.add_line_loop([l41,lv241,-l42,-lv141])
	shfiz= geo.add_surface(bhfiz)
	# Volumen de la Base Fija
	sbhf=geo.add_surface_loop([shfin, shfs, shff, shft, shfd, shfiz])
	vbhf=geo.add_volume(sbhf)

	# Geosintetico
	# InferiorD
	bgsind= geo.add_line_loop([l107,l62,-l207,-l52])
	sgsind= geo.add_surface(bgsind)
	# InferiorM
	bgsinm= geo.add_line_loop([l106,l52,-l206,-l22])
	sgsinm= geo.add_surface(bgsinm)
	# InferiorI
	bgsini= geo.add_line_loop([l105,l22,-l205,-l12])
	sgsini= geo.add_surface(bgsini)
	# SuperiorD
	bgssd= geo.add_line_loop([l110,l63,-l210,-l53])
	sgssd= geo.add_surface(bgssd)
	# SuperiorM
	bgssm= geo.add_line_loop([l109,l53,-l209,-l23])
	sgssm= geo.add_surface(bgssm)
	# SuperiorI
	bgssi= geo.add_line_loop([l108,l23,-l208,-l13])
	sgssi= geo.add_surface(bgssi)
	# FrontalD
	bgsfd= geo.add_line_loop([l107,lv162,-l110,-lv152])
	sgsfd= geo.add_surface(bgsfd)
	# FrontalM
	bgsfm= geo.add_line_loop([l106,lv152,-l109,-lv122])
	sgsfm= geo.add_surface(bgsfm)
	# FrontalI
	bgsfi= geo.add_line_loop([l105,lv122,-l108,-lv112])
	sgsfi= geo.add_surface(bgsfi)
	# TraseraD
	bgstd= geo.add_line_loop([l207,lv262,-l210,-lv252])
	sgstd= geo.add_surface(bgstd)
	# TraseraM
	bgstm= geo.add_line_loop([l206,lv252,-l209,-lv222])
	sgstm= geo.add_surface(bgstm)
	# TraseraI
	bgsti= geo.add_line_loop([l205,lv222,-l208,-lv212])
	sgsti= geo.add_surface(bgsti)
	# Derecha
	bgsd= geo.add_line_loop([l62,lv262,-l63,-lv162])
	sgsd= geo.add_surface(bgsd)
	# Izquierda
	bgsiz= geo.add_line_loop([l12,lv212,-l13,-lv112])
	sgsiz= geo.add_surface(bgsiz)
	# Volumen del Geosintetico
	sgs=geo.add_surface_loop([sgsind, sgsinm, sgsini,
							sgssd, sgssm, sgssi,
							sgsfd, sgsfm, sgsfi,
							sgstd, sgstm, sgsti,
							sgsd,
							sgsiz]
							)
	vgs=geo.add_volume(sgs)

	# Probeta Asfalto
	# InferiorD
	bpaind= geo.add_line_loop([l110,l63,-l210,-l53])
	spaind= geo.add_surface(bpaind)
	# InferiorM
	bpainm= geo.add_line_loop([l109,l53,-l209,-l23])
	spainm= geo.add_surface(bpainm)
	# InferiorI
	bpaini= geo.add_line_loop([l108,l23,-l208,-l13])
	spaini= geo.add_surface(bpaini)
	# SuperiorD
	bpasd= geo.add_line_loop([l113,l64,-l213,-l54])
	spasd= geo.add_surface(bpasd)
	# SuperiorM
	bpasm= geo.add_line_loop([l112,l54,-l212,-l24])
	spasm= geo.add_surface(bpasm)
	# SuperiorI
	bpasi= geo.add_line_loop([l111,l24,-l211,-l14])
	spasi= geo.add_surface(bpasi)
	# FrontalD
	bpafd= geo.add_line_loop([l110,lv163,-l113,-lv153])
	spafd= geo.add_surface(bpafd)
	# FrontalM
	bpafm= geo.add_line_loop([l109,lv153,-l112,-lv123])
	spafm= geo.add_surface(bpafm)
	# FrontalI
	bpafi= geo.add_line_loop([l108,lv123,-l111,-lv113])
	spafi= geo.add_surface(bpafi)
	# TraseraD
	bpatd= geo.add_line_loop([l210,lv263,-l213,-lv253])
	spatd= geo.add_surface(bpatd)
	# TraseraM
	bpatm= geo.add_line_loop([l209,lv253,-l212,-lv223])
	spatm= geo.add_surface(bpatm)
	# TraseraI
	bpati= geo.add_line_loop([l208,lv223,-l211,-lv213])
	spati= geo.add_surface(bpati)
	# Derecha
	bpad= geo.add_line_loop([l63,lv263,-l64,-lv163])
	spad= geo.add_surface(bpad)
	# Izquierda
	bpaiz= geo.add_line_loop([l13,lv213,-l14,-lv113])
	spaiz= geo.add_surface(bpaiz)
	# Volumen del Probeta Asfalto
	spa=geo.add_surface_loop([spaind, spainm, spaini,
							spasd, spasm, spasi,
							spafd, spafm, spafi,
							spatd, spatm, spati,
							spad,
							spaiz,
							])
	vpa=geo.add_volume(spa)

	# Superficies Fisicas
	geo.add_physical(shmin, label=1)
	geo.add_physical(shfin, label=2)

	# Volumenes Fisicos
	geo.add_physical(vbhm, label=3)
	geo.add_physical(vbhf, label=4)
	geo.add_physical(vgs, label=5)
	geo.add_physical(vpa, label=6)

	# Escribimos el archivo .geo
	a=open('%s.geo' %(nombre_archivo), 'w')
	a.write(geo.get_code())
	a.close()

	# Realizar el mashado y el archivo xml
	s= 'gmsh -3  %s.geo -format msh2' %(nombre_archivo)
	os.system(s)
	s= 'dolfin-convert %s.msh %s.xml' %(nombre_archivo, nombre_archivo)
	os.system(s)


	if geometria:
		# Mostrar el resultado en gmsh
		s = 'gmsh %s.geo' %(nombre_archivo) 
		os.system(s)
	
	if malla:
		# Mostrar el resultado en gmsh
		s = 'gmsh %s.msh' %(nombre_archivo) 
		os.system(s)
