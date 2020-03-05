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
	x0 = 0
	x1 = (l-e-d)/2
	x2 = (l-e)/2
	x3 = (l+e)/2
	x4 = (l+e+d)/2
	x5 = l
	
	y0 = 0
	y1 = b
	
	z0 = 0
	z1 = h_base
	z2 = z1+e_geo
	z3 = z2+h_prob
	
	
	# Puntos de la Geometria
	# Fondo
	p001 = geo.add_point([x0, y0, z0], resmax)
	p002 = geo.add_point([x0, y1, z0], resmax)
	p005 = geo.add_point([x2, y0, z0], resmax)
	p006 = geo.add_point([x2, y1, z0], resmax)
	p007 = geo.add_point([x3, y0, z0], resmax)
	p008 = geo.add_point([x3, y1, z0], resmax)
	p011 = geo.add_point([x5, y0, z0], resmax)
	p012 = geo.add_point([x5, y1, z0], resmax)
	# Geosintetico y hormigon
	p101 = geo.add_point([x0, y0, z1], resmax)
	p102 = geo.add_point([x0, y1, z1], resmax)
	p103 = geo.add_point([x1, y0, z1], resmin)
	p104 = geo.add_point([x1, y1, z1], resmin)
	p105 = geo.add_point([x2, y0, z1], resmin)
	p106 = geo.add_point([x2, y1, z1], resmin)
	p107 = geo.add_point([x3, y0, z1], resmin)
	p108 = geo.add_point([x3, y1, z1], resmin)
	p109 = geo.add_point([x4, y0, z1], resmin)
	p110 = geo.add_point([x4, y1, z1], resmin)
	p111 = geo.add_point([x5, y0, z1], resmax)
	p112 = geo.add_point([x5, y1, z1], resmax)
	# Geosintetico y Asfalto
	p201 = geo.add_point([x0, y0, z2], resmax)
	p202 = geo.add_point([x0, y1, z2], resmax)
	p203 = geo.add_point([x1, y0, z2], resmin)
	p204 = geo.add_point([x1, y1, z2], resmin)
	p209 = geo.add_point([x4, y0, z2], resmin)
	p210 = geo.add_point([x4, y1, z2], resmin)
	p211 = geo.add_point([x5, y0, z2], resmax)
	p212 = geo.add_point([x5, y1, z2], resmax)
	# Asfalto superior
	p301 = geo.add_point([x0, y0, z3], resmax)
	p302 = geo.add_point([x0, y1, z3], resmax)
	p303 = geo.add_point([x1, y0, z3], resmin)
	p304 = geo.add_point([x1, y1, z3], resmin)
	p309 = geo.add_point([x4, y0, z3], resmin)
	p310 = geo.add_point([x4, y1, z3], resmin)
	p311 = geo.add_point([x5, y0, z3], resmax)
	p312 = geo.add_point([x5, y1, z3], resmax)
	
	
	# Lineas que conforman la geometria
	# Lineas Horizontales de la cara frontal
	l101 = geo.add_line(p001, p005)
	l102 = geo.add_line(p007, p011)
	l103 = geo.add_line(p101, p103)
	l104 = geo.add_line(p103, p105)
	l105 = geo.add_line(p105, p107)
	l106 = geo.add_line(p107, p109)
	l107 = geo.add_line(p109, p111)
	l108 = geo.add_line(p201, p203)
	l109 = geo.add_line(p203, p209)
	l110 = geo.add_line(p209, p211)
	l111 = geo.add_line(p301, p303)
	l112 = geo.add_line(p303, p309)
	l113 = geo.add_line(p309, p311)
	
	# Lineas Verticales de la Cara Frontal
	lv111 = geo.add_line(p001, p101)
	lv112 = geo.add_line(p101, p201)
	lv113 = geo.add_line(p201, p301)
	lv122 = geo.add_line(p103, p203)
	lv123 = geo.add_line(p203, p303)
	lv131 = geo.add_line(p005, p105)
	lv141 = geo.add_line(p007, p107)
	lv152 = geo.add_line(p109, p209)
	lv153 = geo.add_line(p209, p309)
	lv161 = geo.add_line(p011, p111)
	lv162 = geo.add_line(p111, p211)
	lv163 = geo.add_line(p211, p311)
	
	# Lineas Horizontales de la cara trasera
	l201 = geo.add_line(p002, p006)
	l202 = geo.add_line(p008, p012)
	l203 = geo.add_line(p102, p104)
	l204 = geo.add_line(p104, p106)
	l205 = geo.add_line(p106, p108)
	l206 = geo.add_line(p108, p110)
	l207 = geo.add_line(p110, p112)
	l208 = geo.add_line(p202, p204)
	l209 = geo.add_line(p204, p210)
	l210 = geo.add_line(p210, p212)
	l211 = geo.add_line(p302, p304)
	l212 = geo.add_line(p304, p310)
	l213 = geo.add_line(p310, p312)
	
	# Lineas Verticales de la cara trasera
	lv211 = geo.add_line(p002, p102)
	lv212 = geo.add_line(p102, p202)
	lv213 = geo.add_line(p202, p302)
	lv222 = geo.add_line(p104, p204)
	lv223 = geo.add_line(p204, p304)
	lv231 = geo.add_line(p006, p106)
	lv241 = geo.add_line(p008, p108)
	lv252 = geo.add_line(p110, p210)
	lv253 = geo.add_line(p210, p310)
	lv261 = geo.add_line(p012, p112)
	lv262 = geo.add_line(p112, p212)
	lv263 = geo.add_line(p212, p312)
	
	# Lineas Transeverales
	l11 = geo.add_line(p001, p002)
	l12 = geo.add_line(p101, p102)
	l13 = geo.add_line(p201, p202)
	l14 = geo.add_line(p301, p302)
	l22 = geo.add_line(p103, p104)
	l23 = geo.add_line(p203, p204)
	l24 = geo.add_line(p303, p304)
	l31 = geo.add_line(p005, p006)
	l32 = geo.add_line(p105, p106)
	l41 = geo.add_line(p007, p008)
	l42 = geo.add_line(p107, p108)
	l52 = geo.add_line(p109, p110)
	l53 = geo.add_line(p209, p210)
	l54 = geo.add_line(p309, p310)
	l61 = geo.add_line(p011, p012)
	l62 = geo.add_line(p111, p112)
	l63 = geo.add_line(p211, p212)
	l64 = geo.add_line(p311, p312)
	
	# Superficies de la geometria
	# Superficies horizontales
	# Fondo
	b01 = geo.add_line_loop([l101, l31, -l201, -l11])
	s01 = geo.add_surface(b01)
	b02 = geo.add_line_loop([l102, l61, -l202, -l41])
	s02 = geo.add_surface(b02)
	# Geo-Hormigon
	b11 = geo.add_line_loop([l103, l22, -l203, -l12])
	s11 = geo.add_surface(b11)
	b12 = geo.add_line_loop([l104, l32, -l204, -l22])
	s12 = geo.add_surface(b12)
	b13 = geo.add_line_loop([l105, l42, -l205, -l32])
	s13 = geo.add_surface(b13)
	b14 = geo.add_line_loop([l106, l52, -l206, -l42])
	s14 = geo.add_surface(b14)
	b15 = geo.add_line_loop([l107, l62, -l207, -l52])
	s15 = geo.add_surface(b15)
	#Geo-Asfalto
	b21 = geo.add_line_loop([l108, l23, -l208, -l13])
	s21 = geo.add_surface(b21)
	b22 = geo.add_line_loop([l109, l53, -l209, -l23])
	s22 = geo.add_surface(b22)
	b23 = geo.add_line_loop([l110, l63, -l210, -l53])
	s23 = geo.add_surface(b23)
	#Asfalto
	b31 = geo.add_line_loop([l111, l24, -l211, -l14])
	s31 = geo.add_surface(b31)
	b32 = geo.add_line_loop([l112, l54, -l212, -l24])
	s32 = geo.add_surface(b32)
	b33 = geo.add_line_loop([l113, l64, -l213, -l54])
	s33 = geo.add_surface(b33)
	
	# Superficies Verticales
	# Longitudinales-Frente
	# Hormigon
	b101 = geo.add_line_loop([l101, lv131, -l104, -l103, -lv111])
	s101 = geo.add_plane_surface(b101)
	b102 = geo.add_line_loop([l102, lv161, -l107, -l106, -lv141])
	s102 = geo.add_plane_surface(b102)
	# Geo
	b111 = geo.add_line_loop([l103, lv122, -l108, -lv112])
	s111 = geo.add_surface(b111)
	b112 = geo.add_line_loop([l104, l105, l106, lv152, -l109, -lv122])
	s112 = geo.add_plane_surface(b112)
	b113 = geo.add_line_loop([l107, lv162, -l110, -lv152])
	s113 = geo.add_surface(b113)
	# Asfalto
	b121 = geo.add_line_loop([l108, lv123, -l111, -lv113])
	s121 = geo.add_surface(b121)
	b122 = geo.add_line_loop([l109, lv153, -l112, -lv123])
	s122 = geo.add_surface(b122)
	b123 = geo.add_line_loop([l110, lv163, -l113, -lv153])
	s123 = geo.add_surface(b123)
	# Longitudinales-Fondo
	# Hormigon
	b201 = geo.add_line_loop([l201, lv231, -l204, -l203, -lv211])
	s201 = geo.add_plane_surface(b201)
	b202 = geo.add_line_loop([l202, lv261, -l207, -l206, -lv241])
	s202 = geo.add_plane_surface(b202)
	# Geo
	b211 = geo.add_line_loop([l203, lv222, -l208, -lv212])
	s211 = geo.add_surface(b211)
	b212 = geo.add_line_loop([l204, l205, l206, lv252, -l209, -lv222])
	s212 = geo.add_plane_surface(b212)
	b213 = geo.add_line_loop([l207, lv262, -l210, -lv252])
	s213 = geo.add_surface(b213)
	# Asfalto
	b221 = geo.add_line_loop([l208, lv223, -l211, -lv213])
	s221 = geo.add_surface(b221)
	b222 = geo.add_line_loop([l209, lv253, -l212, -lv223])
	s222 = geo.add_surface(b222)
	b223 = geo.add_line_loop([l210, lv263, -l213, -lv253])
	s223 = geo.add_surface(b223)
	# Transversales
	b301 = geo.add_line_loop([l11, lv211, -l12, -lv111])
	s301 = geo.add_surface(b301)
	b302 = geo.add_line_loop([l12, lv212, -l13, -lv112])
	s302 = geo.add_surface(b302)
	b303 = geo.add_line_loop([l13, lv213, -l14, -lv113])
	s303 = geo.add_surface(b303)
	b311 = geo.add_line_loop([l31, lv231, -l32, -lv131])
	s311 = geo.add_surface(b311)
	b321 = geo.add_line_loop([l41, lv241, -l42, -lv141])
	s321 = geo.add_surface(b321)
	b331 = geo.add_line_loop([l61, lv261, -l62, -lv161])
	s331 = geo.add_surface(b331)
	b332 = geo.add_line_loop([l62, lv262, -l63, -lv162])
	s332 = geo.add_surface(b332)
	b333 = geo.add_line_loop([l63, lv263, -l64, -lv163])
	s333 = geo.add_surface(b333)
	
	# Volumenes de la Geometria
	# Base movil
	sbhm = geo.add_surface_loop([s01, s11, s12, s101, s201, s301, s311])
	vbhm = geo.add_volume(sbhm)
	# Base Fija
	sbhf = geo.add_surface_loop([s02, s14, s15, s102, s202, s321, s331])
	vbhf = geo.add_volume(sbhf)
	# Geo
	sg = geo.add_surface_loop(
		[s11, s12, s13, s14, s15, s21, s22, s23, s302, s332, s211, s212, s213, s111, s112, s113])
	vg = geo.add_volume(sg)
	# Asfalto
	sa = geo.add_surface_loop(
		[s31, s32, s33, s21, s22, s23, s303, s333, s221, s222, s223, s121, s122, s123])
	va = geo.add_volume(sa)
	
	# Superficies Fisicas
	geo.add_physical(s01, label=1)
	geo.add_physical(s02, label=2)
	# Volumenes Fisicos
	geo.add_physical(vbhm, label=3)
	geo.add_physical(vbhf, label=4)
	geo.add_physical(vg, label=5)
	geo.add_physical(va, label=6)
	
	# Escribimos el archivo .geo
	a = open('%s.geo' % (nombre_archivo), 'w')
	a.write(geo.get_code())
	a.close()
	# Realizar el mashado y el archivo xml
	s = 'gmsh -3  %s.geo -format msh2' % (nombre_archivo)
	os.system(s)
	s = 'dolfin-convert %s.msh %s.xml' % (nombre_archivo, nombre_archivo)
	os.system(s)
	if geometria:
		# Mostrar el resultado en gmsh
		s = 'gmsh %s.geo' % (nombre_archivo)
		os.system(s)
	if malla:
		# Mostrar el resultado en gmsh
		s = 'gmsh %s.msh' % (nombre_archivo)
		os.system(s)
	