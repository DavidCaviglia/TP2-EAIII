#TP2 EAIII, calculadora de microtiras, Grupo: Caviglia, Menendez, Velez  

import math as m
import cmath as cm
import numpy as np
#********************************************   DATOS   *************************************************************
#Permitividad relativa del dieléctrico (εr)
ep_rel = 4.4
#Altura o espesor del conductor (cobre) [mm]
t = 0.02
#Altura o espesor del sustrato (dieléctrico) [mm]
h = 1.5

#Ancho punta taper [mm]
taper_punta = 1.5

#Frecuencia
fo = 1.6e+9

#Impedancia de salida
zo = 50
#Impedancias para microtiras de desacople de polarización
zo_C = 25
zo_L = 75   

#Parámetros S
#BFP450  Vce=2 V ; Ic=70 mA ; f=1.6 GHz
s11 = cm.rect(0.7450 , m.radians(158.7))
s21 = cm.rect(4.34 , m.radians(62.1))
s12 = cm.rect(0.0616 , m.radians(54.3))
s22 = cm.rect(0.5332 , m.radians(168.3))

print ("\n", "Frecuencia:")
print ("fo =", format(fo,'.1E'), "[Hz]")
print ("\n", "Impedancia de salida:")
print ("Zo =", zo, "[Ω]")
print ("\n", "Impedancias para microtiras de desacople de polarización:")
print ("Zo_C =", zo_C, "[Ω]")
print ("Zo_L =", zo_L, "[Ω]")
print ("\n", "Parámetros S:")
print ("S11 =", round(abs(s11),6), "|__", round(np.angle(s11, deg=True),2), "°") 
print ("S12 =", round(abs(s12),6), "|__", round(np.angle(s12, deg=True),2), "°")
print ("S21 =", round(abs(s21),6), "|__", round(np.angle(s21, deg=True),2), "°")
print ("S22 =", round(abs(s22),6), "|__", round(np.angle(s22, deg=True),2), "°")
print ("\n", "Permitividad relativa del dieléctrico:")
print("εr =", ep_rel)
print ("\n", "Altura o espesor del conductor (cobre):")
print("t =", t, "[mm]")
print ("\n", "Altura o espesor del sustrato (dieléctrico):")
print("h =", h, "[mm]")
print ("\n", "Ancho de punta de los taper:")
print("W2 =", taper_punta, "[mm]")
#***********************************************************************************************************************



#Verificación de estabilidad incondicional:
deltaRollet = (s11*s22)-(s12*s21)
factorRollet = (1 - abs(s11)**2 - abs(s22)**2 + abs(deltaRollet)**2) / (2*abs(s12*s21))

#Diseño para máxima ganancia
b1 = 1 + abs(s11)**2 - abs(s22)**2 - abs(deltaRollet)**2
b2 = 1 + abs(s22)**2 - abs(s11)**2 - abs(deltaRollet)**2
c1 = s11 - (deltaRollet * s22.conjugate())
c2 = s22 - (deltaRollet * s11.conjugate())

#Coeficientes de reflexión
rms = (b1 - np.sqrt(b1**2 - 4*(abs(c1)**2)))/(2*(abs(c1)**2)) * c1.conjugate()
rml = (b2 - np.sqrt(b2**2 - 4*(abs(c2)**2)))/(2*(abs(c2)**2)) * c2.conjugate()

#Impedancias serie
zin_s = (zo * ((1+rms)/(1-rms))).conjugate()
zout_s = (zo * ((1+rml)/(1-rml))).conjugate()

#Adaptación a la entrada mediante pasaje a modelo paralelo
Rin_p = zin_s.real * (1 + (zin_s.imag/zin_s.real)**2)
Xin_p = zin_s.imag * (1 + (zin_s.real/zin_s.imag)**2)
#Para adaptar la parte imaginaria
Cin = 1 / (2 * m.pi * fo * Xin_p)
#Para adaptar la parte real
zo_in = np.sqrt(Rin_p*zo)

#Adaptación a la salida
#Para adaptar la parte real
zo_out = np.sqrt(zout_s.real*zo)
#Para adaptar la parte imaginaria
xlp = zo_out**2 / -zout_s.imag
Cout = 1 / (2 * m.pi * fo * xlp)

#Ganancias
#Máxima estable
gmsg = 10*m.log10(abs(s21)/abs(s12))
#De potencia de transconductancia máxima
gt_max = 10*m.log10((abs(s21)/abs(s12)) * (factorRollet - np.sqrt((factorRollet**2)-1)))


print ("\n", "Verificación de estabilidad incondicional:")
print ("|Δ| =", round(abs(deltaRollet),3))
print ("K =", round(factorRollet,3))

print ("\n", "Diseño para máxima ganancia:")
print ("B1 =", round(b1,3))
print ("B2 =", round(b2,3))
print ("C1 =", round(abs(c1),3), "|__", round(np.angle(c1, deg=True),2), "°")
print ("C2 =", round(abs(c2),3), "|__", round(np.angle(c2, deg=True),2), "°")

print ("\n", "Coeficientes de reflexión:")
print ("ΓmS =", round(abs(rms),3), "|__", round(np.angle(rms, deg=True),2), "°")
print ("ΓmL =", round(abs(rml),3), "|__", round(np.angle(rml, deg=True),2), "°")

print ("\n", "Impedancias serie:")
print ("Zin (serie) =", np.round(zin_s,3), "[Ω]")
print ("Zout (serie) =", np.round(zout_s,3), "[Ω]")

print ("\n", "Adaptación a la entrada mediante pasaje a modelo paralelo:")
print ("Rin (paralelo) =", np.round(Rin_p,3), "[Ω]")
print ("Xin (paralelo) =", np.round(Xin_p,3), "[Ω]")
print ("Para adaptar la parte imaginaria:")
print ("Cin =", format(Cin,'.1E'), "[F]")
print ("Para adaptar la parte real:")
print ("Zoin =", np.round(zo_in,3), "[Ω]")

print ("\n", "Adaptación a la salida:")
print ("Para adaptar la parte real:")
print ("Zo out =", np.round(zo_out,3), "[Ω]")
print ("Para adaptar la parte imaginaria:")
# print ("XLp =", np.round(xlp,3))
print ("Cout =", format(Cout,'.1E'), "[F]")

print ("\n", "Ganancia máxima estable:")
print ("Gmsg =", round(gmsg,3),"[dB]")
print ("Ganancia máxima de potencia de transconductancia:")
print ("Gt máx =", round(gt_max,3),"[dB]")

#********************************************   MICROTIRAS   *************************************************************
def microtira_WsobreH (Zo):
    #Ecuación 7 Hammerstad
    A = (Zo/60)*np.sqrt((ep_rel+1)/2) + ((ep_rel-1)/(ep_rel+1))*(0.23 + (0.11/ep_rel))
    #Ecuación 8 Hammerstad
    B = (377 * m.pi) / (2 * Zo * np.sqrt(ep_rel))
    #Ecuación 5 Hammerstad
    WsobreH_ec5 = (8 * m.e**A) / (m.e**(2*A) - 2)
    #Ecuación 6 Hammerstad
    WsobreH_ec6 = (2/m.pi) * (B - 1 - m.log(2*B - 1) + (ep_rel-1)/(2*ep_rel)*(m.log(B-1)+0.39-(0.61/ep_rel)))
    #Descarte de la ecuación 5 o 6
    WsobreH = 0
    if (WsobreH_ec5 <= 2):
        WsobreH = WsobreH_ec5
    if (WsobreH_ec6 >= 2):
        WsobreH = WsobreH_ec6
    return WsobreH

def microtira_we (Zo):
    WsobreH = microtira_WsobreH(Zo)
    #Ancho de la microtira
    w =  WsobreH * h
    #Ancho efectivo de la microtira
    if (WsobreH <= (1/(2*m.pi))):
        we = w + (t/m.pi)*(1+m.log((4*m.pi*w)/t)) #Ec. 9 Hammerstad
    if (WsobreH >= (1/(2*m.pi))):
        we = w + (t/m.pi)*(1+m.log((2*h)/t)) #Ec. 10 Hammerstad
    return we

def microtira_d (Zo,Zd):
    WsobreH = microtira_WsobreH(Zo)
    #Ancho de la microtira
    w =  WsobreH * h
    #Permetividad relativa efectiva del dielectrico (εr)
    if (WsobreH <= 1):
        ep_rel_ef = (ep_rel+1)/2 + ((ep_rel-1)/2)*(1/(np.sqrt(1+((12*h)/w))) + 0.04*(1-WsobreH**2)) #Ec. 2 Hammerstad
    if (WsobreH >= 1):
         ep_rel_ef = (ep_rel+1)/2 + ((ep_rel-1)/2)*(1/(np.sqrt(1+((12*h)/w)))) #Ec. 4 Hammerstad
    #Longitud de onda (λo)
    lambda_0 =  300/(fo/1e6) 
    #Longitud de onda efectiva (λe)
    lambda_ef =  lambda_0/np.sqrt(ep_rel_ef)
    #Constante de fase (β)
    beta = (2*m.pi)/lambda_ef
    d = (1/beta)*np.arctan(Zo/Zd)
    return d*1000

#Largo de microtira para los transformadores λ/4
def microtira_d_tr(Zo):
    WsobreH = microtira_WsobreH(Zo)
    #Ancho de la microtira
    w =  WsobreH * h
    #Permetividad relativa efectiva del dielectrico (εr)
    if (WsobreH <= 1):
        ep_rel_ef = (ep_rel+1)/2 + ((ep_rel-1)/2)*(1/(np.sqrt(1+((12*h)/w))) + 0.04*(1-WsobreH**2)) #Ec. 2 Hammerstad
    if (WsobreH >= 1):
        ep_rel_ef = (ep_rel+1)/2 + ((ep_rel-1)/2)*(1/(np.sqrt(1+((12*h)/w)))) #Ec. 4 Hammerstad
    #Longitud de onda (λo)
    lambda_0 =  300/(fo/1e6) 
    #Longitud de onda efectiva (λe)
    lambda_ef =  lambda_0/np.sqrt(ep_rel_ef)
    d_tr = lambda_ef/4
    return d_tr*1000

# Cálculo de taper
#      *  ****  *
#     **  *  *  **
#    * *  *  *  * *
#   ****  ****  **** 
def taper(w_microtira, l_microtira):
    base_triangulo = (w_microtira-taper_punta)/2
    altura_taper = base_triangulo * (m.tan((60*m.pi)/180))
    area_taper = (base_triangulo*altura_taper) + (taper_punta*altura_taper)
    area_microtira = (w_microtira * l_microtira)-area_taper
    l_microtira = area_microtira / w_microtira
    return (altura_taper, l_microtira)


print ("\n", "Diseño con microtiras:")
print ("Transformador λ/4 entrada:")
print ("We =", round(microtira_we(zo_in),3), "[mm]")
print ("d =", round(microtira_d_tr(zo_in),3), "[mm]")
print ("Transformador λ/4 salida:")
print ("We =", round(microtira_we(zo_out),3), "[mm]")
print ("d =", round(microtira_d_tr(zo_out),3), "[mm]")
print ("Síntesis capacitor adaptador entrada:")
w_c_in = microtira_we(zo)
l_c_in = microtira_d(zo,Xin_p)
taper_c_in = taper(w_c_in, l_c_in)
print ("W microtira =", round(w_c_in,3), "[mm]")
print ("L microtira =", round(taper_c_in[1],3), "[mm]")
print ("L taper =", round(taper_c_in[0],3), "[mm]")
print ("Síntesis capacitor adaptador salida:")
w_c_out = microtira_we(zo)
l_c_out = microtira_d(zo,xlp)
taper_c_out = taper(w_c_out, l_c_out)
print ("W microtira =", round(w_c_out,3), "[mm]")
print ("L microtira =", round(taper_c_out[1],3), "[mm]")
print ("L taper =", round(taper_c_out[0],3), "[mm]")
print ("Transformadores λ/4 desacopladores de polarización:")
print ("We =", round(microtira_we(zo_L),3), "[mm]")
print ("d =", round(microtira_d_tr(zo_L),3), "[mm]")
print ("Síntesis de capacitores desacopladores de polarización:")
w_c_pol = microtira_we(zo)
l_c_pol = microtira_d(zo,zo_C)
taper_c_pol = taper(w_c_pol, l_c_pol)
print ("W microtira =", round(w_c_pol,3), "[mm]")
print ("L microtira =", round(taper_c_pol[1],3), "[mm]")
print ("L taper =", round(taper_c_pol[0],3), "[mm]")
#*************************************************************************************************************************