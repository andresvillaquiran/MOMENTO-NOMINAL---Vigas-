# librerias
import tkinter as tk 
from math import sqrt
from numpy import array,around
from math import sqrt
from sympy import symbols, solve, Eq
from scipy.optimize import fsolve, minimize
import matplotlib.pyplot as plt

phi_simbol = '\u03C6'

# creamos ventana principal

def mostrar_frame(frame):
    frame.tkraise() ## trae al frame que sea llamado

'''
def acero_colocado(combi):
    as_colocado = 0

    combinacioanes = combi.split('+') # Separamos por el simbolo + para obtener cada combinacion de acero

    for combinacion in combinacioanes:
        combinacion = combinacion.strip() # Eliminamos espacios en blanco
        if '#' in combinacion:
            num_barras , calibre = combinacion.split('#')
            as_colocado += float(num_barras)*lista_refuerzo[int(calibre),2]/1000000
        else:
            as_colocado += eval(combinacion)/10000  # Si no hay #, evaluamos la expresión directamente (cambiamos a m2 para realizar comparativa)
    return as_colocado
'''
## LEEMOS EL ACERO COLOCADO
def acero_colocado(combi):
    as_colocado = 0
    calis = [] # lista de los calibres en las combinaciones
    combinaciones = combi.split('+') # Separamos por el simbolo + para obtener cada combinacion de acero

    for combinacion in combinaciones:
        combinacion = combinacion.strip() # Eliminamos espacios en blanco
        if '#' in combinacion:
            num_barras , calibre = combinacion.split('#')
            as_colocado += float(num_barras)*lista_refuerzo[int(calibre),2]/1000000
            calis.append(int(calibre))
        else:
            as_colocado += eval(combinacion)/10000  # Si no hay #, evaluamos la expresión directamente (cambiamos a m2 para realizar comparativa)
    if '#' in combinacion:
        db = max(calis)/1000 # definimos el db [m] como el diametro mayor de las combinaciones.
    else:
        db = 0


    return as_colocado,db

# DEFINICION DE LA FUNCION DE AREA QUE REPRESENTARA LA SECCION

def area_compresion (c,B1,b):
    Ac = B1*c*b
    return Ac

## DEFINIMOS EL EJE NEUTRO DE LA SECCION  
def eje_neutro(dp,euc,d,Es,fy,fc,B1,Asp,As,tipo,sup,b): 

    def Eq1(c):
        
        cv =c[0]
        if tipo == "sub":
            fs = fy
        else:
            fs = Es*euc*(d-cv)/cv
        if sup == "fluye":
            fsp = fy
        elif sup == "no_fluye":
            fsp = Es*euc*(cv-dp)/cv

        # Planteamos ecuación
        return (0.85*fc*(area_compresion(cv,B1,b)-Asp) + Asp*fsp - As*fs)**2
    # Solucion del sistema de ecuaciones
    restricciones = [ 
                {'type': 'ineq', 'fun': lambda c: c[0] + d },
                {'type': 'ineq', 'fun': lambda c: c[0] - 0.0 },
                {'type': 'ineq', 'fun': lambda c: c[0] - 0.005 },
    ]
    resultado = minimize(Eq1, x0=[ 0.06], constraints=restricciones)
    c = resultado.x

    return c



def Analisis_viga(fc,fy,b,h,r,dint,dpint,ref_tracc, ref_comp):
    Es = 200000 #[MPa]
    euc = 0.003
    ey = fy/Es

    As,dbt = acero_colocado(ref_tracc)
    Asp,dbc = acero_colocado(ref_comp)

    # Altura efectiva a compresion y traccion

    if dint == '0' or dint == 'N/A':
        d = h - r - dbt/2  # Calculo altura efectiva de la viga
    else:
        d = eval(dint)

    if dpint == '0' or dpint == 'N/A':
        dp = r + dbc/2  # Calculo altura efectiva de la viga
    else:
        dp = eval(dpint)

    # VERIFICO ACERO MINIMO
    As_m1 =0.25*sqrt(fc)*b*d/fy 
    As_m2 = 1.4*b*d/fy
    As_min = max(As_m1, As_m2) # Acero minimo requerido

    # ¿Cumple acero minimo?
    if As < As_min:
        print("NO cumple acero minimo")
    else: 
        print("CUMPLE acero minimo requerido")

    # Calculamos B1 segun f'c
    if fc >= 17:
        if fc < 28:
            B1 = 0.85
        elif fc < 56:
            B1 = 0.85-0.05*(fc-28)/7
        else:
            B1 = 0.65
    else: 
        print("Defina un fc valido")        

    #CONDICION BALANCEADA
    cb = euc*d/(euc+ey)  # Eje neutro balanceado
    espb = euc*(cb-dp)/cb # def acero comp en balanceada
    if espb >= ey: # defino esfuerzo de As' en balanceada
        fspb = fy
    else:
        fspb = Es*espb

    Acb = area_compresion(cb, B1,b)
    Asb = (0.85*fc*(Acb-Asp)+Asp*fspb)/fy  # Area acero balanceado
    Asmax = 0.85*fc*B1*b*(d*euc/(euc+(ey+euc)))/fy # Acero maximo para es=0.0051

    # Definimos el tipo de viga segun el acero de refuerzo
    if As < Asb:
        print("VIGA SUB_REFORZADA")
        tipo = "sub"

    else:
        print("VIGA SOBRE_REFORZADA")
        tipo = "sobre"    

    ############################################################################
    ############################################################################
    ## PROCESO ITERATIVO PARA DETERMINAR SU ESTADO DE ESFUERZO

    #SUPONEMOS QUE As' fluye y hacemos el calculo
    # As' fluye
    sup = "fluye"
    #-----------------------------------------------------------------
    # Solucion del sistema de ecuaciones
    resultado = eje_neutro(dp,euc,d,Es,fy,fc,B1,Asp, As, tipo,sup,b)
    c = float(resultado[0])
    #  ---------------------------------------------------------------
    # Verificamos deformación segun suposición
    esp = euc*(c-dp)/c

    if esp >= ey :
        print("FLUYE acero compresion")
        fsp = fy
    else:
        sup = "no_fluye"  
        # ------------------------------------------------------------
        # Solucion del sistema de ecuaciones 
        resultado = eje_neutro(dp,euc,d,Es,fy,fc,B1,Asp, As, tipo,sup,b)
        c = float(resultado[0])
        #-------------------------------------------------------------
        esp = euc*(c-dp)/c
        fsp = Es*esp
        print("NO FLUYE acero compresion ")

    # Calculamos altura de el cubo de whitney
    a = B1*c

    es = euc*(d-c)/c
   

    if es > (ey+euc):
        condicion = "traccion"
        print("Viga controlada por ", condicion)
    elif es < (ey+euc) and es > ey:
        condicion= "transicion"
        print("Viga controlada por ", condicion)
    else:
        condicion = "compresion"
        print("Viga controlada por ", condicion)

    if As < Asb:
        tipo = "sub"
        fs = fy
    else:
        tipo = "sobre"
        fs = Es*euc*(d-c)/c
    
    Ac = area_compresion(c,B1,b)

    # Calculamos el curvatura 
    curvatura = euc/c
    print("curvatura = ", curvatura)  
    print(f'd = {d}')
    print(f'c = {c}')  
    print(f'fs = {fs}')
    print(f'es = {es}')
    print(f'As = {As}') 
    print(f'esp = {esp}') 

    # CALCULAMOS EL ESFUERZO EN EL ACERO Y EL MOMENTO NOMINAL RESISTENTE DE LA VIGA
    # PARA CADA TIPO DE VIGA Y DE ESTRIBO DEFINIMOS EL FACTOR DE REDUCCION DE RESISTENCIA PHI
    # VIGA SOBREFORZADA

    if tipo == "sobre":

        phi = 0.65    

    ## VIGA SUB REFORZADA   
    elif tipo == "sub":
        fs = fy 
        if condicion == "traccion":
            phi = 0.9
        elif condicion == "transicion":

            phi = 0.65 + 0.25*(es-ey)/euc

    ycentroidal = a/2
    Mn = 0.85*fc*1000*(area_compresion(c,B1,b)-Asp)*(d-ycentroidal) + Asp*fsp*1000*(d-dp) #[kN-m]
    phi_Mn = phi*Mn


    return phi,phi_Mn, As_min,Asb ,Asmax, es, esp, fs,fsp, d, dp, condicion, curvatura,c,As,Asp



def diagrama(es,d,dp,esp,fy,c,phi_Mn):
    euc = 0.003
    ey = fy/200000
    dyp = c - c*ey/euc


    ## CREACION DEL DIAGRAMA DE CURVATURA DEFORMACION DE LA SECCION
    # IMPRESION DEL DIAGRAMA CURVATURA DEFORMACION
    coord = array([
        [  0 ,  0 ],
        [ euc , 0 ],
        [ -es , -d],
        [  0 ,  -d],
        [   0 , 0 ],
    ])

    coord2 = array([
        [0, -dp],
        [float(esp), -dp],
    ])  
    fluy_comp = array([ 
                [0 , - dyp],
                [ey, -dyp],
    ])

    fig , ax = plt.subplots()
    ax.plot(coord[:,0], coord[:,1], marker = "o") # Dibuja los puntos
    ax.plot(coord2[:,0], coord2[:,1], marker = "o", color = "green", label = "d'")
    ax.plot(fluy_comp[:,0], fluy_comp[:,1], marker = "o", color = "red", label = "dy'")
    ax.legend()
    plt.title(label = f" DIAGRAMA CURVATURA - DEFORMACIÓN \u03C6Mn = {round(phi_Mn, 2)} [kN-m]")
    plt.xlabel("Deformación [m/m]") # Titulo eje x
    plt.ylabel("altura efectiva [m]") # Titulo eje y

    plt.grid(True) # Dibuje la malla

    plt.show()
    #ax.set_aspect("equal") # Ejes misma escala
    
    return


#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Base datos refuerzo col1 = tipo de acero, col 2 = diametro [mm], col 3 = area[m]
lista_refuerzo = array([[0,0,0],
                        [1,0,0],
                        [2,6.4,32],
                        [3,9.5,71],
                        [4,12.7,129],
                        [5,15.9,199],
                        [6,19.1,284],
                        [7,22.2,387],
                        [8,25.4,510],
                        [9,28.7,645],
                        [10,32.3,819],
                        [11,35.8,1006],
                        [12,38.1,1134],
                        [13,41.3,1335],
                        [14,43,1452],
                        [18,57.3,2581]])


ventana1 = tk.Tk()
ventana1.geometry("1200x600") ## tamaño de la ventana general
ventana1.title("Calculo de momento nominal resistente") ## titulo de la ventana general

# Creamos los frames para organizar la información

frame1 = tk.Frame(ventana1)
frame2 = tk.Frame(ventana1)

for frame in (frame1,frame2):  ## posicionamos los frames en el mismo lugar
    frame.grid(row=0, column=0, sticky='nsew')

##################################### CONTENIDO FRAME 1 #################################
#####################################                   #################################

etiq_inicio = tk.Label(frame1, text='Bienvenido a la sufridera DCR1 -- Diseñado por Andres Villaquirán ')
etiq_inicio.grid(row=0, column=0, columnspan=3, padx=5, pady=3)


etiq_ing_valores = tk.Label(frame1, text='Ingreso de parametro de la viga')
etiq_raviso = tk.Label(frame1, text='Si hay estribos sumar diametro a el recubrimiento')
etiq_ing_valores.grid(row=1, column=0, padx=5, pady=3)
etiq_raviso.grid(row=1, column=1, columnspan=2 , padx=5, pady=3) # Ubicación en la cuadrícula

etiq_fc = tk.Label(frame1, text='fc [MPa]')
entrada_fc = tk.Entry(frame1)
etiq_fc.grid(row=2, column=0, padx=5, pady=3) # Ubicación en la cuadrícula
entrada_fc.grid(row=2, column=2, padx=3, pady=3) # Ubicación en la cuadrícula

etiq_fy = tk.Label(frame1, text='fy [MPa]')
entrada_fy = tk.Entry(frame1)
etiq_fy.grid(row=3, column=0, padx=5, pady=3) # Ubicación en la cuadrícula
entrada_fy.grid(row=3, column=2, padx=3, pady=3) # Ubicación en la cuadrícula

etiq_b = tk.Label(frame1, text='b [m]')
entrada_b  = tk.Entry(frame1)
etiq_b.grid(row=4, column=0, padx=10, pady=3) # Ubicación en la cuadrícula
entrada_b.grid(row=4, column=2, padx=3, pady=3) # Ubicación en la cuadrícula

etiq_h = tk.Label(frame1, text='h [m]')
entrada_h  = tk.Entry(frame1)
etiq_h.grid(row=5, column=0, padx=5, pady=3) # Ubicación en la cuadrícula
entrada_h.grid(row=5, column=2, padx=3, pady=3) # Ubicación en la cuadrícula

etiq_r = tk.Label(frame1, text='r [m]')
entrada_r  = tk.Entry(frame1)
etiq_r.grid(row=6, column=0, padx=10, pady=3) # Ubicación en la cuadrícula
entrada_r.grid(row=6, column=2, padx=3, pady=3) # Ubicación en la cuadrícula

etiq_dint = tk.Label(frame1, text='d (opcional)[m]')
entrada_dint  = tk.Entry(frame1)
etiq_dint.grid(row=7, column=0, padx=10, pady=3) # Ubicación en la cuadrícula
entrada_dint.grid(row=7, column=2, padx=3, pady=3) # Ubicación en la cuadrícula

etiq_dpint = tk.Label(frame1, text='dp (opcional) [m]')
entrada_dpint  = tk.Entry(frame1)
etiq_dpint.grid(row=8, column=0, padx=10, pady=3) # Ubicación en la cuadrícula
entrada_dpint.grid(row=8, column=2, padx=3, pady=3) # Ubicación en la cuadrícula

etiq_ref_tracc = tk.Label(frame1, text='Refuerzo tracción [cm2]')
entrada_ref_tracc  = tk.Entry(frame1)
etiq_ref_tracc.grid(row=9, column=0, padx=10, pady=3) # Ubicación en la cuadrícula
entrada_ref_tracc.grid(row=9, column=2, padx=3, pady=3) # Ubicación en la cuadrícula

etiq_ref_comp = tk.Label(frame1, text='Refuerzo compresión [cm2]')
entrada_ref_comp  = tk.Entry(frame1)
etiq_ref_comp.grid(row=10, column=0, padx=10, pady=3) # Ubicación en la cuadrícula
entrada_ref_comp.grid(row=10, column=2, padx=3, pady=3) # Ubicación en la cuadrícula


## ventanas de resultados
tk.Label(frame1,text =' RESULTADOS ANALISIS').grid(row=1, column=3, columnspan=3, padx=1, pady=1) # Ubicación en la cuadrícula


etiq_As_min = tk.Label(frame1,text='As min [cm2]')
etiq_As_min.grid(row=2, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_Asb = tk.Label(frame1,text='As balanceado [cm2]')
etiq_Asb.grid(row=3, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_Amax = tk.Label(frame1,text ='As max [cm2]')         
etiq_Amax.grid(row=4, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_As = tk.Label(frame1,text ='As colocado [cm2]')         
etiq_As.grid(row=5, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_Asp = tk.Label(frame1,text ='Asp colocado [cm2]')         
etiq_Asp.grid(row=6, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_d = tk.Label(frame1,text ='d [m]')
etiq_d.grid(row=7, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_dp = tk.Label(frame1,text ='dp [m]')
etiq_dp.grid(row=8, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_c = tk.Label(frame1,text ='c [m]')
etiq_c.grid(row=9, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_es = tk.Label(frame1,text ='\u03B5s [m/m]')
etiq_es.grid(row=10, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_fs = tk.Label(frame1,text ='Fs [MPa]')
etiq_fs.grid(row=11, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_esp = tk.Label(frame1,text ='\u03B5sp [m/m]')
etiq_esp.grid(row=12, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_fsp = tk.Label(frame1,text ='Fsp [MPa]')
etiq_fsp.grid(row=13, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_curv = tk.Label(frame1,text ='Curvatura \u03A6 [1/m]')
etiq_curv.grid(row=14, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_tipo = tk.Label(frame1,text='Controlada por:')      
etiq_tipo.grid(row=15, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_phi = tk.Label(frame1,text = f'{phi_simbol} factor')
etiq_phi.grid(row=16, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

etiq_phi_Mn = tk.Label(frame1,text = f'{phi_simbol}Mn [kN-m]')
etiq_phi_Mn.grid(row=17, column=4, padx=5, pady=3) # Ubicación en la cuadrícula

## Etiqueta de resultados

result_As_min = tk.Label(frame1)
result_As_min.grid(row=2, column=5, padx=5, pady=3) # Ubicación en la cuadrí

result_Asb = tk.Label(frame1)
result_Asb.grid(row=3, column=5, padx=5, pady=3) # Ubicación en la cuadríc

result_Amax = tk.Label(frame1)         
result_Amax.grid(row=4, column=5, padx=5, pady=3) # Ubicación en la cuadrí

result_As = tk.Label(frame1)         
result_As.grid(row=5, column=5, padx=5, pady=3) # Ubicación en la cuadrí

result_Asp = tk.Label(frame1)         
result_Asp.grid(row=6, column=5, padx=5, pady=3) # Ubicación en la cuadrí

result_d = tk.Label(frame1)
result_d.grid(row=7, column=5, padx=5, pady=3) # Ubicación en la cuadríc

result_dp = tk.Label(frame1)
result_dp.grid(row=8, column=5, padx=5, pady=3) # Ubicación en la cuadrícula

result_c = tk.Label(frame1)
result_c.grid(row=9, column=5, padx=5, pady=3) # Ubicación en la cuadrícula

result_es = tk.Label(frame1)
result_es.grid(row=10, column=5, padx=5, pady=3) # Ubicación en la cuadríc

result_fs = tk.Label(frame1)
result_fs.grid(row=11, column=5, padx=5, pady=3) # Ubicación en la cuadríc

result_esp = tk.Label(frame1)
result_esp.grid(row=12, column=5, padx=5, pady=3) # Ubicación en la cuadríc

result_fsp = tk.Label(frame1)
result_fsp.grid(row=13, column=5, padx=5, pady=3) # Ubicación en la cuadríc

result_curv = tk.Label(frame1)
result_curv.grid(row=14, column=5, padx=5, pady=3) # Ubicación en la cuadríc

result_tipo = tk.Label(frame1)      
result_tipo.grid(row=15, column=5, padx=5, pady=3) # Ubicación en la cuadrícula

result_phi = tk.Label(frame1)
result_phi.grid(row=16, column=5, padx=5, pady=3) # Ubicación en la cuadríc

result_phi_Mn = tk.Label(frame1)
result_phi_Mn.grid(row=17, column=5, padx=5, pady=3) # Ubicación en la cuadríc


def calculo_entradas():   # Leemos las variables que entraron al principio

    fc = float(entrada_fc.get())
    fy = float(entrada_fy.get()) 
    b = float(entrada_b.get())
    h = float(entrada_h.get())
    r = float(entrada_r.get())
    dint = entrada_dint.get()
    dpint = entrada_dpint.get()
    ref_tracc = entrada_ref_tracc.get() 
    ref_comp =  entrada_ref_comp.get() 
    
    
    #Analisis_viga(fc,fy,b,h,r,dint,dpint,ref_tracc, ref_comp)
    phi,phi_Mn, As_min,Asb ,Asmax, es, esp, fs,fsp, d, dp, condicion, curvatura, c,As,Asp = Analisis_viga(fc,fy,b,h,r,dint,dpint,ref_tracc, ref_comp)

    result_As_min['text'] = round(As_min*10000,4)
    result_Asb['text'] = round(Asb*10000,4)
    result_Amax['text'] = round(Asmax*10000,4)
    result_As['text'] = round(As*10000,4)
    result_Asp['text'] = round(Asp*10000,4)
    result_d['text'] = round(d,5)
    result_dp['text'] = round(dp,5)
    result_c['text'] = round(c,6)
    result_es['text']= round(es,5)
    result_fs['text']= round(fs,2)
    result_esp['text']= round(esp,5)
    result_fsp['text']= round(fsp,2)
    result_curv['text']= round(curvatura,5)
    result_tipo['text']= condicion
    result_phi['text'] = round(phi,3)
    result_phi_Mn['text'] = round(phi_Mn,3)


    diagrama(es,d,dp,esp,fy,c,phi_Mn)

    global gAsmax , gAsmin, gAsb , gfc, gfy, gb, gh, gr
    gAsmax = Asmax
    gAsmin = As_min
    gAsb = Asb
    gfc = fc
    gfy = fy
    gb = b
    gh = h
    gr = r



boton_calcu1=tk.Button(frame1,text = 'Calcular', command=lambda: calculo_entradas()) # Botón para ejecutar la función de lectura de variables
boton_calcu1.grid(row=11, column=1, padx=5, pady=3) # Ubicación en la cuadrícula

boton_cambio_f2=tk.Button(frame1, text='Ayuda/información', command= lambda: mostrar_frame(frame2))
boton_cambio_f2.grid(row=19, column=1, padx=5, pady=3) # Ubicación en la cuadrícula del boton




boton_regresar=tk.Button(frame2, text='Regresar', command= lambda: mostrar_frame(frame1)) 
#boton_regresar.grid(row = pos_boton, column=4, padx=5, pady=3) # Ubicación en la cuadrícula del boton

mostrar_frame(frame1)

ventana1.mainloop()






