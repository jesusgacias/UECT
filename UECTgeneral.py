import numpy as np
from FuncionesUECT_V2 import HallarLugar, MiAngulo
import matplotlib.pyplot as plt
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from PIL import Image

#Definición del complejo P=(V,f,A,C) que representa a una pajarita

def Pajarita():
    """
    Genera el cpp de una pajarita.
    """
    V=range(12)
    A=[
        [0,1],[0,2],[0,3],
        [1,4],[1,5],[1,8],
        [2,4],[2,6],
        [3,5],[3,7],
        [4,6],[4,8],
        [5,7],[5,8],
        [6,10],
        [7,9],
        [8,9],[8,10],[8,11],
        [9,11],
        [10,11]
    ]
    C=[
        [0,1,4,2],[0,3,5,1],
        [1,5,8],[1,8,4],
        [2,4,6],
        [3,7,5],
        [4,8,10,6],
        [5,7,9,8],
        [8,9,11],
        [8,11,10]
    ]
    f=[
        np.array([0,0]),
        np.array([1,0]),
        np.array([1,1]),
        np.array([1,1]),
        np.array([2,1]),
        np.array([2,1]),
        np.array([2,0]),
        np.array([2,0]),
        np.array([1,2]),
        np.array([1,1]),
        np.array([1,1]),
        np.array([2,2])
    ]
    return [V,f,A,C]

def Mitad():
    """
    Genera el cpp de una hoja cuadrada plegada por la mitad.
    """
    V=range(6)
    A=[
        [0,1], [0,5],
        [1,2], [1,4],
        [2,3],
        [3,4],
        [4,5]
    ]
    C=[
        [0,1,4,5],
        [1,2,3,4]
    ]
    f=[
        np.array([1,2]),
        np.array([0,2]),
        np.array([1,2]),
        np.array([1,0]),
        np.array([0,0]),
        np.array([1,0])
    ]
    return [V,f,A,C]

def Cuatro():
    """
    Genera el cpp de una hoja cuadrada plegada en cuatro.
    """
    V=range(9)
    A=[
        [0,1],[0,3],[0,5],[0,7],
        [1,2],[1,8],
        [2,3],
        [3,4],
        [4,5],
        [5,6],
        [6,7],
        [7,8]
    ]
    C=[
        [0,1,8,7],
        [0,1,2,3],
        [0,3,4,5],
        [0,5,6,7]
    ]
    f=[
        np.array([0,0]),
        np.array([0,1]),
        np.array([1,1]),
        np.array([1,0]),
        np.array([1,1]),
        np.array([0,1]),
        np.array([1,1]),
        np.array([1,0]),
        np.array([1,1])
    ]
    return [V,f,A,C]

def PuntosPajarita():
    """
    Genera el cpp de tan solo los vértices de una pajarita.
    """
    V=range(12)
    A=[]
    C=[]
    f=[
        np.array([0,0]),
        np.array([1,0]),
        np.array([1,1]),
        np.array([1,1]),
        np.array([2,1]),
        np.array([2,1]),
        np.array([2,0]),
        np.array([2,0]),
        np.array([1,2]),
        np.array([1,1]),
        np.array([1,1]),
        np.array([2,2])
    ]
    return [V,f,A,C]

def Varita():
    """
    Genera el cpp de una varita rígida plegada por la mitad.
    """
    V=range(3)
    A=[
        [0,1],[1,2]
    ]
    C=[]
    f=[
        np.array([2,0]),
        np.array([1,0]),
        np.array([2,0])
    ]
    return [V,f,A,C]

def ImposibleHull(R): #R es el exceso de papel 
    """
    Genera el cpp del plegado imposible de Hull. Se recomienda utilizar HullGen sobre este.
    """
    V=range(12)
    A=[
        [0,1],[0,8],
        [1,2],[1,9],
        [2,3],[2,10],
        [3,4],
        [4,5],[4,10],
        [5,6],[5,11],
        [6,7],
        [7,8],[7,11],
        [8,9],
        [9,10],[9,11],
        [10,11]
    ]
    C=[
        [0,1,9,8],
        [1,2,10,9],
        [2,3,4,10],
        [4,5,11,10],
        [5,6,7,11],
        [8,9,11,7],
        [9,10,11]
    ]

    rot1=[[np.cos(2*np.pi/3), np.sin(2*np.pi/3)],[-np.sin(2*np.pi/3), np.cos(2*np.pi/3)]]
    rot2=[[np.cos(2*np.pi/3), -np.sin(2*np.pi/3)],[np.sin(2*np.pi/3), np.cos(2*np.pi/3)]]
    desp=np.array([0.5, np.sqrt(3)/6])
    f=[
        np.array([R*np.sqrt(3),R]),#0
        np.matmul(rot1,(np.array([1,R])-desp))+desp, #1
        np.matmul(rot1,(np.array([0,R])-desp))+desp, #2
        np.matmul(rot1,(np.array([R*np.sqrt(3),R])-desp))+desp, #3
        np.matmul(rot2,(np.array([1,R])-desp))+desp, #4
        np.matmul(rot2,(np.array([0,R])-desp))+desp, #5
        np.array([1-np.sqrt(3)*R,R]), #6
        np.array([1,R]), #7
        np.array([0,R]), #8
        np.array([0,0]), #9
        np.array([0.5,np.sqrt(3)/2]), #10
        np.array([1,0]), #11
    ]

    return(V,f,A,C)

def HullGen(R, theta, cen=False):
    """
    Devuelve el cpp de una variación del plegado de Hull con un ángulo theta variable (ver anexo C).
    R es la diferencia de grosor entre el triángulo equilátero que forman los pliegues y el de los bordes del papel.
    theta admite valores entre PI/6 y 5*PI/6, ambos excluidos. El plegado es posible para theta en los rangos (PI/6, PI/3] y [2*PI/3,5*PI/6), e imposible en el rango (PI/3, 2*PI/3).
    Poniendo cen=True, se centra el modelo en el origen, dando una UECT mucho más simétrica.
    """
    V=range(12)
    A=[
        [0,1],[0,8],
        [1,2],[1,9],
        [2,3],[2,10],
        [3,4],
        [4,5],[4,10],
        [5,6],[5,11],
        [6,7],
        [7,8],[7,11],
        [8,9],
        [9,10],[9,11],
        [10,11]
    ]
    C=[
        [0,1,9,8],
        [1,2,10,9],
        [2,3,4,10],
        [4,5,11,10],
        [5,6,7,11],
        [8,9,11,7],
        [9,10,11]
    ]

    rot1=[[np.cos(2*np.pi/3), np.sin(2*np.pi/3)],[-np.sin(2*np.pi/3), np.cos(2*np.pi/3)]]
    rot2=[[np.cos(2*np.pi/3), -np.sin(2*np.pi/3)],[np.sin(2*np.pi/3), np.cos(2*np.pi/3)]]
    desp=np.array([0.5, np.sqrt(3)/6])
    rottheta=[[np.cos(np.pi/2-theta), -np.sin(np.pi/2-theta)], [np.sin(np.pi/2-theta),np.cos(np.pi/2-theta)]]
    rotado=np.matmul(rottheta, np.array([-np.sqrt(3)*R,-R]))
    reflejado=[-rotado[0], rotado[1]]
    pos0=np.array([reflejado[0], -reflejado[1]])
    f=[
        pos0,#0
        np.matmul(rot1,(np.array([1-R/np.tan(theta),R])-desp))+desp, #1
        np.matmul(rot1,(np.array([-R/np.tan(theta),R])-desp))+desp, #2
        np.matmul(rot1,(pos0-desp))+desp, #3
        np.matmul(rot2,(np.array([1-R/np.tan(theta),R])-desp))+desp, #4
        np.matmul(rot2,(np.array([-R/np.tan(theta),R])-desp))+desp, #5
        np.matmul(rot2,(pos0-desp))+desp, #6
        np.array([1-R/np.tan(theta),R]), #7
        np.array([-R/np.tan(theta),R]), #8
        np.array([0,0]), #9
        np.array([0.5,np.sqrt(3)/2]), #10
        np.array([1,0]), #11
    ]
    if cen==True:
        for i in range(len(f)):
            f[i]=f[i]-desp

    return(V,f,A,C)


def Euler(V,A,C):
    """
    Calcula la característica de Euler de un cpp.
    """
    return len(V)-len(A)+len(C)


#Ordenado por alturas
def Orden(V,f,v):
    """
    Dado un conjunto de vértices V, una función geométrica f y una dirección v, ordena los vértices según el orden de choque en la dirección v.
    """
    alts=[]
    for i in range(len(V)):
        alts.append(np.matmul(v,f[i]))
    Listasunidas=sorted(zip(alts,V))
    alts_ord=[x[0] for x in Listasunidas]
    V_ord=[x[1] for x in Listasunidas]

    return([alts_ord,V_ord])

def AlturasIguales(O):
    """
    Dado un orden de alturas, detecta qué alturas son iguales entre sí y cuántos vértices hay coincidiendo en cada.
    """
    Ig=[0,1]
    for i in range(len(O[0])-1):
        if O[0][i]==O[0][i+1]:
            Ig[len(Ig)-1]=Ig[len(Ig)-1]+1
        else:
            Ig.append(Ig[len(Ig)-1]+1)
    return Ig

def CambiosEuler(V,f,A,C,v):
    """
    Dado un cpp (V,f,A,C) y una dirección v, devuelve las alturas en las que se producen cambios de la característica de Euler y cuáles son dichos cambios.
    """
    O=Orden(V,f,v)
    Alts=AlturasIguales(O)
    Pasos=len(Alts)-1
    Criticas=sorted(list(set(O[0]))) #Para quitar las alturas repetidas
    Eu=[0]
    V_M=[]

    for j in range(Pasos):
        V_M.extend(O[1][Alts[j]:Alts[j+1]])
        #print(V_M)
        A_M=[]
        C_M=[]
        for a in A:
            if a[0] in V_M and a[1] in V_M:
                A_M.append(a)
        #print(A_M)

        for c in C:
            C_M.append(c)
            N=len(c)
            for j in range(N):
                if c[j] not in V_M:
                    C_M.pop()
                    break
        #print(C_M)
        
        Eu.append(Euler(V_M,A_M,C_M))
        #print(Euler(V_M,A_M,C_M))

    

    return [Criticas, Eu]


def ValorEuler(V,f,A,C,v,t):
    """
    Devuelve el valor de la UECT del cpp (V,f,A,C) en dirección v y a altura t.
    """
    Cam=CambiosEuler(V,f,A,C,v)
    j=HallarLugar(t,Cam[0])

    return Cam[1][j]


def CurvaEuler(V,f,A,C,v, Paso=300, Dibujar=True, R=4):
    """
    Devuelve la curva de Euler de la UECT del cpp (V,f,A,C) en dirección v y con el Paso aportado, en el rango [-R,R] aportado.
    """
    X=np.arange(-R,R,2*R/Paso)
    D=[0 for _ in range(Paso)]
    i=0
    
    for x in X:
        D[i]=ValorEuler(V,f,A,C,v,x)
        i=i+1

    if Dibujar==True:
        fig, ax1=plt.subplots()
        ax1.plot(X,D)
        ax1.set_title('UECT en dirección')
        plt.gca().set_ylim(bottom=-2, top=4)
        show()
    return D


def CurvaSECT(V,f,A,C,v, Paso=300, Dibujar=True, R=4):
    """
    Devuelve la curva de Euler de la SUECT del cpp (V,f,A,C) en dirección v y con el Paso aportado, en el rango [-R,R] aportado.
    (LA INTEGRACION ES ALGO INEXACTA, CORREGIR SI SE QUIEREN HACER DISTANCIAS)
    """
    ECT=CurvaEuler(V=V, f=f, A=A, C=C, v=v, Paso=Paso, Dibujar=False, R=R)


    A=len(ECT)
    Media=np.mean(ECT)
    ECT=ECT-Media*np.array([1 for _ in range(A)])
    SECT=[0 for _ in range(A-1)]
    for i in range(A-1):
        SECT[i]=sum(ECT[0:i+1])/Paso
    SECT.append(0)

    if Dibujar==True:
        X=np.arange(-R,R,2*R/Paso)
        fig, ax1=plt.subplots()
        ax1.plot(X,SECT)
        show()
    return SECT

      

def Cilindrica(V,f,A,C,Paso=300, SECT=False,  R=4, Dibujar=True):
    """
    Devuelve como matriz la imagen de la UECT cilíndrica. Al poner SECT como True, hace la SECT en vez de la UECT
    (AHORA MISMO UN POCO INEXACTA LA SECT, CORREGIR SI SE QUIEREN HACER DISTANCIAS)
    """

    T=np.arange(-np.pi/2,3*np.pi/2,2*np.pi/Paso)
    j=0
    if SECT==False:
        D=[[0 for _ in range(Paso)] for _ in range(Paso)]
        for t in T:
            D[j]=CurvaEuler(V,f,A,C,[np.cos(t),np.sin(t)],Paso=Paso, Dibujar=False, R=R)
            j=j+1
        if Dibujar==True:
            im = imshow(np.flipud(D),cmap=cm.RdBu, extent=[-R,R,-np.pi/2,3*np.pi/2], aspect='auto') # drawing the function
            # adding the Contour lines with labels
            #cset = contour(Zdata.reshape(Paso,Paso),np.arange(0,8,1),linewidths=2,cmap=cm.Set2)
            #clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
            colorbar(im)
            show()
        return D
    if SECT==True:
        D=[[0 for _ in range(Paso)] for _ in range(Paso)]
        for t in T:
            D[j]=CurvaSECT(V=V,f=f,A=A,C=C, v=[np.cos(t), np.sin(t)], Paso=Paso, Dibujar=False, R=R)
            j=j+1
        if Dibujar==True:
            im = imshow(np.flipud(D),cmap=cm.RdBu, extent=[-R,R,-np.pi/2,3*np.pi/2], aspect='auto') # drawing the function
            # adding the Contour lines with labels
            #cset = contour(Zdata.reshape(Paso,Paso),np.arange(0,8,1),linewidths=2,cmap=cm.Set2)
            #clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
            colorbar(im)
            show()
        return D


def Silgen(f,A, R=4, Dibujar=False, Corte=None): 
    """
    Dibuja la silueta del complejo poligonal plano (V,f,A,C) y la guarda en siluGEN.png (solo hay que aportarle A y f). Puede dibujar un corte si se le indica en Corte un punto [x,y].
    """
    for a in A: 
        Pto=f[a[0]]
        Pto2=f[a[1]]
        Seg=np.transpose([Pto,Pto2])
        plt.plot(Seg[0],Seg[1],'-k')
    plt.gca().set_aspect('equal')

    if Corte != None:
        Pto=Corte
        Pto2=[Corte[1]+Corte[0],Corte[1]-Corte[0]]
        x=(Pto2[0]-Pto[0])*np.linspace(-4,4,2)+Pto[0]
        y=(Pto2[1]-Pto[1])*np.linspace(-4,4,2)+Pto[1]
        plt.plot(x,y, '-g')
        plt.plot(Corte[0], Corte[1], 'go')

    plt.xlim(-R,R)
    plt.ylim(-R,R)

    if Dibujar==True:
        plt.show()
    plt.savefig('siluGEN.png', transparent=True)


def Radial(V,f,A,C, Paso=200, R=4, Dibujar=True, Silueta=True, SECT=False, Corte=None):
    """
    Dibuja la UECT en proyección radial del cpp (V,f,A,C), con el paso y el radio R indicado. Se puede dibujar y optar por dibujar también la silueta del origami y un corte.
    """

    Abs=np.arange(-R,R,2*R/Paso)
    Ord=Abs
    
    ZDOS=[[0 for _ in range(Paso)] for _ in range(Paso)]
    ZTRES=[[0 for _ in range(Paso)] for _ in range(Paso)]

    if SECT==False:
        for i in range(Paso):
            for j in range(Paso):
                norm=np.sqrt(Abs[j]**2+Ord[Paso-1-i]**2)
                ZDOS[i][j]=[ValorEuler(V,f,A,C,np.array([Abs[j],Ord[Paso-1-i]])/norm,norm)]
                ZTRES[i][j]=[ValorEuler(V,f,A,C,-1*np.array([Abs[j],Ord[Paso-1-i]])/norm,-norm)]

    if SECT==True:
            Normal=Cilindrica(V,f,A,C, Paso=Paso, SECT=True, R=R, Dibujar=False)
            for i in range(Paso):
                for j in range(Paso):
                    Avance=np.sqrt(Abs[j]**2+Ord[Paso-1-i]**2)
                    if Avance>=R:
                        ZDOS[i][j]=0
                        ZTRES[i][j]=0
                    else:
                        Angulo=MiAngulo(complex(Abs[j],Ord[Paso-1-i]))
                        IndAngulo=int((Paso*(Angulo/(2*np.pi)+1/4)))
                        IndAngulo2=int((IndAngulo+Paso/2)%Paso)
                        IndAvance=int((Paso/2+(Paso-2)/(2*R)*Avance))
                        IndAvance2=int((Paso/2)+Paso/(2*R)*(-Avance))

                        ZDOS[i][j]=Normal[IndAngulo][IndAvance]
                        ZTRES[i][j]=Normal[IndAngulo2][IndAvance2]
    

    if Dibujar==True:
        m=np.amin(np.array([ZDOS,ZTRES]))
        M=np.amax(np.array([ZDOS,ZTRES]))
        plt.figure()
        fig,axarr=plt.subplots(1,2)
        #im = imshow(ZDOS,cmap=cm.RdBu,extent=[-R,R,-R,R]) # drawing the function
        #im2 =imshow(ZTRES,cmap=cm.RdBu,extent=[-R,R,-R,R])
        # adding the Contour lines with labels
        #cset = contour(ZdataDOS.reshape(Paso,Paso),np.arange(0,8,1),linewidths=2,cmap=cm.Set2)
        #clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
        #colorbar(im)
        #colorbar(im2)
        im1=axarr[0].imshow(ZDOS,cmap=cm.RdBu, vmin=m, vmax=M, extent=[-R,R,-R,R])
        im2=axarr[1].imshow(ZTRES,cmap=cm.RdBu, vmin=m, vmax=M, extent=[-R,R,-R,R])

        plt.colorbar(im1, ax=[axarr[0], axarr[1]])


        if Silueta==True:
            Silgen(f,A, R=R, Corte=Corte)
            background = Image.fromarray(np.uint8( im2.get_cmap()(im2.get_array())*255))
            overlay=Image.open("siluGEN.png")
            overlay=overlay.resize(background.size)

            background = background.convert("RGBA")
            overlay = overlay.convert("RGBA")
            
            background.paste(overlay)
            background.save("completaGEN.png")
            imgplot = plt.imshow(background)

            


            plt.show()
        else:
            show()
    return [ZDOS, ZTRES]




###EJEMPLOS###


#INICIALIZAR CPP

#[V,f,A,C]=Pajarita()
#[V,f,A,C]=HullGen(R=0.2, theta=np.pi/2, cen=True)

#VISUALIZACION CILINDRICA DE LAS UECT
#Cilindrica(V,f,A,C,R=3,Paso=100, SECT=False)

#Se recomienda visualizar la sigiuente proyección radial habiendo inicializado el segundo cpp 
#Radial(V,f,A,C, Paso=300, SECT=True, Silueta=False, R=0.6)

