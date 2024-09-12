from FuncionesUECT_V2 import *
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image


def CurvaEuler(Angulos, theta, Paso=300, Tipo='Todos', Dibujar=True):
    """
    Dado un patrón de un único vértice definido por Angulos y un ángulo theta, devuelve como vector la curva de Euler del ángulo theta con el paso aportado.
    Introducir como tipo 'Circ', 'Poli', 'Iso' o 'Todos' para obetener curvas de diferentes tipos de patrones. Si no se quieren obtener las curvas impresas por pantalla, establecer Dibujar=False.
    """
    if Tipo not in ['Todos', 'Circ', 'Poli', 'Iso']:
        print("Elegir un tipo de entre los siguientes: 'Todos', 'Circ', 'Poli', 'Iso'")
        return 0
    SumAlt=SumaAlterna(Angulos)  
    [Post, Ind, Signos]=OrdenPost(SumAlt)
    X=np.arange(0,1,1/Paso)
    if Dibujar==True:
        if Tipo=='Todos':
            fig, (ax1,ax2,ax3)=plt.subplots(3)
        else:
            fig, ax1=plt.subplots()


    if Tipo=='Circ' or Tipo=='Todos':
        D1=[0 for _ in range(Paso)]
        i=0
        for x in X:
            D1[i]=CalculoUECT(theta,x,Post, Signos)
            i=i+1
        if Dibujar==True:
            ax1.plot(X,D1)
            ax1.set_title('Tipo Circular')
            

    if Tipo=='Poli' or Tipo=='Todos':
        D2=[0 for _ in range(Paso)]
        i=0
        for x in X:
            D2[i]=CalculoUECTPoli(theta,x,Post, Signos)
            i=i+1

        if Dibujar==True:
            if Tipo=='Poli':
                ax1.plot(X,D2)
                ax1.set_title('Tipo Poligonal')
            else:
                ax2.plot(X,D2)
                ax2.set_title('Tipo Poligonal')
    if Tipo=='Iso' or Tipo=='Todos':
        D3=[0 for _ in range(Paso)]
        i=0
        for x in X:
            D3[i]=CalculoUECTIso(theta,x,Post, Signos)
            i=i+1
        if Dibujar==True:
            if Tipo=='Iso':
                ax1.plot(X,D3)
                ax1.set_title('Tipo Isósceles')
            else:
                ax3.plot(X,D3)
                ax3.set_title('Tipo Isósceles')
    
    if Dibujar==True:
        plt.gca().set_ylim(bottom=0)
        show()
    if Tipo=='Todos':
        return D1,D2,D3
    if Tipo=='Circ':
        return D1
    if Tipo=='Poli':
        return D2
    if Tipo=='Iso':
        return D3

def CurvaSECT(Angulos, theta, Paso=300, Tipo='Circ', R=1, Dibujar=True):
    """
    Devuelve la curva de la SECT en un ángulo concreto theta para el intervalo [-R,R].
    (LA INTEGRACION ES ALGO INEXACTA, CORREGIR SI SE QUIEREN HACER DISTANCIAS)
    """
    if Tipo not in ['Circ', 'Poli', 'Iso']:
        print("Elegir un tipo de entre los siguientes: 'Circ', 'Poli', 'Iso'")
        return 0
    if R<1:
        print("R ha de ser igual o mayor que 1")
        return 0
    ECT=CurvaEuler(Angulos=Angulos, theta=theta, Paso=Paso, Tipo=Tipo, Dibujar=False)



    ECT=list(ECT)
    Despues=list([0 for _ in range(int((R-1)*Paso))])
    Antes=list([0 for _ in range(int(R*Paso))])
    ECT=np.asarray(Antes+ECT+Despues)

    A=len(ECT)
    Media=np.mean(ECT)
    ECT=ECT-Media*np.array([1 for _ in range(A)])
    SECT=[0 for _ in range(A-1)]
    for i in range(A-1):
        SECT[i]=sum(ECT[0:i+1])/Paso
    SECT.append(0)

    if Dibujar==True:
        X=np.arange(-R,R,1/Paso)
        fig, ax1=plt.subplots()
        ax1.plot(X,SECT)
        ax1.set_title('SECT de tipo '+Tipo)
        show()
    return SECT

def CurvaSECT_ANTIGUA(Angulos, theta, Paso=300, Tipo='Circ', Dibujar=True):
    """
    (ANTIGUA)Devuelve la SECT circular en un ángulo concreto (ANTIGUA)
    (AHORA MISMO UN POCO INEXACTA LA INTEGRACION, CORREGIR SI SE QUIEREN HACER DISTANCIAS)
    """
    if Tipo not in ['Circ', 'Poli', 'Iso']:
        print("Elegir un tipo de entre los siguientes: 'Circ', 'Poli', 'Iso'")
        return 0
    ECT=CurvaEuler(Angulos=Angulos, theta=theta, Paso=Paso, Tipo=Tipo, Dibujar=False)

    if ECT[0:len(ECT)-1]==[0 for _ in range(len(ECT)-1)]:
        return ECT

    inicio=0
    for i in range(len(ECT)):
        if ECT[i]!=0:
            break
        else:
            inicio=inicio+1
    Media=sum(ECT)/(len(ECT)-inicio)

    ECT=ECT-Media*np.array([1 for _ in range(Paso)])
    SECT=[0 for _ in range(Paso)]
    valores=np.arange(inicio, len(ECT), step=1)
    for i in valores:
        SECT[i]=sum(ECT[inicio:i+1])/Paso
    

    if Dibujar==True:
        X=np.arange(-1,0,1/Paso)
        fig, ax1=plt.subplots()
        ax1.plot(X,SECT)
        ax1.set_title('SECT de tipo '+Tipo)
        show()
    return SECT


def UECTCompleta(Angulos,Paso=300, SECT=False, Tipo='Circ', R=1, Dibujar=True):
    """
    Devuelve como matriz la imagen cilíndrica de la UECT. Al establecer SECT como True, hace la SECT en vez de la UECT.
    """
    if Tipo not in ['Circ', 'Poli', 'Iso']:
       print("Elegir un tipo de entre los siguientes: 'Circ', 'Poli', 'Iso'")
       return 0
    T=np.arange(-pi/2,3*pi/2,2*pi/Paso)
    j=0
    if SECT==False:
        D=[[0 for _ in range(Paso)] for _ in range(Paso)]
        for t in T:
            D[j]=CurvaEuler(Angulos,t,Paso=Paso, Tipo=Tipo, Dibujar=False)
            j=j+1
        if Dibujar==True:
            im = imshow(np.flipud(D),cmap=cm.RdBu, extent=[0,1,-pi/2,3*pi/2], aspect='auto') # drawing the function
            # adding the Contour lines with labels
            #cset = contour(Zdata.reshape(Paso,Paso),np.arange(0,8,1),linewidths=2,cmap=cm.Set2)
            #clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
            colorbar(im)
            show()
        return D
    if SECT==True:
        D=[[0 for _ in range(int(2*R*Paso))] for _ in range(Paso)]
        for t in T:
            D[j]=CurvaSECT(Angulos,t,Paso=Paso, Tipo=Tipo, R=R, Dibujar=False)
            j=j+1
        if Dibujar==True:
            im = imshow(np.flipud(D),cmap=cm.RdBu, extent=[-R,R,-pi/2,3*pi/2], aspect='auto') # drawing the function
            # adding the Contour lines with labels
            #cset = contour(Zdata.reshape(Paso,Paso),np.arange(0,8,1),linewidths=2,cmap=cm.Set2)
            #clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
            colorbar(im)
            show()
        return D



def Sil(Angulos, Tipo='Circ', Dibujar=False, Corte=None): 
    """
    Guarda en silu.png la silueta del origami aportado por Angulos y Tipo. Se puede precisar también que dibuje una recta de corte, pasando un punto [x,y] a Corte.
    """
    if Tipo not in ['Circ', 'Poli', 'Iso']:
       print("Elegir un tipo de entre los siguientes: 'Circ', 'Poli', 'Iso'")
       return 0
    SumAlt=SumaAlterna(Angulos)  
    [Post, Ind, Signos]=OrdenPost(SumAlt)

    RotacionLabels=[]
    for j in range(len(Post)):
        RotacionLabels.append([np.cos(j*2*pi/len(Post))/30,np.sin(j*2*pi/len(Post))/30])

    if Tipo=='Circ':
        for i in range(len(Post)): 
            Pto=[np.cos(Post[i]),np.sin(Post[i])]
            x=Pto[0]*np.linspace(0,1,2)
            y=Pto[1]*np.linspace(0,1,2)
            plt.plot(x,y, '-k')
            plt.text(Pto[0]+RotacionLabels[i][0],Pto[1]+RotacionLabels[i][1], Ind[i])
        M=max(Post)
        t=np.linspace(0,M,100)
        xcir=np.cos(t) 
        ycir=np.sin(t)
        plt.gca().set_aspect('equal')
        plt.plot(xcir,ycir, '-k')

    if Tipo=='Poli':
        for i in range(len(Post)): 
            Pto=[np.cos(Post[i]),np.sin(Post[i])]
            x=Pto[0]*np.linspace(0,1,2)
            y=Pto[1]*np.linspace(0,1,2)
            plt.plot(x,y, '-k')
            plt.text(Pto[0]+RotacionLabels[i][0],Pto[1]+RotacionLabels[i][1], Ind[i])
            if i<len(Post)-1:
                Pto2=[np.cos(Post[i+1]),np.sin(Post[i+1])]
                x=(Pto2[0]-Pto[0])*np.linspace(0,1,2)+Pto[0]
                y=(Pto2[1]-Pto[1])*np.linspace(0,1,2)+Pto[1]
                plt.plot(x,y, '-k')
        plt.gca().set_aspect('equal')

    if Tipo=='Iso':
        M=max(Post)
        Cte=np.cos(M/2)
        for i in range(len(Post)): 
            Mult=Cte/np.cos(M/2-Post[i])
            Pto=[np.cos(Post[i])*Mult,np.sin(Post[i])*Mult]
            x=Pto[0]*np.linspace(0,1,2)
            y=Pto[1]*np.linspace(0,1,2)
            plt.plot(x,y, '-k')
            plt.text(Pto[0]+RotacionLabels[i][0],Pto[1]+RotacionLabels[i][1], Ind[i])
        Pto2=[1,0]
        x=(Pto2[0]-Pto[0])*np.linspace(0,1,2)+Pto[0]
        y=(Pto2[1]-Pto[1])*np.linspace(0,1,2)+Pto[1]
        plt.plot(x,y, '-k')
        plt.gca().set_aspect('equal')
    
    if Corte != None:
        Pto=Corte
        Pto2=[Corte[1]+Corte[0],Corte[1]-Corte[0]]
        x=(Pto2[0]-Pto[0])*np.linspace(-4,4,2)+Pto[0]
        y=(Pto2[1]-Pto[1])*np.linspace(-4,4,2)+Pto[1]
        plt.plot(x,y, '-g')
        plt.plot(Corte[0], Corte[1], 'go')

    plt.xlim(-1,1)
    plt.ylim(-1,1)

    if Dibujar==True:
        plt.show()
    plt.savefig('silu.png', transparent=True)


def UECTRadial(Angulos, Paso=300, Tipo='Circ', SECT=False, Dibujar=True, Silueta=True, Corte=None, R=1):
    """
    Devuelve la imagen radial de la UECT y la guarda en completa.png. Si se le introduce un valor [x,y] en Corte crea la recta de corte que representa ese punto. 
    Al poner SECT=True, calcula la SECT en vez de la UECT.
    """
    if Tipo not in ['Circ', 'Poli', 'Iso']:
       print("Elegir un tipo de entre los siguientes: 'Circ', 'Poli', 'Iso'")
       return 0
    Abs=np.arange(-1,1,2/Paso)
    Ord=Abs
    SumAlt=SumaAlterna(Angulos)  
    [Post, Ind, Signos]=OrdenPost(SumAlt)

    ZDOS=[[0 for _ in range(Paso)] for _ in range(Paso)]

    if SECT==False:
        if Tipo=='Circ':
            for i in range(Paso):
                for j in range(Paso):
                    ZDOS[i][j]=[CalculoUECT(MiAngulo(complex(Abs[j],Ord[Paso-1-i])),np.sqrt(Abs[j]**2+Ord[Paso-1-i]**2), Post, Signos)]
        if Tipo=='Poli':
            for i in range(Paso):
                for j in range(Paso):
                    ZDOS[i][j]=[CalculoUECTPoli(MiAngulo(complex(Abs[j],Ord[Paso-1-i])),np.sqrt(Abs[j]**2+Ord[Paso-1-i]**2), Post, Signos)]
        if Tipo=='Iso':
            for i in range(Paso):
                for j in range(Paso):
                    ZDOS[i][j]=[CalculoUECTIso(MiAngulo(complex(Abs[j],Ord[Paso-1-i])),np.sqrt(Abs[j]**2+Ord[Paso-1-i]**2), Post, Signos)]
        
    if SECT==True:
        Normal=UECTCompleta(Angulos, Paso=Paso, Tipo=Tipo, SECT=True, R=R, Dibujar=False)
        for i in range(Paso):
            for j in range(Paso):
                Avance=np.sqrt(Abs[j]**2+Ord[Paso-1-i]**2)
                if Avance>=R:
                    ZDOS[i][j]=0
                else:
                    Angulo=MiAngulo(complex(Abs[j],Ord[Paso-1-i]))
                    IndAngulo=int((Paso*(Angulo/(2*np.pi)+1/4)))
                    IndAvance=int((Paso*R*(Avance+1)))

                    ZDOS[i][j]=Normal[IndAngulo][IndAvance]


        #vals=np.arange(-1,0,1/Paso).tolist()
        #for i in range(Paso):
        #    for j in range(Paso):
        #        D=CurvaSECT(Angulos, theta=MiAngulo(complex(Abs[j],Ord[Paso-1-i])), Dibujar=False, Tipo=Tipo)
        #        ZDOS[i][j]=D[vals.index(min(vals, key=lambda x:abs(x+np.sqrt(Abs[j]**2+Ord[Paso-1-i]**2))))]

    
    if Dibujar==True:
        im = imshow(ZDOS,cmap=cm.RdBu,extent=[-1,1,-1,1]) # drawing the function
        # adding the Contour lines with labels
        #cset = contour(ZdataDOS.reshape(Paso,Paso),np.arange(0,8,1),linewidths=2,cmap=cm.Set2)
        #clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
        colorbar(im)
        
        if Silueta==True:
            Sil(Angulos, Tipo=Tipo, Corte=Corte)
            background = Image.fromarray(np.uint8( im.get_cmap()(im.get_array())*255))
            overlay=Image.open("silu.png")
            overlay=overlay.resize(background.size)

            background = background.convert("RGBA")
            overlay = overlay.convert("RGBA")
            
            background.paste(overlay)
            background.save("completa.png")
            imgplot = plt.imshow(background)
            plt.show()
        else:
            show()
    return ZDOS


def DistanciaVertice(D, p=1, q=1, Paso=300):
    """
    Devuelve la distancia (p,q) al vértice dada una UECT, D, y el paso empleado para el cálculo.
    """
    
    if p!='inf':
        norma_p = lambda t: t ** p
        norma_todo = np.vectorize(norma_p)
        D=norma_todo(D)
        Int1=[(sum(Fila)/Paso)**(1/p) for Fila in D ]
    else:
        Int1=[max(Fila) for Fila in D]

    if q!='inf':
        Int1q=[el**q for el in Int1]
        Int2=(2*pi*sum(Int1q)/Paso)**(1/q)
    else:
        Int2=max(Int1)

    return Int2

#def DistanciaInf(D, ):


def ImagenDistancias(Angulos, Paso=300, Tipo='Circ', Limsup=10):
    """
    Dado un origami Angulos, calcula una buena cantidad de distancias p,q a su vértice. La fila de arriba del todo es con q infinito.
    """
    if Tipo not in ['Circ', 'Poli', 'Iso']:
       print("Elegir un tipo de entre los siguientes: 'Circ', 'Poli', 'Iso'")
       return 0
    D=UECTCompleta(Angulos, Paso=Paso, Tipo=Tipo, Dibujar=False)

    Ps=np.linspace(1,Limsup,2*Limsup-1)
    Qs=Ps

    Dists=[[DistanciaVertice(D,p=pe, q=qu, Paso=Paso) for pe in Ps] for qu in Qs]

    A=[]
    for pe in Ps:
        A.append(DistanciaVertice(D, p=pe, q='inf', Paso=Paso))
    
    Dists.append(A)

    im = imshow(np.flipud(Dists),cmap=cm.RdBu,extent=[0.75,Limsup+0.25,0.75, Limsup+0.75])
    colorbar(im)
    show()




###Algunos ejemplos###

#PATRONES DE ORIGAMI (recordar que la suma alterna de los ángulos ha de ser 0):

#Angulos=[pi/6,pi/4,pi/6,pi/2,2*pi/3,pi/4]
#Angulos=[pi/3,pi/6,pi/6,pi/3,pi/4,pi/6,pi/4,pi/3]
#Angulos=[pi/3, 4*pi/5, 2*pi/3, pi/5]

#Angulos=[pi/3,pi/6,pi/6,pi/3,pi/4,pi/4] #Ejemplo cónico


#CURVAS DE EULER

#CurvaEuler(Angulos, 1.256, Tipo='Circ')
#CurvaSECT(Angulos, 1.256, Tipo='Circ', R=2)


#IMAGEN DE UNA UECT
#UECTCompleta(Angulos, Tipo='Circ')


#IMAGEN DE UNA SUECT
#UECTCompleta(Angulos,Tipo='Iso', SECT=True, R=1, Paso=500)


#IMAGEN RADIAL DE UNA UECT 
UECTRadial(Angulos, Tipo='Iso', Paso=500, SECT=False, Silueta=True, R=1)


#CÁLCULO DE DISTANCIAS AL VÉRTICE
#ImagenDistancias(Angulos, Paso=300, Limsup=10)

