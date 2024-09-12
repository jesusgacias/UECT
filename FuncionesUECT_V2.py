import numpy as np
import matplotlib.pyplot as plt
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from scipy import array
import math 

pi=np.pi


def SumaAlterna(Alphas):
    """
    Dado un conjunto de ángulos Alphas, calcula su suma alterna y después los rota para que el menor ángulo se encuentre alineado con el cero
    """
    S=[0]*len(Alphas)
    S[0]=0
    for i in range(len(Alphas)-1):
        S[i+1]=S[i]+(-1)**(i)*Alphas[i]
    S=np.sum([S,[-min(S)]*len(Alphas)], axis=0)
    return S

def OrdenPost(Sumas):
    """
    Ordena los ángulos finales de cada pliegue. Devuelve: 1-Los ángulos ordenados. 2-Los índices de los ángulos según el orden post-plegado. 3-Las signaturas de los vértices ordenados.
    """
    Post=sorted(Sumas)
    IndicesPost=sorted(range(len(Sumas)), key=lambda k: Sumas[k])
    Signos=[(-1)**n for n in IndicesPost]
    return Post, IndicesPost, Signos

def HallarLugar(Theta, Post):
    """
    Dado un ángulo y los ángulos en orden post-plegado de un patrón, devuelve la posición del ángulo entre todos ellos.
    """
    j=0
    for i in range(len(Post)):
        if Theta>Post[i]:
            j=j+1
        else:
            break
    return j

def rindex(lst, value):
    """
    Devuelve la última posición de un elemento en una lista.
    """
    lst.reverse()
    i = lst.index(value)
    lst.reverse()
    return len(lst) - i - 1

def Capas(Theta, Post, Signos):
    """
    Calcula la función de número de capas en un ángulo Theta. Requiere el ángulo, el orden post-plegado de los ángulos, y las signaturas de los pliegues.
    """
    if Theta in Post:
        j=Post.index(Theta)
        k=rindex(Post, Theta)
        C=sum(Signos[:j])-sum(Signos[k+1:])
        return C
    else:
        j=HallarLugar(Theta,Post)
        if j==0 or j==len(Post):
            return 0
        C=sum(Signos[:j])-sum(Signos[j:])
        return C

def CapasPoli(Theta, j, Post, Signos):
    """
    Devuelve el valor de la función de capas en la esquina más cercana a Theta, para su uso en el patrón poligonal. j es el lugar que ocupa Theta en el array Post, de ángulos en orden post-plegado.
    """
    if j==0 or j==len(Post):
        return 0
    return Capas(Post[MasCercana(Theta,j,Post)], Post, Signos)

def MasCercana(Theta, j, Post):
    """
    Calcula el pliegue más cercano a Theta, que se encuentra en posición j en el array Post de orden post-plegado.
    """
    if j==0:
        return 0
    if j==len(Post):
        return len(Post)-1
    k=rindex(Post, Post[j-1])
    if Theta-Post[j-1]<Post[k+1]-Theta:
        return j-1
    else:
        return k+1

def Indicadora(Theta, Limite):
    """
    Función indicadora {Theta>=Limite}
    """
    if Theta<Limite:
        return 0
    if Theta>=Limite:
        return 1

def Indicadora2(x, LimiteInf):
    """
    Función indicadora x>LimiteInf
    """
    if x>LimiteInf:
        return 1
    else:
        return 0

def SignosRel(Theta, Post, Signos):
    """
    Devuelve las signaturas relativas de unos ángulos en orden Post y con signaturas Signos, con respecto a Theta. 
    """
    SigRel=np.array(Signos)
    if Theta in Post:
        i=Post.index(Theta)
        j=rindex(Post, Theta)
        SigRel[i:j+1]=[0]*(j+1-i)
        SigRel[j+1:]=-SigRel[j+1:]
        return SigRel
    else:
        k=HallarLugar(Theta, Post)
        SigRel[k:]=-SigRel[k:]
        return SigRel

def MiAngulo(z):
    """
    Toma un ángulo z entre -pi y pi y lo redefine entre -pi/2 y 3*pi/2.
    """
    A=np.angle(z)
    if A>=-pi/2:
        return A
    else:
        return A+2*pi

def CalculoUECT(theta, x, Post, Signos):
    """
    Calcula la UECT en ángulo theta (entre -pi/2 y 3pi/2), y a altura x real de un patrón circular con orden post-plegado Post y signaturas Signos. Devuelve 1 para cualquier x negativo.
    """
    if x<=0:
        return 1
    if x>1:
        return 0
    M=max(Post)
    if theta<-pi or theta>M+pi/2:
        return 0
    SRel=SignosRel(theta,Post,Signos)
    res=0
    #if theta >-pi/2:
    for i in range(len(Post)):
        res=res+SRel[i]*Indicadora2(x, np.cos(theta-Post[i]))
    #else:
    #    for i in range(len(Post)):
    #        res=res+Signos[i]*DaLaVuelta(theta,Post[i])*Indicadora(x, -np.cos(theta-Post[i]))
    return res


def CalculoUECTPoli(theta,x,Post,Signos):
    """
    Calcula la UECT en theta (entre -pi/2 y 3pi/2), y x real de un patrón poligonal con orden post-plegado Post y signaturas Signos. Devuelve 1 para cualquier x negativo.
    """
    if x<=0:
        return 1
    if x>1:
        return 0
    M=max(Post)
    if theta<-pi or theta>M+pi/2:
        return 0
    j1=HallarLugar(theta,Post)
    #if j1==0 or j1==len(Post):
    #    return CalculoUECT(theta,x,Post,Signos)

    p=MasCercana(theta,j1,Post)
    if Indicadora(np.cos(theta-Post[p]),x)==0:
        return 0
    SRel=SignosRel(theta, Post, Signos)

    res=0
    for i in range(len(Post)):
        res=res+SRel[i]*Indicadora2(x, np.cos(theta-Post[i]))

    return res


def CalculoUECTIso(theta,x,Post,Signos):
    """
    Calcula la UECT en theta (entre -pi/2 y 3pi/2), y x real de un patrón isósceles con orden post-plegado Post y signaturas Signos. Devuelve 1 para cualquier x negativo.
    """
    if x<=0:
        return 1
    if x>1:
        return 0
    M=max(Post)
    if theta<-pi or theta>M+pi/2:
        return 0
    Mitad=M/2
    Cte=np.cos(Mitad)
    res=0
    for i in range(len(Post)):
        res=res+Signos[i]*Indicadora(Cte*np.cos(Post[i]-theta)/np.cos(Mitad-Post[i]),x)
    if theta>Mitad:
        res=-res
    return res