import matplotlib.pyplot as plt
import numpy as np


def fftIndgen(n):
    a = range(0, int(n/2+1))
    b = range(1, int(n/2))
    b =list(b)
    a =list(a)
    b.reverse()
    b = [-i for i in b]
    return a + b


rady=10.
radx=10.0
from math import *
def Pk_break(kx,ky,beta):
    rad=np.sqrt((kx/1.75)**2+(1.75*ky)**2)
    theta=-1.5*np.pi/4*1e-3
    theta=0
    kxr=(kx*cos(theta)+ky*sin(theta))
    kyr=(-kx*sin(theta)+ky*cos(theta))
    rad=np.sqrt((kxr)**2+(kyr)**2)
    xf=radx**(-3)/radx**-(5./3)
    xf=0.5/radx**-(5./3)
    if rad<0.1:
        return 0
    if rad<radx:
        return (0.5)
    else:
        return (xf*(rad)**(beta))

def Pk3(kx,ky,kz,beta):
    rad=np.sqrt((kx/1.75)**2+(1.75*ky)**2)
    theta=-1.5*np.pi/4-(kz/100.)*pi/2
    kxr=(kx*cos(theta)+ky*sin(theta))
    kyr=(-kx*sin(theta)+ky*cos(theta))
    rad=np.sqrt((kxr/1.25)**2+(1.25*kyr)**2+(kz*6.)**2)
    xf=radx**(-3)/radx**-(5./3)
    xf=0.5/radx**-(5./3)
    if rad<0.1:
        return 0.5
    if rad<radx:
        #return rad**(-3)
        return (0.5)
    else:
        return (xf*(rad)**beta)
        
def gaussian_random_field(beta=-4, size = 100):
    def Pk2(kx, ky):
        if kx == 0 and ky == 0:
            return 0.0
        #return np.sqrt(Pk(np.sqrt(ky**2)))
        return np.sqrt(Pk_break(kx,ky,beta))
    nreal=np.random.normal(size = (size, size))
    noise = np.fft.fft2(nreal)
    #print nreal.max(), nreal.min()
    nreal[nreal>5.5]=5.5
    nreal[nreal<-5.5]=5.5
    amplitude = np.zeros((size,size))
    for i, kx in enumerate(fftIndgen(size)):
        for j, ky in enumerate(fftIndgen(size)):            
            amplitude[i, j] = Pk2(kx, ky)
    return np.fft.ifft2(noise * amplitude)
