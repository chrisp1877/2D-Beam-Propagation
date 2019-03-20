from pylab import *
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def BPM(width):
    x_dom = 100
    y_dom = 50
    dx = 1
    dy = 0.1
    steps = x_dom/dx
    
    lam_0 = 1 #um
    n_air = 1
    n_wg = 1.53
    eps_ave = 1 * sqrt(n_air) + 1 * sqrt(n_wg)
    deps = sqrt(n_wg) - eps_ave
    lam = lam_0 / eps_ave**2
    d_wg = 0.2
    x = arange(0, x_dom, dx) ### 0 : dx : x_dom
    y = arange(-0.5 * y_dom, 0.5 * y_dom, dy) ### -0.5 * y_dom : dy : 0.5 * y_dom
    l_y = len(y)
    
    k_0 = 2 * pi / lam
    dky = pi / y_dom * (l_y - 1) / l_y
    k_y = dky * arange(-l_y / 2, l_y/2 ) #l_y/2 - 1?
    l_wg = int(rint(d_wg / (y_dom * dy)))
    l_air = int(rint(0.5 * (y_dom - d_wg) / (y_dom * dy)))
    k_y = k_y * append(ones((1, l_air)), append(ones((1, l_wg)), ones((1, int(l_y-l_air-l_wg)))))
    A = 1
    w = width
    E = A * exp(-(y/w)**2) 
    E_final = []
    
    for j in range(int(steps)):
        Eft = fftshift(fft(fftshift(E)))
        D = exp(1j * sqrt(k_0**2 - k_y**2) * dx * 0.5)
        Eft_D = Eft * D
        E = ifftshift(ifft(ifftshift(Eft_D)))
        
        W = exp(1j * k_0 * deps * dx / (2 * eps_ave))
        E_W = E * W
        
        Eft_W = fftshift(fft(fftshift(E_W)))
        D = exp(1j * sqrt(k_0**2 - k_y**2) * dx * 0.5)
        Eft_WD = Eft_W * D
        E = ifftshift(ifft(ifftshift(Eft_WD)))
        
        E[0] = 0
        E[-1] = 0
        if j == 0:
            E_final = E
        else:
            E_final = vstack((E_final, E))
#    return E_final

    [Y,X] = meshgrid(y, x)
    pcolormesh(Y, X, abs(E_final)**2)
        
    
if __name__ == "__main__":
    BPM(1)