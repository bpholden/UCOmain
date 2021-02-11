#!/usr/bin/env python
from __future__ import print_function

import sys
try:
    sys.path.append('/u/user/bjscripts/lib/python')
except:
    print("Error appending to path")
    
import numpy as np
import time
import pyfits
from scipy.optimize import leastsq

win = [14,14]
try:
    lines = np.genfromtxt('apflines.txt', dtype=float)
except:
    lines = np.genfromtxt('/mir4/apflines.txt', dtype=float)

def gauss(x, p):
    A, mu, sigma, C = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + C

errfunc = lambda p, x, y: (gauss(x, p) - y)**2

def tharfwhm(image,lines=lines, doplot=False, silent=False, pltname='fwhm.png'):
    
    imdat = pyfits.getdata(image).transpose()

    if not silent: print( "Pixel xy guess  FWHM    peak ADU")

    allfwhm = []
    allpeaks = []
    x = []
    y = []
    for i in range(lines.shape[0]):
        lines[i] = np.round(lines[i])
        imbox = imdat[lines[i,0]-win[0]/2:lines[i,0]+win[0]/2, lines[i,1]-1-win[1]/2:lines[i,1]-1+win[1]/2]

        if np.max(imbox) < 500:
            if not silent: print( "Line flux too small!")
            continue
        if np.max(imbox) > 40000:
            if not silent: print( "Line too close to saturation!")
            continue
        
        bins, edges = np.histogram(imbox, bins=50)
        mode = edges[np.argmax(bins)]
        imflat = np.sum(imbox-mode,axis=0)
        
        xfit = np.arange(-win[0]/2, win[0]/2, 1)
        p0 = [np.max(imflat)-np.min(imflat), 0.0, 2.5, np.min(imflat)]
        #pfin,covar = curve_fit(gauss, xfit, imflat, p0=p0)
        pfin = leastsq(errfunc, p0, args=(xfit, imflat), maxfev=10000)[0]
        yfit = gauss(xfit, pfin)

        fwhm = 2*np.sqrt(2*np.log(2))*pfin[2]
        if fwhm < 0:
            if not silent: print( "Negative FWHM, bad fit!")
            continue
        
        allfwhm.append(fwhm)
        x.append(lines[i,0])
        y.append(lines[i,1])
        allpeaks.append(np.max(imbox))
        
        if False and doplot:
            import pylab as pl
            pl.ion()
            pl.clf()
            xplot = np.linspace(min(xfit), max(xfit), 200)
            yfit = yfit = gauss(xplot, pfin)
            pl.plot(xfit, imflat,'k-o', markersize=8)
            pl.plot(xplot, yfit, 'b-')
            pl.ylabel('ADU')
            pl.xlabel('$\\Delta$pixel')
            pl.draw()
            time.sleep(0.2)
        
        if not silent: print( "%s\t%6.4f\t%d" % (lines[i], fwhm, np.max(imbox)))

    if not silent: print( "-------------------------------\n")

    if doplot:
        import pylab as pl
        pl.ioff()
        from scipy.interpolate import griddata
        x = np.array(x)
        y = np.array(y)
        z = np.array(allfwhm)
        xi = np.linspace(min(x), max(x), 1000)
        yi = np.linspace(min(y), max(y), 1000)
        zi = griddata((x,y),z, (xi[None,:],yi[:,None]), method='nearest')

        pl.imshow(np.log10(imdat), vmin=np.log10(1000), vmax=np.log10(1500), aspect=7)
        #pl.contourf(yi,xi,zi, 50)
        pl.scatter(y,x,c=z,s=100, vmin=2, vmax=5)
        pl.xlim(min(y), max(y))
        pl.ylim(min(x), max(x))
        pl.colorbar()
        pl.title("median FWHM = %4.2f" % np.median(allfwhm))
        pl.xlabel('row [px]')
        pl.ylabel('col [px]')
        pl.savefig(pltname)
        #pl.show()
        pl.clf()
    
    return np.median(allfwhm), np.max(allpeaks)

if __name__ == '__main__':
    import sys
    
    try:
        image = sys.argv[1]
    except IndexError:
        print( "usage: ./tharfwhm.py [image file]")

    fwhm,flux = tharfwhm(image, lines, doplot=False, pltname=sys.argv[1]+'.foc.png')
    print( "Median FWHM = %6.4f" % fwhm)
