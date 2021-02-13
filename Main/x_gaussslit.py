from __future__ import print_function

import numpy
from scipy import interpolate

def x_gaussslit(width,height,xcenter,ycenter):
    

    psf = numpy.array([1.000,   .995,   .985,   .971,   .954, .933,   .911,   .886,   .860,   .833, .804,   .774,   .743,   .713,   .682, .651,   .620,   .594,   .559,   .529, .500,   .471,   .443,   .417,   .391, .366,   .342,   .319,   .297,   .276, .256,   .237,   .218,   .202,   .187, .172,   .158,   .145,   .132,   .122, .113,   .104,   .097,   .089,   .082, .077,   .072,   .065,   .059,   .057, .052,   .049,   .046,   .042,   .039, .037,   .034,   .032,   .029,   .027, .026,   .024,   .023,   .021,   .019, .018,   .017,   .017,   .016,   .016, .015,   .014,   .013,   .012,   .011, .010,   .010,   .009,   .009,   .008, .008,   .007,   .007,   .006,   .006, .005,   .005,   .005,   .004,   .004, .004,   .004,   .003,   .003,   .003, .003,   .003,   .002,   .002,   .002])
    psfrads = numpy.linspace(0,psf.size-1,num=psf.size)

    
    width = width * 20
    height = height * 20
    xoff = 40*xcenter
    yoff = 40*ycenter

    yvs, xvs = numpy.mgrid[-99:100,-99:100]
    
    xvs -= xoff
    yvs -= yoff

    rs = xvs**2 + yvs**2
    rs = numpy.sqrt(rs)
    
    dx = numpy.abs(xvs)
    dy = numpy.abs(yvs)

    inside = numpy.where((dy < height) & (dx < width))
    outside = numpy.where((dy > height) | (dx > width))
    edgey = numpy.where((dy == height) & (dx <= width))
    edgex = numpy.where((dy <= height) & (dx == width))
    
    fluxes = numpy.zeros_like(rs)
    mod = interpolate.splrep(psfrads,psf,s=0)
    fluxes[numpy.where(rs<=99)] = interpolate.splev(rs[numpy.where(rs<=99)],mod,der=0)

    insideflux = fluxes[inside].sum()
    outsideflux = fluxes[outside].sum()
    edgexflux = fluxes[edgex].sum()
    edgeyflux = fluxes[edgey].sum()

    insideflux += 0.5*(edgexflux + edgeyflux)
    outsideflux += 0.5*(edgexflux + edgeyflux)
    
    fraction = insideflux / (insideflux + outsideflux)
    
    return fraction

if __name__ == "__main__":

    fwhm = 14.6
    width = 9.26
    height = 27.78
    #    width = 74.07
    #    height = 74.07
    print (x_gaussslit(width/fwhm,height/fwhm,0,0))
    
