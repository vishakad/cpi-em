#!/usr/bin/env python2.7
#import numpy as np
import argparse
from scipy.stats import gaussian_kde, lognorm, norm, gamma
from scipy.special import gammaln
from scipy.optimize import minimize
from pandas import DataFrame, read_table
from numpy import append, array, copy, exp, isnan, log, median, sqrt, sum, zeros, inf, pi
from numpy.random import random as np_random
from os.path import isfile

def indLognormalPdf( logX, logY, halfLogSqX, halfLogSqY, logConstX, logConstY ):
    logxpdf = -logConstX - logX - halfLogSqX
    logypdf = -logConstY - logY - halfLogSqY
    logxpdf[isnan( logxpdf )] = -80
    logypdf[isnan( logypdf )] = -80

    return logxpdf + logypdf

def logGaussianPdf( x, muX, sigmaX, y, muY, sigmaY  ):
    logxpdf = -0.5 * log( 2*pi*sigmaX*sigmaX ) - 0.5 * (x - muX) * (x - muX)/(sigmaX*sigmaX) 
    logypdf = -0.5 * log( 2*pi*sigmaY*sigmaY ) - 0.5 * (y - muY) * (y - muY)/(sigmaY*sigmaY)
    logxpdf[isnan( logxpdf )] = -80
    logypdf[isnan( logypdf )] = -80

    return logxpdf + logypdf

def logGammaPdf( x, logx, shapeX, scaleX, y, logy, shapeY, scaleY  ):
    logxpdf = -gammaln( shapeX ) - shapeX*log(scaleX) + (shapeX - 1)*logx - x/scaleX
    logypdf = -gammaln( shapeY ) - shapeY*log(scaleY) + (shapeY - 1)*logy - y/scaleY
    logxpdf[isnan( logxpdf )] = -80
    logypdf[isnan( logypdf )] = -80

    return logxpdf + logypdf

def gammaEstep( pi0, xlist, ylist, shapeXind, scaleXind, shapeXint, scaleXint, shapeYind, scaleYind, shapeYint, scaleYint ):
    x = xlist[0]
    logx = xlist[1]
    y = ylist[0]
    logy = ylist[1]
    weights = zeros( (2,len(logx)) )

    logIndPdfVals = logGammaPdf( x, logx, shapeXind, scaleXind, y, logy, shapeYind, scaleYind )
    logIntPdfVals = logGammaPdf( x, logx, shapeXint, scaleXint, y, logy, shapeYint, scaleYint )

    lpi0 = log( pi0 )
    lpi1 = log( 1 - pi0 )
    indPopnMask = (lpi0 + logIndPdfVals - lpi1 - logIntPdfVals) > 8
    intPopnMask = (lpi0 + logIndPdfVals - lpi1 - logIntPdfVals) < -8
    neither = ~(indPopnMask | intPopnMask)
    relIndPdfVals = exp( logIndPdfVals[neither] - logIntPdfVals[neither] )

    weights[0][indPopnMask] = 1
    weights[0][intPopnMask] = 0

    weights[0][neither] = pi0 * relIndPdfVals / (pi0 * relIndPdfVals + (1-pi0) )
    #temp = pi0 * exp( logIndPdfVals )/ (pi0 * exp( logIndPdfVals ) + (1-pi0)*exp( logIntPdfVals ) )
    weights[1] = 1 - weights[0]
    return weights

def gammaMstep( params, weights, xlist, ylist ):
    shapeXind, scaleXind, shapeXint, scaleXint, shapeYind, scaleYind, shapeYint, scaleYint = params
    x = xlist[0]
    logx = xlist[1]
    y = ylist[0]
    logy = ylist[1]

    logIndPdfVals = logGammaPdf( x, logx, shapeXind, scaleXind, y, logy, shapeYind, scaleYind )
    logIntPdfVals = logGammaPdf( x, logx, shapeXint, scaleXint, y, logy, shapeYint, scaleYint )

    logLikelihood = -sum( weights[0] * logIndPdfVals + weights[1] * logIntPdfVals )
    return logLikelihood

def lognormalEstep( pi0, logx, logy, shapeXind, scaleXind, shapeXint, scaleXint, shapeYind, scaleYind, shapeYint, scaleYint ):
    weights = zeros( (2,len(logx)) )
    logXind = logx - log(scaleXind)
    logYind = logy - log(scaleYind)
    logXint = logx - log(scaleXint)
    logYint = logy - log(scaleYint)

    halfLogSqXind = 0.5 * logXind * logXind/(shapeXind*shapeXind)
    halfLogSqYind = 0.5 * logYind * logYind/(shapeYind*shapeYind)
    halfLogSqXint = 0.5 * logXint * logXint/(shapeXint*shapeXint)
    halfLogSqYint = 0.5 * logYint * logYint/(shapeYint*shapeYint)

    logConstXind = log( scaleXind * shapeXind * sqrt( 2*pi ) )
    logConstYind = log( scaleYind * shapeYind * sqrt( 2*pi ) )
    logConstXint = log( scaleXint * shapeXint * sqrt( 2*pi ) )
    logConstYint = log( scaleYint * shapeYint * sqrt( 2*pi ) )

    logIndPdfVals = indLognormalPdf( logXind, logYind, halfLogSqXind, halfLogSqYind, logConstXind, logConstYind )
    logIntPdfVals = indLognormalPdf( logXint, logYint, halfLogSqXint, halfLogSqYint, logConstXint, logConstYint )

    lpi0 = log( pi0 )
    lpi1 = log( 1 - pi0 )
    indPopnMask = (lpi0 + logIndPdfVals - lpi1 - logIntPdfVals) > 8
    intPopnMask = (lpi0 + logIndPdfVals - lpi1 - logIntPdfVals) < -8
    neither = ~(indPopnMask | intPopnMask)
    relIndPdfVals = exp( logIndPdfVals[neither] - logIntPdfVals[neither] )

    weights[0][indPopnMask] = 1
    weights[0][intPopnMask] = 0

    weights[0][neither] = pi0 * relIndPdfVals / (pi0 * relIndPdfVals + (1-pi0) )
    weights[1] = 1 - weights[0]

    return weights

def lognormalMstep( params, weights, logx, logy ):
    shapeXind, scaleXind, shapeXint, scaleXint, shapeYind, scaleYind, shapeYint, scaleYint = params
 
    logXind = logx - log(scaleXind)
    logYind = logy - log(scaleYind)
    logXint = logx - log(scaleXint)
    logYint = logy - log(scaleYint)
 
    halfLogSqXind = 0.5 * logXind * logXind/(shapeXind*shapeXind)
    halfLogSqYind = 0.5 * logYind * logYind/(shapeYind*shapeYind)
    halfLogSqXint = 0.5 * logXint * logXint/(shapeXint*shapeXint)
    halfLogSqYint = 0.5 * logYint * logYint/(shapeYint*shapeYint)
 
    logConstXind = log( scaleXind * shapeXind * sqrt( 2*pi ) )
    logConstYind = log( scaleYind * shapeYind * sqrt( 2*pi ) )
    logConstXint = log( scaleXint * shapeXint * sqrt( 2*pi ) )
    logConstYint = log( scaleYint * shapeYint * sqrt( 2*pi ) )
 
    logIndPdfVals = indLognormalPdf( logXind, logYind, halfLogSqXind, halfLogSqYind, logConstXind, logConstYind )
    logIntPdfVals = indLognormalPdf( logXint, logYint, halfLogSqXint, halfLogSqYint, logConstXint, logConstYint )
 
    #Have to multiply by -1 in order to maximize.
    logLikelihood = -sum( weights[0] * logIndPdfVals + weights[1] * logIntPdfVals )
 
    return logLikelihood
    #return logLikelihood, max(logIndPdfVals), max(logIntPdfVals)

def gaussianEstep( pi0, x, y, muXind, sigmaXind, muXint, sigmaXint, muYind, sigmaYind, muYint, sigmaYint ):
    weights = zeros( (2,len(x)) )
    logIntPdfVals = logGaussianPdf( x, muXint, sigmaXint, y, muYint, sigmaYint )
    logIndPdfVals = logGaussianPdf( x, muXind, sigmaXind, y, muYind, sigmaYind )

    lpi0 = log( pi0 )
    lpi1 = log( 1 - pi0 )
    indPopnMask = (lpi0 + logIndPdfVals - lpi1 - logIntPdfVals) > 8
    intPopnMask = (lpi0 + logIndPdfVals - lpi1 - logIntPdfVals) < -8
    neither = ~(indPopnMask | intPopnMask)
    relIndPdfVals = exp( logIndPdfVals[neither] - logIntPdfVals[neither] )

    weights[0][indPopnMask] = 1
    weights[0][intPopnMask] = 0

    weights[0][neither] = pi0 * relIndPdfVals / (pi0 * relIndPdfVals + (1-pi0) )
    weights[1] = 1 - weights[0]

    return weights

def gaussianMstep( params, weights, x, y, ):
    muXind, sigmaXind, muXint, sigmaXint, muYind, sigmaYind, muYint, sigmaYint  = params
    logIntPdfVals = logGaussianPdf( x, muXint, sigmaXint, y, muYint, sigmaYint )
    logIndPdfVals = logGaussianPdf( x, muXind, sigmaXind, y, muYind, sigmaYind )

    #-1 because we use scipy.optimize.minimize
    logLikelihood = -sum( weights[0] * logIndPdfVals + weights[1] * logIntPdfVals )

    return logLikelihood

    #return logLikelihood, max(logIndPdfVals), max(logIntPdfVals)

def deterministicMstep( x, y, lnXind, lnYind, lnXint, lnYint, weights, ltype='lognormal' ):
    logx = log( x )
    logy = log( y )

    params = array([])
    for p in [lnXind,lnYind,lnXint,lnYint]:
        params = append( params, p )

    if ltype == 'lognormal':
        mStep = lognormalMstep
        xdata = logx
        ydata = logy
    elif ltype == 'gamma':
        mStep = gammaMstep
        xdata = [x,logx]
        ydata = [y,logy]
    elif ltype == 'gaussian':
        mStep = gaussianMstep
        xdata = x
        ydata = y

    res = minimize( mStep, params.tolist(), args=(weights,xdata,ydata), tol=1e-6, method='Powell' )

    qVal = -res['fun'] 
    newParams = res['x']

    return newParams, qVal

def getEmInitialParams( x, y, ltype ):
    if ltype == 'lognormal':
        func = lognorm
    elif ltype == 'gaussian':
        func = norm
    elif ltype == 'gamma':
        func = gamma

    if ltype in ['lognormal','gamma']: 
        shapeX, locX, scaleX = func.fit( x )
        shapeY, locY, scaleY = func.fit( y )
    elif ltype == 'gaussian':
        shapeX, scaleX = func.fit( x )
        shapeY, scaleY = func.fit( y )

    return array([shapeX, scaleX]), array([shapeY, scaleY])

def validateInput( inputFileName, ltype, nEMitr ):
    emptyInput = DataFrame()
    if isfile( inputFileName ):
        input = read_table( inputFileName, sep="\t", header=None )
        numColumns = input.shape[1]
    else:
        print('Input file {} does not exist'.format( inputFileName)) 
        return emptyInput
        
    if numColumns > 2:
        print('Number of columns recognized : {}'.format( numColumns ))
        print('Input file must contain exactly two columns.')
        return emptyInput

    if numColumns < 2:
        print('First entry in file processed as : {}'.format( input[0][0] ))
        print('Number of columns recognized : {}'.format( numColumns ))
        print('Ensure that\n(1)Peak intensity entries are delimited by <TAB>\n(2)The first line contains no header information\n(3)Input file contains exactly two columns.')
        return emptyInput

    if input.dtypes[0] == 'O':
        print('First column has some non-numeric data. Ensure that all entries input are numeric.')
        return emptyInput

    if input.dtypes[1] == 'O':
        print('Second column has some non-numeric data. Ensure that all entries input are numeric.')
        return emptyInput

    if ltype not in ['lognormal','gaussian','gamma']:
        print('Mixture type specified : {}'.format( ltype ))
        print('Distribution mixture type must be lognormal, gaussian or gamma.')
        return emptyInput

    if nEMitr < 0 :
        print('Max. number of EM iterations specified : {}'.format( nEMitr ))
        print('Max. number of EM iterations must be a positive  number.')
        return emptyInput

    return input

def run( inputFileName, ltype="lognormal", outputFileName="", nEMitr=10000 ):
    xy = validateInput( inputFileName, ltype, nEMitr )
    if len(xy) == 0:
        print('No output produced.')
        return None

    x = xy[0].values
    y = xy[1].values

    pX, pY = getEmInitialParams( x, y, ltype )
    pi0 = np_random()
    params = []

    logx = log( x )
    logy = log( y )
    if ltype == 'lognormal':
        evalEstep = lognormalEstep
        xEstep = logx
        yEstep = logy
    elif ltype =='gamma':
        evalEstep = gammaEstep
        xEstep = [x,logx]
        yEstep = [y,logy]
    elif ltype == 'gaussian':
        evalEstep = gaussianEstep
        xEstep = x
        yEstep = y

    theta = {}
    theta['pi0'] = pi0      #Fraction of non-cooperatively bound peak pairs.
    theta['lnXA'] = pX      #Parameters of f_A^X
    theta['lnXB'] = pX      #Parameters of f_B^X
    theta['lnYA'] = pY      #Parameters of f_A^Y
    theta['lnYB'] = pY      #Parameters of f_B_Y
    qValPrev = -inf
    initTheta = {}
    for key in theta.keys():
        initTheta[key] = copy( theta[key] )

    emItr = 1
    #initTheta = copy.deepcopy( theta )
    print('Fitting a {} mixture to data'.format( ltype ))
    print('Maximum number of EM iterations to run : {}'.format( nEMitr ))
    converged = False
    while emItr <= nEMitr:
        weights = evalEstep( pi0, xEstep, yEstep, theta['lnXA'][0], theta['lnXA'][1], theta['lnXB'][0], theta['lnXB'][1], theta['lnYA'][0], theta['lnYB'][1], theta['lnYB'][0], theta['lnYB'][1] )

        params, qVal = deterministicMstep( x, y, theta['lnXA'],theta['lnYA'],theta['lnXB'],theta['lnYB'],weights, ltype=ltype) 
        theta['lnXA'] = params[:2]
        theta['lnYA'] = params[2:4]
        theta['lnXB'] = params[4:6]
        theta['lnYB'] = params[6:8]

        pi0 = sum( weights[0] )*1.0/x.shape[0]
        if abs( (qVal - qValPrev)/qValPrev ) < 1e-6 and emItr > 1: #Convergence criterion
            converged = True
            break

        emItr += 1
        qValPrev = qVal

    if converged:
        print('EM converged in {} iterations. Q_2 = {}'.format( emItr, qVal ))
    else:
        print('EM did not converge. Try increasing the maximum number of iterations allowed. Output will still be provided.')

    """Medians of f_A^X and f_B^X are computed. Since cooperatively bound
    primary TFs have a lower peak intensity on average, if med( f_A^X ) < med(
    f_B^X), then, f_1^X = f_A^X and f_1^Y = f_A^Y. If med( f_A^X ) > med( f_B^X
    ), then, f_1^X = f_B^X and f_1^Y = f_B^Y. 
    """
    if ltype == 'lognormal':
        pXindMedian = theta['lnXA'][1]
        pXintMedian = theta['lnXB'][1] 
    elif ltype == 'gamma':
        xIndSamps = gamma.rvs( theta['lnXA'][0], 0, theta['lnXA'][1], size=100000 )
        pXindMedian = median( xIndSamps )
        xIntSamps = gamma.rvs( theta['lnXB'][0], 0, theta['lnXB'][1], size=100000 )
        pXintMedian = median( xIntSamps )
    elif ltype == 'gaussian':
        pXindMedian = theta['lnXA'][0]
        pXintMedian = theta['lnXB'][0]

    if pXintMedian < pXindMedian:
        probInt = weights[1]
    else:
        probInt = weights[0]

    if len(outputFileName) == 0:
        outputFileName = inputFileName + ".cpi-em"

    print('Writing output to {}'.format( outputFileName ))
    print('The third column is the probability of a cooperative interaction being present for each entry in the input file')
    xy.loc[:,2] = probInt
    xy.to_csv( outputFileName, sep="\t", index=False, header=None )

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument( "-i", "--input-file", help="Name of input file. Input should be a tab-separated file with two columns consisting of primary TF peak intensities (first column) and partner TF peak intensities (second column). The partner TF is assumed to be TF that cooperatively aids in binding the primary TF to DNA.", nargs=1, metavar="<input file>", required=True )
    parser.add_argument( "-m","--mixture", help="Distribution to be used to model peak intensities. Can be lognormal, gaussian or gamma.  Default : lognormal", nargs=1, metavar="<dist>", default="lognormal", choices=["lognormal","gamma","gaussian"])
    parser.add_argument( "-o", "--output-file", help="File to which output will be written. Default : <input file>.cpi-em ", nargs=1, metavar="<output file>", default="" )
    parser.add_argument( "-n", "--num-iterations", help="Maximum number of EM iterations to execute. Default : 10000", nargs=1, metavar="<NUM>", default=[10000], type=int )

    args = vars(parser.parse_args())
    inputFile = args['input_file'][0]
    ltype = args['mixture'][0]
    nEMitr = args['num_iterations'][0]
    outputFile = args['output_file'][0]

    run( inputFile, ltype=ltype, outputFileName=outputFile, nEMitr=nEMitr )

if __name__ == "__main__":
    main()
