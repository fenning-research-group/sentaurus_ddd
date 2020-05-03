import numpy as np
from scipy.optimize import OptimizeResult
from typing import Callable
import scipy.linalg as LA
from scipy.stats.distributions import t
import warnings


def confint(n: int, pars: np.ndarray, pcov: np.ndarray, confidence: float = 0.95,
            **kwargs):
    """
    This function returns the confidence interval for each parameter

    Parameters
    ----------
    n : int
        The number of data points
    pars : [double]
        The array with the fitted parameters
    pcov : [double]
        The covariance matrix
    confidence : float
        The confidence interval
    
    Returns
    -------
    np.ndarray
        ci: The matrix with the confindence intervals for the parameters
        
    Note:
        Adapted from 
        http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/
        Copyright (C) 2013 by John Kitchin.
        https://kite.com/python/examples/702/scipy-compute-a-confidence-interval-from-a-dataset
    
    """
    is_log = kwargs.get('is_log',False)
    from scipy.stats.distributions import t
    
    if is_log:
        p = np.power(10,pars)
        pcov = np.power(10, pcov)

    p = len(pars)  # number of data points
    dof = max(0, n - p)  # number of degrees of freedom

    # Quantile of Student's t distribution for p=(1 - alpha/2)
    # tval = t.ppf((1.0 + confidence)/2.0, dof) 
    alpha = 1.0 - confidence
    tval = t.ppf(1.0 - alpha / 2.0, dof)

    ci = np.zeros((p, 2), dtype=np.float64)

    for i, p, var in zip(range(n), pars, np.diag(pcov)):
        sigma = var ** 0.5
        ci[i, :] = [p - sigma * tval, p + sigma * tval]

    return ci


def confidence_interval(res: OptimizeResult, **kwargs):
    """
    This function estimates the confidence interval for the optimized parameters
    from the fit.
    
    Parameters
    ----------
    res: OptimizeResult
        The optimized result from least_squares minimization
    **kwargs
        confidence: float
            The confidence level (default 0.95)
    Returns
    -------
    np.ndarray
        ci: The confidence interval
    """
    if not isinstance(res, OptimizeResult):
        raise ValueError('Argument \'res\' should be an instance of \'scipy.optimize.OptimizeResult\'')

    confidence = kwargs.get('confidence', 0.95)

    # The vector of residuals at the solution
    residuals = res.fun
    # The number of data points
    n = len(residuals)
    # The number of parameters
    p = len(res.x)
    # The degrees of freedom
    dfe = n - p
    # Get MSE. The degrees of freedom when J is full rank is v = n-p and n-rank(J) otherwise
    mse = (LA.norm(residuals)) ** 2 / dfe

    # Needs to estimate the jacobian at the predictor point!!!
    # ypred = func(x,res.x)
    # delta = np.zeros((len(ypred),p));
    # fdiffstep       = np.amax(np.spacing(res.x)**(1/3));
    # for i in range(p):
    #     change = np.zeros(p)
    #     if res.x[i] == 0:
    #         nb = np.sqrt(LA.norm(res.x))
    #         change[i] = fdiffstep * (nb + (nb == 0))
    #     else:
    #         change[i] = fdiffstep * res.x[i]
    #
    #     predplus    = func(x,res.x+change)
    #     delta[:,i]  = (predplus - ypred)/change[i]

    # Find R to get the variance
    _, R = LA.qr(res.jac)
    # Get the rank of jac
    Rinv = LA.pinv(R)

    v = np.sum(Rinv ** 2, axis=1) * mse
    alpha = 1.0 - confidence
    tval = t.ppf(1.0 - alpha / 2.0, dfe)
    delta = np.sqrt(v) * tval
    ci = np.zeros((p, 2), dtype=np.float64)

    for i, p, d in zip(range(n), res.x, delta):
        ci[i, :] = [p - d, p + d]

    return ci


def get_rsquared(x: np.ndarray, y: np.ndarray, popt: np.ndarray, func: Callable[[np.ndarray, np.ndarray], np.ndarray]):
    """
    This function estimates R^2 for the fitting

    Parameters
    ----------
    x : np.ndarray
        The experimetnal x points
    y : np.ndarray
        The experimental y points
    popt: np.ndarray
        The best fit parameters
    func: Callable[[np.ndarray, np.ndarray]
        The fitted function
    
    Returns
    -------
    float
        rsquared: The value of R^2
        
    Reference:
        http://bagrow.info/dsv/LEC10_notes_2014-02-13.html
    """

    # Get the sum of the residuals from the linear function
    slin = np.sum((y - func(x, *popt)) ** 2)
    # Get the sum of the residuals from the constant function
    scon = np.sum((y - y.mean()) ** 2)
    # Get r-squared
    return 1.0 - slin / scon


def predband(x: np.ndarray, xd: np.ndarray, yd: np.ndarray, p: np.ndarray,
             func: Callable[[np.ndarray, np.ndarray], np.ndarray], conf: float = 0.95):
    """
    This function estimates the prediction bands for the specified function without using the jacobian of the fit
    https://codereview.stackexchange.com/questions/84414/obtaining-prediction-bands-for-regression-model

    Parameters
    ----------
    x: np.ndarray
        The requested data points for the prediction bands
    xd: np.ndarray
        The experimental values for x
    yd: np.ndarray
        The experimental values for y
    p: np.ndarray
        The fitted parameters
    func: Callable[[np.ndarray, np.ndarray], np.ndarray]
        The optimized function
    conf: float
        The confidence level

    Returns
    -------
    np.ndarray:
        yp: The value of the function at the requested points (x)
    np.ndarray:
        lpb: The lower prediction band
    np.ndarray:
        upb: The upper prediction band
    """
    alpha = 1.0 - conf  # significance
    npoints = len(xd)  # data sample size
    var_n = len(p)  # number of parameters
    # Quantile of Student's t distribution for p=(1-alpha/2)
    from scipy.stats.distributions import t
    q = t.ppf(1.0 - alpha / 2.0, npoints - var_n)
    # Stdev of an individual measurement
    se = np.sqrt(1. / (npoints - var_n) * np.sum((yd - func(xd, *p)) ** 2))
    # Auxiliary definitions
    sx = (x - xd.mean()) ** 2
    sxd = np.sum((xd - xd.mean()) ** 2)
    # Predicted values (best-fit model)
    yp = func(x, *p)
    # Prediction band
    dy = q * se * np.sqrt(1.0 + (1.0 / npoints) + (sx / sxd))
    # Upper & lower prediction bands.
    lpb, upb = yp - dy, yp + dy
    return yp, lpb, upb


def mean_squared_error(yd: np.ndarray, ym: np.ndarray):
    """
    This function estimates the mean squared error of a fitting

    Parameters
    ----------
    yd: np.ndarray
        The observed data points
    ym: np.ndarray
        The datapoints from the model

    Returns
    -------
    double
        mse: The mean squared error
    """
    if len(yd) != len(ym):
        raise ValueError('The length of the observations should be the same ' +
                         'as the length of the predictions.')
    if len(yd) <= 1:
        raise ValueError('Too few datapoints')
    n = len(yd)
    mse = np.sum((yd - ym) ** 2) / n
    return mse


def predint(x: np.ndarray, xd: np.ndarray, yd: np.ndarray, func: Callable[[np.ndarray, np.ndarray], np.ndarray],
            res: OptimizeResult, **kwargs):
    """
    This function estimates the prediction bands for the fit
    (see: https://www.mathworks.com/help/curvefit/confidence-and-prediction-bounds.html)
    Parameters 
    ----------
    x: np.ndarray
        The requested x points for the bands
    xd: np.ndarray
        The x datapoints
    yd: np.ndarray
        The y datapoints
    func: Callable[[np.ndarray, np.ndarray]
        The fitted function
    res: OptimizeResult
        The optimzied result from least_squares minimization
    **kwargs
        confidence: float
            The confidence level (default 0.95)
        simulateneous: bool
            True if the bound type is simultaneous, false otherwise
        mode: [functional, observation]
            Default observation        
    """

    if len(yd) != len(xd):
        raise ValueError('The length of the observations should be the same ' +
                         'as the length of the predictions.')
    if len(yd) <= 1:
        raise ValueError('Too few datapoints')
    from scipy.optimize import optimize

    if not isinstance(res, optimize.OptimizeResult):
        raise ValueError('Argument \'res\' should be an instance of \'scipy.optimize.OptimizeResult\'')

    simultaneous = kwargs.get('simultaneous', True)
    mode = kwargs.get('mode', 'observation')
    confidence = kwargs.get('confidence', 0.95)

    p = len(res.x)

    # Needs to estimate the jacobian at the predictor point!!!
    ypred = func(x, res.x)
    if callable(res.jac):
        delta = res.jac(x)
    else:
        delta = np.zeros((len(ypred), p))
        fdiffstep = np.spacing(np.abs(res.x)) ** (1 / 3)
        #    print('diff_step = {0}'.format(fdiffstep))
        #    print('popt = {0}'.format(res.x))
        for i in range(p):
            change = np.zeros(p)
            if res.x[i] == 0:
                nb = np.sqrt(LA.norm(res.x))
                change[i] = fdiffstep[i] * (nb + (nb == 0))
            else:
                change[i] = fdiffstep[i] * res.x[i]

            predplus = func(x, res.x + change)
            delta[:, i] = (predplus - ypred) / change[i]
    #    print('delta = {0}'.format(delta))

    # Find R to get the variance
    _, R = LA.qr(res.jac)
    # Get the rank of jac
    rankJ = res.jac.shape[1]
    Rinv = LA.pinv(R)
    pinvJTJ = np.dot(Rinv, Rinv.T)

    # The residual
    resid = res.fun
    n = len(resid)
    # Get MSE. The degrees of freedom when J is full rank is v = n-p and n-rank(J) otherwise
    mse = (LA.norm(resid)) ** 2 / (n - rankJ)
    # Calculate Sigma if usingJ 
    Sigma = mse * pinvJTJ

    # Compute varpred
    varpred = np.sum(np.dot(delta, Sigma) * delta, axis=1)
    #    print('varpred = {0}, len: '.format(varpred,len(varpred)))
    alpha = 1.0 - confidence
    if mode == 'observation':
        # Assume a constant variance model if errorModelInfo and weights are 
        # not supplied.
        errorVar = mse * np.ones(delta.shape[0])
        #        print('errorVar = {0}, len: '.format(errorVar,len(errorVar)))
        varpred += errorVar
    # The significance
    if simultaneous:
        from scipy.stats.distributions import f
        sch = [rankJ + 1]
        crit = f.ppf(1.0 - alpha, sch, n - rankJ)
    else:
        from scipy.stats.distributions import t
        crit = t.ppf(1.0 - alpha / 2.0, n - rankJ)

    delta = np.sqrt(varpred) * crit

    lpb = ypred - delta
    upb = ypred + delta

    return ypred, lpb, upb

def predint_multi(x: np.ndarray, xd: np.ndarray, yd: np.ndarray, 
                  func: Callable[[np.ndarray, np.ndarray], np.ndarray],
                  res: OptimizeResult, **kwargs):
    """
    This function estimates the prediction bands for the fit
    (see: https://www.mathworks.com/help/curvefit/confidence-and-prediction-bounds.html)
    Parameters 
    ----------
    x: np.ndarray
        The requested x points for the bands
    xd: np.ndarray
        The x datapoints
    yd: np.ndarray
        The y datapoints
    func: Callable[[np.ndarray, np.ndarray]
        The fitted function
    res: OptimizeResult
        The optimzied result from least_squares minimization
    **kwargs
        confidence: float
            The confidence level (default 0.95)
        simulateneous: bool
            True if the bound type is simultaneous, false otherwise
        mode: [functional, observation]
            Default observation        
    """

    if len(yd) != len(xd):
        raise ValueError('The length of the observations should be the same ' +
                         'as the length of the predictions.')
    if len(yd) <= 1:
        raise ValueError('Too few datapoints')
    from scipy.optimize import optimize

    if not isinstance(res, optimize.OptimizeResult):
        raise ValueError('Argument \'res\' should be an instance of \'scipy.optimize.OptimizeResult\'')

    simultaneous = kwargs.get('simultaneous', True)
    mode = kwargs.get('mode', 'observation')
    confidence = kwargs.get('confidence', 0.95)

    p = len(res.x)
    
    

    # Needs to estimate the jacobian at the predictor point!!!
    ypred = func(x, res.x)
    cols = ypred.shape[1]
    rows = len(x)
    if callable(res.jac):
        delta = res.jac(x)
    else:
        delta = np.zeros((cols*rows, p))
        fdiffstep = np.spacing(np.abs(res.x)) ** (1 / 3)
        #    print('diff_step = {0}'.format(fdiffstep))
        #    print('popt = {0}'.format(res.x))
        for i in range(p):
            change = np.zeros(p)
            if res.x[i] == 0:
                nb = np.sqrt(LA.norm(res.x))
                change[i] = fdiffstep[i] * (nb + (nb == 0))
            else:
                change[i] = fdiffstep[i] * res.x[i]

            predplus = func(x, res.x + change)
            for j in range(cols):
                for k in range(rows):
                    n = int(j*rows + k) 
                    with warnings.catch_warnings():
                        warnings.filterwarnings('error')
                        try:
                            delta[n, i] = (predplus[k,j] - ypred[k,j]) / change[i]
                        except Warning as e:
                            print(e)
                            print('change:')
                            print(change)
                            print('res.x:')
                            print(res.x)
    #    print('delta = {0}'.format(delta))

    # Find R to get the variance
    _, R = LA.qr(res.jac)
    # Get the rank of jac
    rankJ = res.jac.shape[1]
    Rinv = LA.pinv(R)
    pinvJTJ = np.dot(Rinv, Rinv.T)

    # The residual
    resid = res.fun
    n = len(resid)
    # Get MSE. The degrees of freedom when J is full rank is v = n-p and n-rank(J) otherwise
    mse = (LA.norm(resid)) ** 2 / (n - rankJ)
    # Calculate Sigma if usingJ 
    Sigma = mse * pinvJTJ

    # Compute varpred
    varpred = np.sum(np.dot(delta, Sigma) * delta, axis=1)
    #    print('varpred = {0}, len: '.format(varpred,len(varpred)))
    alpha = 1.0 - confidence
    if mode == 'observation':
        # Assume a constant variance model if errorModelInfo and weights are 
        # not supplied.
        errorVar = mse * np.ones(delta.shape[0])
        #        print('errorVar = {0}, len: '.format(errorVar,len(errorVar)))
        varpred += errorVar
    # The significance
    if simultaneous:
        from scipy.stats.distributions import f
        sch = [rankJ + 1]
        crit = f.ppf(1.0 - alpha, sch, n - rankJ)
    else:
        from scipy.stats.distributions import t
        crit = t.ppf(1.0 - alpha / 2.0, n - rankJ)

    delta = np.sqrt(varpred) * crit
        
    lpb = np.empty((rows,cols), dtype=np.float)
    upb = np.empty((rows,cols), dtype=np.float)
        
    for j in range(cols):
        for k in range(rows):
            n = int(j*rows + k) 
            lpb[k,j] = ypred[k,j] - delta[n]
            upb[k,j] = ypred[k,j] + delta[n]

#    lpb = ypred - delta
#    upb = ypred + delta

    return ypred, lpb, upb

# References:
# - Statistics in Geography by David Ebdon (ISBN: 978-0631136880)
# - Reliability Engineering Resource Website:
# - http://www.weibull.com/DOEWeb/confidence_intervals_in_simple_linear_regression.htm
# - http://reliawiki.com/index.php/Simple_Linear_Regression_Analysis#Confidence_Intervals_in_Simple_Linear_Regression
# - University of Glascow, Department of Statistics:
# - http://www.stats.gla.ac.uk/steps/glossary/confidence_intervals.html#conflim
