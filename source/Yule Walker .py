#!/usr/bin/env python
# coding: utf-8

# In[3]:


import statsmodels.api as sm #this library contains yule-Walker Algorithm
X=[] #a sequence X using Yule-Walker equation
rho, sigma = sm.regression.yule_walker(X, order=4,method="mle")
"""
    Estimate AR(p) parameters from a sequence X using Yule-Walker equation.

    Unbiased or maximum-likelihood estimator (mle)

    Parameters
    ----------
    X : array-like
        1d array
    order : integer, optional
        The order of the autoregressive process.  Default is 1.
    method : string, optional
       Method can be "unbiased" or "mle" and this determines denominator in
       estimate of autocorrelation function (ACF) at lag k. If "mle", the
       denominator is n=X.shape[0], if "unbiased" the denominator is n-k.
       The default is unbiased.
    df : integer, optional
       Specifies the degrees of freedom. If `df` is supplied, then it is assumed
       the X has `df` degrees of freedom rather than `n`.  Default is None.
    inv : bool
        If inv is True the inverse of R is also returned.  Default is False.
    demean : bool
        True, the mean is subtracted from `X` before estimation.

    Returns
    -------
    rho
        The autoregressive coefficients
"""

