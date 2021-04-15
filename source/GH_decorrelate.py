import numpy as np
import statsmodels.api as sm
from scipy.linalg import toeplitz

### Erasing noise and decorrelate
## Generalized Least Squares


data = sm.datasets.longley.load(as_pandas=False)
data.exog = sm.add_constant(data.exog)
print(data.exog[5])
ols_resid = sm.OLS(data.endog, data.exog).fit().resid

resid_fit = sm.OLS(ols_resid[1:], sm.add_constant(ols_resid[:-1])).fit()
print(resid_fit.tvalues[1])
print(resid_fit.pvalues[1])

rho = resid_fit.params[1]

order = toeplitz(range(len(ols_resid)))

sigma = rho**order
gls_model = sm.GLS(data.endog, data.exog, sigma=sigma)
gls_results = gls_model.fit()
print(gls_results.summary())

