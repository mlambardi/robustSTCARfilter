# Robust STCAR Filter

This is an [Rcpp code](https://www.rcpp.org/) that can be loaded in the `R` console via the command
```
Rcpp::sourceCpp("stcar.cpp")
```

The code allows to estimate a Spatio-Temporal Conditional Auto-Regressive (STCAR) model, with slowly time-varying parameters, via an adaptive filter that is robust to outliers.

This code serves as a supplementary material for the reviewers of a manuscript submitted to Statistical Models & Applications.

Moreover, it serves as a (late) supplement to two conference papers:

[1] M Lambardi di San Miniato, R Bellio, L Grassetti, P Vidoni (2022) *Adaptive filters for time-varying correlation parameters*. In: Proceedings of the 51st Scientific Meeting of the Italian Statistical Society, Pearson. pp 1389-1394

[2] M Lambardi di San Miniato, R Bellio, L Grassetti, P Vidoni (2022) *Robust regression and adaptive filtering*. In: Proceedings of the 36th International Workshop on Statistical Modelling, EUT Edizioni Università di Trieste. pp 211-216

* The parameters are assumed to vary over time and their estimates are only slightly updated in the event of new data input, to cope with regime shifts and the obsolescence of prediction rules. So an **adaptive filter** is used to re-estimate the parameters without parsing again all past data. See more in [1].

* The estimation can be made robust, or, less sensitive to outliers, by resorting to **median** regression. In this way, $\hat{Y}$ will over-predict $Y$ roughly half of the times. See more in [1] and [2].

## Model

The function allows to track the time-varying parameters of a Spatio-Temporal Conditional Auto-Regressive model. This model can be used to predict a spatio-temporally referenced real-valued response variable $Y_{st}$ based on a covariate vector $x_{st}$, exploiting the autoregressive correlation of the response over time and a conditional autoregressive correlation structure in space.

One can define, in sequence, a detrended response

$$U_{st} = Y_{st} - \beta^\top x_{st'} ,$$

with $\beta$ a regression coefficient vector, $t'=t-\Delta$ and $\Delta>0$ a time horizon for prediction. Then, a temporal correlation structure is removed, as

$$V_{st} = U_{st} - \phi U_{st'} ,$$

with $\phi$ a temporal correlation parameter. One can additionally remove seasonal-type correlations as

$$V_{st} \leftarrow V_{st} - \phi_S V_{s(t-\Delta_S)} ,$$

with $\Delta_S$ the period of seasonality.

Then, one can define a set of spatial weights $w_{s,s'}$ that help to remove correlation in space as

$$\epsilon_{st} = V_{st} - \rho \bar{V}_{st'} ,$$

with

$$\bar{z}_{st} = \sum_{s'\neq s} w_{s,s'} z_{s't}$$

and $\rho$ a spatial correlation parameter.

Ideally, with suitable values of $\beta$, $\phi$ and $\rho$, to be estimated, $\epsilon_{st}$ should be purely white noise, so $Y$ can be predicted as

$$\hat{Y} = \rho \bar{Y}_{st'} + \phi Y_{st'} + \beta^\top x_{st'} - \rho \phi \bar{Y}_{st''} - \rho \beta^\top \bar{x}_{st''} - \phi \beta^\top x_{st''} + \rho \phi \beta^\top \bar{x}_{st'''} ,$$

so the model predicts $Y$ using lags of $x$ and $Y$ itself, in a structured fashion that makes $\beta$, $\rho$ and $\phi$ interpretable as correlation parameters *net of* all the other correlation structures.
