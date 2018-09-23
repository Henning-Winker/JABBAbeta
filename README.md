## JABBAbeta
Development repository for JABBA (https://github.com/jabbamodel)

<B> New JABBA beta version [`JABBAv1.2beta.R`](https://github.com/Henning-Winker/JABBAbeta/blob/master/SWO_SA_prime_v1.2.R) is now available!</b>

This new beta version has been developed and tested during stock assessments of:
+ [ICCAT Atlantic blue marlin (BUM)](https://www.iccat.int/Documents/Meetings/Docs/2018/REPORTS/2018_BUM_SA_ENG.pdf)
+ [ICCAT Atlantic bigeye tuna (BET)](https://www.iccat.int/Documents/Meetings/Docs/2018/REPORTS/2018_BET_SA_ENG.pdf)
+ [IOTC Indian Ocean striped marlin (MLS)](http://www.iotc.org/documents/WPB/16/16-MLS_JABBA)
+ [IOTC Indian Ocean black marlin (BLM)](http://www.iotc.org/documents/WPB/16/15-BLM_JABBA)
+ [NOAA Hawaii Kona crab benchmark assessment (KONA)](https://www.fisheries.noaa.gov/pacific-islands/population-assessments/western-pacific-stock-assessment-review#2018-kona-crab-in-the-main-hawaiian-islands)

[New Features](https://github.com/Henning-Winker/JABBAbeta/tree/master/V1.2_NewFeatures) include:
+ Plotting code is outsouced in [`JABBA_plots_v1.2.R`](https://github.com/Henning-Winker/JABBAbeta/blob/master/JABBA_plots_v1.2beta.R) to facilitate debugging
+ Settings.txt saved for reference in Input folder
+ Preliminary estimate shape m option with informative 
+ Catch.CV option: Allows addimitting uncertainty about the catch
+ CatchOnly option: JABBA run with catch and priors, but without fitting any abundance indices
+ Lower and upper values of P_bound, K_bound, q_bound can be set manually to enforce "soft boundaries (CV=0.1)     
+ Option to manually set starting values for r, q and K

See examples [`SWO_SA_NewFeatures_v1.2.R`](https://github.com/Henning-Winker/JABBAbeta/blob/master/V1.2_NewFeatures/SWO_SA_NewFeatures_v1.2.R)


## JABBA: Just Another Bayesian Biomass Assessment
The materials in this repository present the stock assessment tool ‘Just Another Bayesian Biomass Assessment’ JABBA. The motivation for developing JABBA was to provide a user-friendly R to JAGS (Plummer) interface for fitting generalized Bayesian State-Space SPMs with the aim to generate reproducible stock status estimates and diagnostics. Building on recent advances in optimizing the fitting procedures through the development of Bayesian state-space modelling approaches, JABBA originates from a continuous development process of a Bayesian State-Space SPM tool that has been applied and tested in many assessments across oceans. JABBA was conceived in the Summer of 2015 as a collaboration between the South Africa Department of Agriculture, Forestry and Fisheries and the Pacific Islands Fisheries Science Center (NOAA) in Honolulu, HI USA. The goal was to provide a bridge between age-structured and biomass dynamic models, which are still widely used. JABBA runs quickly and by default generates many useful plots and diagnosic tools for stock assessments.

Inbuilt JABBA features include:

+ Integrated state-space tool for averaging multiple CPUE series (+SE) for optional use in assessments
+ Automatic fitting of multiple CPUE time series and associated standard errors
+ Fox, Schaefer or Pella Tomlinson production function (optional as input Bmsy/K)
+ Kobe-type biplot plotting functions 
+ Forecasting for alternative TACs 
+ Residual and MCMC diagnostics 
+ Estimating or fixing the process variance
+ Optional estimation additional observation variance for individual or grouped CPUE time series
+ Easy implementation of time-block changes in selectivity

**Reference**

[Winker, H., Carvalho, F., Kapur, M. (2018) <U>JABBA: Just Another Bayesian Biomass Assessment.</U> *Fisheries Research* **204**: 275-288.](https://www.sciencedirect.com/science/article/pii/S0165783618300845)   

<B>A self-contained R package of JABBA is forthcoming.</b>
