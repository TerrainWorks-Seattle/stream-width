# Review of some recent attempts at estimating stream width and depth

## Summary
Several studies using the USGS [HYDRoSWOT](https://www.usgs.gov/data/usgs-hydroacoustic-dataset-support-surface-water-oceanographic-topography-satellite-mission) 
as a training dataset to develop ML models to predict stream width and depth. 
I don't see many significant differences between these different studies. 
All published datasets predicting steam width and depth for all NHDPlusV2 
reaches. 
Important predictors to consider: 
  - discharge (always by far the most important)
  - vegetation (measured as EVI)
  - elevation
  - aridity
  - geology (sand popped out a couple times)
Found pretty good success from random forest and 
xgboost. 

## Some ideas: 
- attempt to line up published predictions along NHDPlus reaches with 
  synthetic channels. Use these predictions as a covariate. We may 
  also want to include other covariates they excluded like slope metrics.
  We could still include spatial autocorrelation term if there is 
  still spatial autocorrelation (which I would guess there would be).
- This won't work in Alberta because NHD is only for US.
- Fit a RF model to predict class or width and use this as a predictor in
  the full model which will include spatial autocorrelation. 

## Notes
* Chang et al., 2024 (Water Resources Research) (([STREAM-geo](https://figshare.com/articles/dataset/Stream_Reach_Evaluation_and_Metrics_-_Geometry_STREAM-geo_/24463240/1)): 
  - Create independent models for median width and depth 
    using HYDRoSWOT as training dataset, testing 
    several different ML algorithms: 
      - Linear regression with stepwise variable selection
      - Random forest
      - XGBoost
      - Multilayer perceptron (neural net)
  - published dataset of predictions for all NHDPLUS
    reaches 
  - RF and XGBoost worked best. 
  - Important predictors: 
    - discharge
    - contributing area
    - elevation
    - EVI (only for width)
    - aridity (only for depth)
  - In linear regression, found negative 
    relationship between watershed vegetation and width 
    and depth. 
* Modaresi Rad et al., 2024 (Journal of Geophysical Research) ([predictions](https://www.hydroshare.org/resource/d147fcf554a54b2aaa4f146f85da0e03/))
  - Used ensemble ML to separately predict bankfull width, 
    in-channel width, bankfull depth, and in-channel depth
    HYDRoSWOT as training dataset
  - Most important predictors from SHAP values: 
    - discharge
    - length of upstream flowlines (width)
    - catchment area (width)
    - elevation (depth)
    - slope (depth)
* Zarrabi et al., 2025 (Water Resources Reseach) ([predictions](https://www.hydroshare.org/resource/63ae139ccd2445959470d0e5a2ebf6a5/))
  - create independent models for bankfull and mean width 
    and depth using HYDRoSWOT as training dataset
  - Most important predictors from MDI values: 
    - discharge
    - drainage area
    - stream order
    - EVI (width only)
    - aridity index (width only)
    - minimum elevation (depth only)
    - catchment avg percent sand (depth only)

