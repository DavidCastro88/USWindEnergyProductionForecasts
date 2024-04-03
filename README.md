![image](https://github.com/DavidCastro88/USWindEnergyProductionForecasts/assets/91480088/45833281-2829-4ce3-9e89-2de30310df30)# USWindEnergyProductionForecasts


![image](https://github.com/DavidCastro88/USWindEnergyProductionForecasts/assets/91480088/023a84da-7899-4141-8557-b7f31a0864a8)

The purpose is to analyze the time series of wind energy production in the United States for the available periods ranging from 2001 to 2023. With this analysis, different models for time series will be proposed, first of all, that they are valid models and see your ability to forecast.

**Potential Uses**:

- Conducting time series analysis to forecast wind energy production and capacity factors
- Performing exploratory data analysis to identify trends and patterns in wind energy production
- Comparing wind energy production to other electricity generation sources to inform policy decisions
- Modeling wind energy production and capacity factors for forecasting and planning purposes
- Evaluating the impact of policy changes on wind energy production in the United States

**Data set Descriptiom**:
This dataset, provided by the U.S. Energy Information Administration (EIA) in the Electric Power Monthly report, contains monthly data on wind energy production and other renewables in the United States.

Here is some other informations about the variables available :

Time Range: January 2001 to the latest month available
Geographic Coverage: United States
Granularity: Monthly
Variables:
"date": Month and year
"wind_state_name" : wind power production for the current state
"other_state_name" : production for all other renewables sources for the current state.

In our case we are only going to be interested in the production of wind energy.

The series presents multiplicative components, so it is important to apply natural logarithm to its values to stabilize the variance, since many models are sensitive to the pattern of the series components.

![image](https://github.com/DavidCastro88/USWindEnergyProductionForecasts/assets/91480088/ce6b9621-47b8-453d-92ac-9a4c83f5b5e7)

In the series, trend components, seasonality, cycles and residuals were found.
![image](https://github.com/DavidCastro88/USWindEnergyProductionForecasts/assets/91480088/b8ff3dd0-a035-4afd-a907-d8a683e8a0c1)

The models that were implemented were ARIMA (model 1), LSTM (model 2) and Prophet (model 3). The validity of the model, its fit and its prognostic capacity were verified. (Python)

Models forecast results:
![image](https://github.com/DavidCastro88/USWindEnergyProductionForecasts/assets/91480088/b862d9c0-d209-4d72-81a1-2f905dbb31a1)

In R, different ARIMA models were tested only, where the parameters p,q,P,Q were obtained from the auto_arima() and armasubsets() functions.

1) SARIMA (1,1,1)(0,1,1)[12]
2) SARIMA (4,1,1)(2,0,0)[12] with drift
3) SARIMA (1,1,7)(1,1,2)[12] ma(q) with q only in 1,2,7.
4) SARIMA (1,1,2)(1,1,1)[12] 

Results adjustement of model 3 (best model):
![image](https://github.com/DavidCastro88/USWindEnergyProductionForecasts/assets/91480088/66f4284d-1a39-4244-b7d6-47b5155eeb8e)

Models forecast results:
![image](https://github.com/DavidCastro88/USWindEnergyProductionForecasts/assets/91480088/b398d67f-8c45-4133-ac76-bcd18778d435)

