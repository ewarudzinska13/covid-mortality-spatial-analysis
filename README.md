# covid-mortality-spatial-analysis
Spatial econometric analysis of COVID-19 mortality across Polish counties using R

# Spatial Econometric Analysis of COVID-19 Mortality in Poland

Econometric analysis examining spatial patterns of COVID-19 mortality across Polish counties (powiaty) using advanced spatial regression models.

## Project Overview
This study investigates the spatial distribution and determinants of COVID-19 mortality in Poland at the county level, employing spatial econometric techniques to account for spatial dependencies.

## Methodology
- **Data preparation**: COVID-19 data merged with demographic indicators
- **Spatial weight matrix**: Contiguity-based (shared borders)
- **Spatial autocorrelation testing**: Moran's I statistic
- **Model comparison**: OLS, SLX, SAR, SEM, SDM, SDEM, SAC, GNS
- **Model selection**: Likelihood Ratio tests

## Key Findings
- Significant spatial autocorrelation detected in mortality patterns
- Spatial Error Model (SEM) identified as optimal specification
- Key determinants: infection rates per capita, quarantine coverage, population density, median age

## Tools & Technologies
- **Language**: R
- **Key packages**: `spdep`, `spatialreg`, `sf`, `sp`
- **Visualization**: Spatial maps with residual analysis

## Variables
- Dependent: Log(deaths)
- Independent: Log(infections per 10k), Log(quarantine), Log(population density), Log(median age)
