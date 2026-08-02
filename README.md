# Global forest fragmentation and extinction risk for forest-dependent raptors

This repository contains the data and R code for the manuscript:

**Li, H., Wittig, C., Buij, R., Buechley, E.R., Xie, Z., Torres-Romero, E.J., Surasinghe, T.D., Foysal, M., Shandilya, A.S., Martens, P., & O'Bryan, C.J. (2025).** *Global forest fragmentation signals declining habitat quality and increased extinction risk for forest-dependent raptors.*

## Key findings

- Forest specialists experienced greater fragmentation increases than generalists, particularly among small-bodied species.
- Fragmentation change showed a hump-shaped relationship with range size but was unrelated to IUCN threat status.
- Landscape configuration explained fragmentation change more strongly than forest-loss trends.
- Nearly half of the highest-risk forest specialists are currently classified as Least Concern, revealing a systematic blind spot in current threat assessments.

## R packages required

```r
library(terra)
library(dplyr)
library(stringr)
library(ggplot2)
library(mgcv)
library(multcompView)
library(agricolae)
library(showtext)
library(scales)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(ggcorrplot)
library(randomForest)
library(caret)
library(sf)
library(exactextractr)
