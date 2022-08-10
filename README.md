# Pneumococcal Conjugate Vaccine Effectiveness Prediction Function
***
##### Predicts pneumococcal conjugate vaccine (PCV) effectiveness for shared serotypes with current PCV using summary-level (mean and 95% CI) current PCV serotype specific effectiveness data, summary-level (geometric mean and 95% CI) current PCV serotype specific IgG data, and summary-level (geometric mean and 95% CI) new PCV serotype specific IgG data.

In the current state the function only uses summary-level data. A future iteration will have the capability to also use individual level immunogenicity data.

### Installation
***
devtools::install_github("jryman17/pneuve")

### How to Use the Function
***
This is an example which shows you how to predict a serotype's effectiveness

A few packages
```{r}
library(tidyverse)
```

Input number of simulations (nsim=1000) and reported summary-level data

Example summary-level data 

Serotype | Current PCV vaccine effectiveness (95% CI) | Placebo GMC (95% CI) | Current PCV GMC (95% CI) | New PCV GMC (95% CI)
:-------:| :----------------------------------------: | :-----------------:  | :----------------------: | :--------------:
   4     | 97% (65-100)                               | 0.03 (0.02-0.04)     | 1.61 (1.40-1.84)         |	1.55 (1.41-1.70)

Number of Subjects in each treatment arm

Placebo | Current PCV | New PCV 
:-----: | :---------: | :--------:
  189   | 123         | 235


* To run prediction enter in:
  + the number of simulations
  + keep seed=realization
  + observed vaccine effectiveness (starting with the mean, followed by the lower, and then the upper 95% CI)
  + the number of subjects in the placebo treated arm
  + the number of subjects in the current PCV treated arm
  + the number of subjects in the next PCV treated arm
  + the serotype being analyzed (in "")
  + the geometric mean concentration in the placebo treated arm
  + the upper 95% CI bound in the placebo treated arm
  + the geometric mean concentration in the current PCV treated arm
  + the upper 95% CI bound in the current PCV treated arm
  + the geometric mean concentration in the next PCV treated arm
  + the upper 95% CI bound in the new PCV treated arm

```{r}
ser4Prediction <- VEpred(nsim=1000,seed=realization,obs.VEs=c(97, 65, 100),nSubPlacebo=189,
                         nSubCurrentPCV=123,nSubNewPCV=235,serotype="pn4",
                         GMC.Placebo=0.03,upperCI.Placebo=0.04,GMC.CurrentPCV=1.61,
                         upperCI.CurrentPCV=1.84,GMC.NewPCV=1.55,upperCI.NewPCV=1.70)
```
Extract estimated protective threshold IgG concentration
```{r}
s4Pt <- ser4Prediction[[2]] 
```

Extract predicted effectiveness
```{r}
s4Eff <- ser4Prediction[[1]]
```

Extract simulation dataframe 
```{r}
df4 <- ser4Prediction[[3]]
```


For help please reach out to: josiah.ryman@merck.com

