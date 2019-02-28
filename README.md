# bayfoxm

A Bayesian, global planktic foraminifera core top calibration to sea-surface temperature (SST), for Octave/MATLAB.

**Please note that these scripts are currently under development.**


## What is bayfoxm?

bayfoxm is a suite of linear Bayesian calibration models for planktic core top foraminiferal δ18O (δ18Oc) and SSTs. These calibrations are especially useful because they capture the uncertainty in the relationship between modern SSTs and core top δ18Oc. This package is a companion to a paper currently under preparation for the journal "Paleoceanography and Paleoclimatology".

These calibration are also available in the [`bayfox` package](https://github.com/brews/bayfox) for Python, and the [`bayfoxr` package](https://github.com/brews/bayfoxr) for the R statistical environment.

## A quick example

If we have a few hypothetical SST observations (°C) and a seawater δ18O measurement (‰ VSMOW) we can then get an ensemble of predictions for foraminiferal δ18Oc (‰ VPDB):

```octave
sst = [15; 15.5; 16];
d18osw = 0.75;

d18oc_ensemble = predict_d18oc(sst, d18osw);
```

where `d18oc_ensemble` is a [n x m] matrix, an m-length δ18Oc posterior distribution 
for each of our n input SST observations in `sst`. We can also infer SSTs from foraminiferal δ18Oc, if we give information for a gaussian prior SST distribution:

```octave
d18oc = [0; 0.15; 0.2];  % (‰ VPDB)
d18osw = 0.75;  % (‰ VSMOW)

% SST gaussian prior parameters
prior_mean = 15.0;  % (°C)
prior_std = 10.0;  % (°C)

sst_ensemble = predict_seatemp(d18oc, d18osw, prior_mean, prior_std);
```

Much like before, `sst_ensemble` is a [n x m] matrix, an m-length posterior SST distribution for each of our n input δ18Oc observations in `d18oc`.

This is the simplest case, using the pooled annual calibration model. There are four calibration models available. These model are described extensively in our paper (see the *Citing bayfoxm in your research* section, below). Also see `help predict_d18oc` and `help predict_seatemp` for more details. In short, we can use the hierarchical annual model by giving a foraminifera species. For example

```octave
ens = predict_d18oc([15; 15.6; 16], 0.75, false, "G. bulloides")
```

and we can use the hierarchical seasonal calibration model by passing in `true`, instead of `false`, like so

```octave
ens = predict_d18oc([15; 15.6; 16], 0.75, true, "G. bulloides")
```

## Citing bayfoxm in your research

Please cite our work if you use bayfoxr in your research. We have a paper currently in preparation and I'll be sure to update this section with the citation as soon as the paper is out.

To cite the code repository directly use:

*Malevich, Steven B., 2019. bayfoxm. \<https://github.com/brews/bayfoxm \>.*

## Installation

These are just a few scripts so no heavy installation is needed. Just [download `bayfoxm`](https://github.com/brews/bayfox/archive/master.zip) and unzip the files in your working directory, or add the directory to your system path.

## Support

Documentation is included in the code and can be viewed in MATLAB/Octave with `help`. Please file issues and requests in the [bug tracker](https://github.com/brews/bayfoxm/issues).
