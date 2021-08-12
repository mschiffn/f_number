# Frequency-dependent F-number for coherent plane-wave compounding

[![GitHub](license-image)](license-url)
[![GitHub](downloads-image)](downloads-url)
![GitHub repo size](https://img.shields.io/github/repo-size/mschiffn/comp_ui_toolbox)
![GitHub All Releases](https://img.shields.io/github/downloads/mschiffn/comp_ui_toolbox/total)

[license-image]: https://img.shields.io/github/license/mschiffn/comp_ui_toolbox
[license-url]: https://github.com/mschiffn/comp_ui_toolbox/COPYING
[downloads-image]: https://img.shields.io/github/downloads/mschiffn/comp_ui_toolbox/total
[downloads-url]: https://npmjs.org/package/ieee754

A simple [MATLAB](mathworks-url) implementation of
the frequency-dependent F-number for
coherent plane-wave compounding (CPWC).

![CIRS040](./figures/f_number_effect.png)

[mathworks-url]: https://mathworks.com/products/matlab.html

## Motivation

The F-number is
an important parameter in
all image formation methods using
the delay-and-sum (DAS) algorithm, such as

- coherent plane-wave compounding, or
- synthetic aperture imaging.

The F-number significantly reduces
image artifacts.

The F-number, for a uniform linear array, equals
the quotient of
the focal length and
the width of
the receive subaperture.
The width of the receive subaperture depends on
the focal length.
The parameter results in
a dynamic receive aperture.

Two methods to compute
the F-number have been proposed.

1. Directivity of the array elements:
This method attributes
the image artifacts to
noise.
The element directivity attenuates
the recorded signals and reduces
the signal-to-noise ratio.

2. Undersampling: This method attributes
the image artifacts to
the grating lobes.
The F-number enforces
a specific grating lobe-to-main lobe ratio.

Both methods, although they yield
similar F-numbers (1 <= F <= 2), are
mutually contradictory.

## What Does the Frequency-Dependent F-number Accomplish?

The proposed Fourier-domain beamforming algorithm accounts for
the strong frequency dependence of
the F-number.
It varies the width of
the receive subaperture with
the frequency and, in contrast to
a fixed F-number, includes
additional frequency components.
These improve both
the contrast and
the spatial resolution.

| Method            | Width of the receive subaperture  | Spatial resolution | Grating lobe suppression |
| ----------------- | --------------------------------- | ------------------ | ------------------------ |
| No F-number       | always full                       | optimal            | none                     |
| Fixed F-number    | position-dependent                | minimal            | exaggerated              |
| Proposed F-number | frequency- and position dependent | almost optimal     | almost optimal           |

## References :notebook:

1. M. F. Schiffner and G. Schmitz, "Frequency-Dependent F-Number Increases the Contrast and the Spatial Resolution in Fast Pulse-Echo Ultrasound Imaging", 2021 IEEE Int. Ultrasonics Symp. (IUS), accepted
