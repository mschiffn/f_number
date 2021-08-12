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

- coherent plane-wave compounding [[2]](#article:MontaldoITUFFC2009), or
- synthetic aperture imaging [[3]](#article:JensenUlt2006).

The F-number significantly reduces
image artifacts.

The F-number, for a uniform linear array, equals
the quotient of
the focal length and
the width of
the receive subaperture.
The usage of
a fixed F-number results in
a dynamic receive aperture whose
width depends on
the current focal length.

Methods to compute
the F-number rely on
two phenomena:

1. Directivity of the array elements [[4]](#article:PerrotUlt2021), [[5]](#book:Szabo2013), [[2]](#article:MontaldoITUFFC2009):
These methods attribute
the image artifacts to
noise.
The noise results from
the attenuation of
the recorded signals by
the element directivity.

2. Undersampling [[6]](#article:DelannoyJAP1979), [[7]](#article:BruneelJAP1978): This method attributes
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

## Getting Started

1. Clone or download the repository to your local hard drive.
2. Add the repository to your MATLAB path using addpath( genpath( 'your/folder/here' ) ).
3. Call the function das_pw to form images.

## Image Formation

Use the function das_pw to form images.

```matlab
[ image, F_number_values ] = das_pw( positions_x, positions_z, data_RF, f_s, e_theta, element_width, element_pitch, ( 1 - N_elements ) / 2, [ f_lb, f_ub ], c_ref, N_samples_shift, window, F_number);
```

The packages +f_numbers and +windows contain an exemplary class hierarchy to manage various types of F-numbers and window functions.
The proposed F-number can be instantiated by f_number = f_numbers.grating.angle_lb().
The directivity-derived F-number is f_number = f_numbers.directivity.perrot( , 3 ).

## References :notebook:

1. M. F. Schiffner and G. Schmitz, "Frequency-Dependent F-Number Increases the Contrast and the Spatial Resolution in Fast Pulse-Echo Ultrasound Imaging", 2021 IEEE Int. Ultrasonics Symp. (IUS), accepted

2. <a name="article:MontaldoITUFFC2009"></a>
G. Montaldo, M. Tanter, J. Bercoff, N. Benech, and M. Fink,
“Coherent plane-wave compounding for very high frame rate ultrasonography and transient elastography,"
IEEE Trans. Ultrason., Ferroelectr., Freq. Control, vol. 56, no. 3, pp. 489–506, Mar. 2009.
[![DOI:10.1109/TUFFC.2009.1067](https://img.shields.io/badge/DOI-10.1109%2FTUFFC.2009.1067-blue)](http://dx.doi.org/10.1109/TUFFC.2009.1067)

3. <a name="article:JensenUlt2006"></a>
J. A. Jensen, S. I. Nikolov, K. L. Gammelmark, and M. H. Pedersen,
“Synthetic aperture ultrasound imaging,” Ultrasonics, vol. 44, Supplement, e5–e15, Dec. 2006.
[![DOI:10.1016/j.ultras.2006.07.017](https://img.shields.io/badge/DOI-10.1016%2Fj.ultras.2006.07.017-blue)](http://dx.doi.org/10.1016/j.ultras.2006.07.017)

4. <a name="article:PerrotUlt2021"></a>
V. Perrot, M. Polichetti, F. Varray, and D. Garcia,
“So you think you can DAS? A viewpoint on delay-and-sum beamforming,”
Ultrasonics, vol. 111, p. 106 309, Mar. 2021.
[![DOI:10.1016/j.ultras.2020.106309](https://img.shields.io/badge/DOI-10.1016%2Fj.ultras.2020.106309-blue)](http://dx.doi.org/10.1016/j.ultras.2020.106309)

5. <a name="book:Szabo2013"></a>
T. L. Szabo, Diagnostic Ultrasound Imaging: Inside Out, 2nd. Elsevier Academic Press, Dec. 2013

6. <a name="article:DelannoyJAP1979"></a>
B. Delannoy, R. Torguet, C. Bruneel, E. Bridoux, J. M. Rouvaen, and H. Lasota,
“Acoustical image reconstruction in parallel-processing analog electronic systems,”
J. Appl. Phys., vol. 50, no. 5, pp. 3153–3159, May 1979.
[![DOI:10.1063/1.326397](https://img.shields.io/badge/DOI-10.1063%2F1.326397-blue)](http://dx.doi.org/10.1063/1.326397)

7. <a name="article:BruneelJAP1978"></a>
C. Bruneel, E. Bridoux, B. Delannoy, B. Nongaillard, J. M. Rouvaen, and R. Torguet,
“Effect of spatial sampling on an acoustical image reconstruction,”
J. Appl. Phys., vol. 49, no. 2, pp. 569–573, Feb. 1978.
[![DOI:10.1063/1.324680](https://img.shields.io/badge/DOI-10.1063%2F1.324680-blue)](http://dx.doi.org/10.1063/1.324680)
