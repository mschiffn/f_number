# Frequency-Dependent F-Number for Coherent Plane-Wave Compounding

[![GitHub](license-image)](license-url)
[![GitHub](downloads-image)](downloads-url)
![GitHub repo size](https://img.shields.io/github/repo-size/mschiffn/f_number)
![GitHub All Releases](https://img.shields.io/github/downloads/mschiffn/f_number/total)

[license-image]: https://img.shields.io/github/license/mschiffn/f_number
[license-url]: https://github.com/mschiffn/f_number/COPYING
[downloads-image]: https://img.shields.io/github/downloads/mschiffn/f_number/total
[downloads-url]: https://npmjs.org/package/ieee754

Simple [MATLAB](mathworks-url) implementation of
the frequency-dependent F-number [[1]](#proc:SchiffnerIUS2021) for
coherent plane-wave compounding (CPWC).

![CIRS040](./figures/f_number_effect.png)

[mathworks-url]: https://mathworks.com/products/matlab.html

## What is the F-Number?

The F-number significantly reduces
image artifacts in
all image formation methods using
the delay-and-sum (DAS) algorithm, such as

- coherent plane-wave compounding [[2]](#article:MontaldoITUFFC2009), or
- synthetic aperture imaging [[3]](#article:JensenUlt2006).

The F-number, for a uniform linear array, equals
the quotient of
the focal length and
the width of
the receive subaperture.

![F-number](./figures/f_number_definition.png)

The usage of
a fixed F-number, thus, results in
a dynamic receive subaperture whose
width increases with
the focal length.

## Established Methods to Compute the Optimal F-Number are Mutually Contradictory

Established methods to compute
the optimal F-number attribute
the image artifacts to
two different phenomena:

1. Noise [[4]](#article:PerrotUlt2021), [[5]](#book:Szabo2013), [[2]](#article:MontaldoITUFFC2009):
The directivity of
the array elements attenuates
the echoes and reduces
the signal-to-noise ratio of
the recorded signals.

2. Grating lobes [[6]](#article:DelannoyJAP1979), [[7]](#article:BruneelJAP1978):
The width of
the receive subaperture determines
the grating lobe-to-main lobe ratio.

Both methods, although yielding
similar F-numbers (0.5 <= F <= 2), are
mutually contradictory.

Wide array elements, for example, exhibit
an increased directivity.
The first method suggests
the usage of
narrow receive subapertures or, equivalently,
large F-numbers for
such elements to avoid
noise problems.
The second method, in contrast, permits
wide receive subapertures or, equivalently,
small F-numbers for
such elements because they attenuate
the grating lobes.

Both methods, furthermore, yield
F-numbers with
a strong frequency dependence.

## What Does the Proposed F-Number Accomplish?

The proposed F-number not only eliminates
image artifacts but also maintains
the spatial resolution of
the full aperture [[1]](#proc:SchiffnerIUS2021).
This F-number, in particular, prevents
the first-order grating lobes from insonifying
reflective image structures by imposing
a lower bound on
the angular distance of
these lobes.

## How Does the Implementation Work?

The proposed Fourier-domain beamforming algorithm accounts for
the strong frequency dependence of
this F-number.
The algorithm not only varies
the width of
the receive subaperture with
the voxel position but also with
the frequency.
This additional frequency dependence, in contrast to
a fixed F-number, includes
additional frequency components to improve both
the contrast and
the spatial resolution.

| Method            | Width of the receive subaperture  | Spatial resolution | Grating lobe suppression |
| ----------------- | --------------------------------- | ------------------ | ------------------------ |
| No F-number       | always full                       | optimal            | none                     |
| Fixed F-number    | position-dependent                | minimal            | exaggerated              |
| Proposed F-number | frequency- and position-dependent | almost optimal     | optimal                  |

## Getting Started

1. Clone the repository or download the release to your local hard drive.

```
git clone https://github.com/mschiffn/f_number
```

2. Add the repository to your MATLAB path using .

```matlab
addpath( genpath( './f_number' ) )
```

## Folder Structure

The repository has the following structure:

    .
    ├── +f_numbers      # classes for various types of F-numbers (e.g., fixed, directivity-derived, proposed)
    ├── +windows        # classes for various window functions (e.g., boxcar, Hann, Tukey)
    ├── das_pw.m        # main function
    ├── LICENSE         # license file
    └── README.md       # this readme

The packages +f_numbers and +windows contain an exemplary class hierarchy to manage various types of F-numbers and window functions.
The proposed F-number can be instantiated by

```matlab
F_number_rx = f_numbers.grating.angle_lb( 45, 3 )
```

The directivity-derived F-number is
```matlab
F_number_rx = f_numbers.directivity.perrot( 0.918 )
```

## Image Formation

Use the function das_pw to form images.

In MATLAB type
```matlab
help das_pw
```

to obtain an explanation of the input and output arguments.


```matlab
[ image, F_number_values ] = das_pw( positions_x, positions_z, data_RF, f_s, e_theta, element_width, element_pitch, ( 1 - N_elements ) / 2, [ f_lb, f_ub ], c_ref, N_samples_shift, window, F_number);
```

## References :notebook:

1. <a name="proc:SchiffnerIUS2021"></a>
M. F. Schiffner and G. Schmitz, "Frequency-Dependent F-Number Increases the Contrast and the Spatial Resolution in Fast Pulse-Echo Ultrasound Imaging", 2021 IEEE Int. Ultrasonics Symp. (IUS), accepted

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
