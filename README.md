Digital Signal Processing Building Blocks
---

![Language](https://img.shields.io/badge/Language-C++17-blue)
[![License](https://img.shields.io/badge/License-MIT-blue)](#license)
![Build](https://github.com/petiaccja/DSPBB/workflows/Build/badge.svg)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=petiaccja_DSPBB&metric=alert_status)](https://sonarcloud.io/dashboard?id=petiaccja_DSPBB)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=petiaccja_DSPBB&metric=coverage)](https://sonarcloud.io/dashboard?id=petiaccja_DSPBB)

DSPBB is a modern C++17/20 library for digital signal processing. DSPBB is not geared towards any specific signal processing domain like audio processing or seismology, it only aims to provide the building blocks for such applications:
- Classes for signals and systems
- Waveform generation
- Arithmetic & math functions
- Statistical functions
- Filtering:
  - FFT
  - IIR filter design
  - FIR flter design
  - Resampling
  - Windowing

For a more detailed list of features, [scroll down](#features).

### Development status

The library now implements most of the features that I was planning and is largely covered by unit tests. It's currently in a consolidation phase where I'm adding smaller features, ironing out problems, and improving ergonomy as I start using it in larger projects. As such, I'm expecting moderate, but breaking changes to the interface.

<a name="features"></a>
## Features

✔️ Supported<br>
❔️ Partially suppported / unknown<br>
❌️ Not supported but planned

### Design

- No dependencies
  - ❌️ XSimd is optional (needed for vectorization)
  - ❌️ Eigen is optional (needed for least-squares FIR design)
  - ✔️ PocketFFT included
- Vectorization
  - ✔️ Most code is vectorized
- Embedded-friendly
  - ❔️ Avoid memory allocation (partial)
  - ❌️ Allocator awareness
  - ✔️ No recursion (3rd parties not verified)

### Functionality

- Primitives:
  - ✔️ Signal
  - ✔️ SignalView
  - ✔️ Arithmetic operators
- Generators
  - Waveform
    - ✔️ Sine
    - ✔️ Square
    - ✔️ Sawtooth (fw, bw, tri, any)
    - ✔️ PWM
    - ✔️ Chirp/Sweep (for all the above types)
  - Space
    - ✔️ Linspace
    - ✔️ Logspace
- Mathematical functions
  - ✔️ Complex (Abs, Arg, Real, Imag, Conj)
  - ✔️ Trigonometric (Sin, Cos, Tan x inverse)
  - ✔️ Hyperbolic (Sinh, Cosh, Tanh x inverse)
  - ✔️ Exponential (Exp, Log, Log2, Log10)
  - ✔️ Polynomial (Pow, Sqrt, Cbrt)
  - ✔️ Erf and gamma
- Vector math basics
  - ✔️ Dot product
  - ✔️ Norm (sqrt(SumSquare))
- Statistics
  - ✔️ Sum
  - ✔️ Mean
  - ✔️ SumSquare
  - ✔️ MeanSquare
  - ✔️ RootMeanSquare
  - ✔️ Min/Max
  - ✔️ Central moments
  - ✔️ Standardized moments
  - ✔️ Standard deviation (popultion & corrected)
  - ✔️ Variance (popultion & corrected)
  - ✔️ Skewness (popultion & corrected)
  - ✔️ Kurtosis (popultion & corrected)
  - ✔️ Covariance (popultion & corrected)
  - ✔️ Correlation
- Filtering
  - Convolution
    - ✔️ Regular
    - ✔️ Overlap-add
  - FFT
    - ✔️ R->C, C->C, C->C, C->R
    - ✔️ FFT shift
    - ✔️ Bin <-> Frequency conversions
  - FIR filtering
    - Methods:
      - ✔️ Window method
      - ✔️ Least-squares method
      - ❌️ Parks-McClellan method
    - Types:
      - ✔️ Low-pass
      - ✔️ High-pass
      - ✔️ Band-pass
      - ✔️ Band-stop
      - ✔️ Arbitrary response
      - ✔️ Hilbert
    - Realizations:
      - ✔️ Convolution
      - ✔️ Overlap-add
  - IIR filtering
    - Methods:
      - ✔️ Butterworth
      - ✔️ Chebyshev I
      - ✔️ Chebyshev II
      - ✔️ Elliptic
    - Types:
      - ✔️ Low-pass
      - ✔️ High-pass
      - ✔️ Band-pass
      - ✔️ Band-stop
      - ❌️ Notch
    - Realizations:
      - ✔️ Direct form I.
      - ✔️ Direct form II.
      - ✔️ Cascaded biquad
  - Filter response analysis
    - ✔️ Compute amplitude & phase response
    - ✔️ Classify amplitude response: LP/HP/BP/BS
    - Measure amplitude parameters for LP/HP/BP/BS
      - ✔️ Transition edges
      - ✔️ Stopband attenuation
      - ✔️ Passband ripple
  - Polyphase FIR decomposition
  - Resampling
    - ✔️ Decimation (every n-th)
    - ✔️ Expansion (zero-fill)
    - ✔️ Interpolation (polyphase)
    - ✔️ Arbitrary resampling (polyphase)
  - Windowing
    - Derived properties
      - ✔️ Gain
      - ✔️ Energy
    - Functions
      - ✔️ Rectangular
      - ✔️ Triangular
      - ✔️ Hamming
      - ✔️ Blackman
      - ✔️ Blackman-Harris
      - ✔️ Flat top
      - ✔️ Kaiser
      - ✔️ Gaussian
      - ✔️ Dolph-Chebyshev
      - ✔️ Lanczos
  - Wavelets
    - ❌️ DWT
    - ❌️ CWT
    - ❌️ Wavelet design
- Literals
  - ❌️ Hz, kHz, MHz, GHz, THz
  - ❌️ dB, dB10

<a name="license"></a>
## License

DSPBB is distributed under the MIT license, therefore can be used in commercial and non-commercial projects alike with very few restrictions. The dependencies of DSPBB are also distributed under premissive MIT or BSD licenses.