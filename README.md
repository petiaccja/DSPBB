Digital Signal Processing Building Blocks
---

![Build](https://github.com/petiaccja/DSPBB/workflows/Build/badge.svg)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=petiaccja_DSPBB&metric=alert_status)](https://sonarcloud.io/dashboard?id=petiaccja_DSPBB)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=petiaccja_DSPBB&metric=coverage)](https://sonarcloud.io/dashboard?id=petiaccja_DSPBB)

DSPBB is a modern C++17/20 library for digital signal processing. Unlike other libraries, DSPBB is not geared towards any specific field. Instead, it aims to provide the raw essentials only: signals, basic math, statistics, filtering (FFT, FIR, IIR, resampling, wavelets, etc.) and filter design. As such, you won't find a Meier crossfeed implementation in DSPBB, but you can use DSPBB to quickly implement it in your code or even to build a full-fledged library for music processing.

### Status of the library

While the library is already quite capable, I'm planning to change the interfaces frequently as well as to refactor and fine-tune the code, therefore it's not quite ready for production just yet. However, if you fancy the features listed below, stay tuned.

## Features

✔️ Supported
❔️ Partially suppported or unknown
❌️ Not supported but planned

### Design

- ❌️ No dependencies
  - ❌️ XSimd is optional
  - ❌️ Eigen is optional
  - ✔️ PocketFFT included
- ✔️ Vectorization
- ❔️ Embedded-friendly
  - ❔️ No malloc

### Functionality

- ✔️ Primitives:
  - ✔️ Signal
    - ✔️ Real
    - ✔️ Complex
    - ✔️ Arithmetic operators
  - ✔️ SignalView
    - ✔️ Arithmetic operators
- ✔️ Generators
  - ✔️ Wave
    - ✔️ Sine
    - ✔️ Square
    - ✔️ Sawtooth (fw, bw, tri, any)
    - ✔️ PWM
  - ✔️ Chirp/Sweep
  - ✔️ Space
    - ✔️ Linspace
    - ✔️ Logspace
- ✔️ Common math functions
  - ✔️ Complex (Abs, Arg, Real, Imag, Conj)
  - ✔️ Trigonometric (Sin, Cos, Tan x inverse)
  - ✔️ Hyperbolic (Sinh, Cosh, Tanh x inverse)
  - ✔️ Exponential (Exp, Log, Log2, Log10)
  - ✔️ Polynomial (Pow, Sqrt, Cbrt)
  - ✔️ Erf and gamma
- ✔️ Vector math basics
  - ✔️ Dot product
  - ✔️ Norm (sqrt(SumSquare))
- ✔️ Statistics
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
- ❔️ Filtering
  - ✔️ Convolution
  - ✔️ FFT
    - ✔️ R->C, C->C, C->C, C->R
    - ✔️ FFT shift
    - ✔️ Bin <-> Frequency conversions
  - ❔️ FIR filtering
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
  - ❔️ IIR filtering
    - Methods:
      - ✔️ Butterworth
      - ❌️ Chebyshev I
      - ❌️ Chebyshev II
      - ❌️ Elliptic
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
  - ✔️ Filter response analysis
    - ✔️ Compute amplitude & phase response
    - ✔️ Classify amplitude response: LP/HP/BP/BS
    - ✔️ Measure amplitude parameters for LP/HP/BP/BS
      - ✔️ Transition edges
      - ✔️ Stopband attenuation
      - ✔️ Passband ripple
  - ✔️ Polyphase FIR decomposition
  - ✔️ Resampling
    - ✔️ Decimation (every n-th)
    - ✔️ Expansion (zero-fill)
    - ✔️ Interpolation (polyphase)
    - ✔️ Arbitrary resampling (polyphase)
  - ✔️ Windowing
    - ✔️ Derived properties
      - ✔️ Gain
      - ✔️ Energy
    - ✔️ Functions
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
  - ❌️ Wavelets
    - ❌️ DWT
    - ❌️ CWT
    - ❌️ Wavelet design
- ❌️ Literals
  - ❌️ Hz, kHz, MHz, GHz, THz
  - ❌️ dB, dB10, dB20

## License

DSPBB is distributed under the MIT license, therefore can be used in commercial and non-commercial projects alike with very few restrictions. The dependencies of DSPBB are also distributed under premissive MIT or BSD licenses.