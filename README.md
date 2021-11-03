Digital Signal Processing Building Blocks
---

![Build](https://github.com/petiaccja/DSPBB/workflows/Build/badge.svg)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=petiaccja_DSPBB&metric=alert_status)](https://sonarcloud.io/dashboard?id=petiaccja_DSPBB)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=petiaccja_DSPBB&metric=coverage)](https://sonarcloud.io/dashboard?id=petiaccja_DSPBB)

I couldn't really find a **low-level digital signal processing** (DSP) library that is 
- **free** (BSD, MIT)
- **modern** (C\++17 - C\++20)
- **maintained**
- **nice to use**

... so I decided to make one.

At this point, the interfaces change frequently, features come and go, so I don't suggest you use this library for anything other than your own entertainment. However, if the features below appeal to you, stay tuned!




## Features

More like notes to myself as to what I want to have.

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
- ❔️ Common math functions
  - ✔️ Complex (Abs, Arg, Real, Imag)
  - ✔️ Trigonometric (Sin, Cos, Tan x inverse)
  - ✔️ Hyperbolic (Sinh, Cosh, Tanh x inverse)
  - ✔️ Exponential (Exp, Log, Log2, Log10)
  - ✔️ Polynomial (Pow, Sqrt, Cbrt)
  - ❌️ Special math (Erf, Gamma, etc.)
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
  - ❌️ IIR filtering
    - Methods:
      - ❌️ Butterworth
      - ❌️ Chebyshev I
      - ❌️ Chebyshev II
      - ❌️ Elliptic
    - Types:
      - ❌️ Low-pass
      - ❌️ High-pass
      - ❌️ Band-pass
      - ❌️ Band-stop
      - ❌️ Notch
    - Realizations:
      - ❌️ Direct form I.
      - ❌️ Direct form II.
      - ❌️ Cascaded biquad (I. & II.)
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
  - ❌️ Wavelets
    - ❌️ DWT
    - ❌️ CWT
    - ❌️ Wavelet design
  - ❔️ Windowing
    - ✔️ Derived properties
      - ✔️ Gain
      - ✔️ Energy
    - ❔️ Functions
      - ✔️ Rectangular
      - ✔️ Triangular
      - ✔️ Hamming
      - ✔️ Blackman
      - ✔️ Blackman-Harris
      - ✔️ Flat top
      - ❌️ Kaiser
      - ✔️ Gaussian
      - ❌️ Dolph-Chebyshev
      - ❌️ Lanczos
- ❌️ Literals
  - ❌️ Hz, kHz, MHz, GHz, THz
  - ❌️ dB, dB10, dB20
