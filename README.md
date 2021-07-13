Digital Signal Processing Building Blocks
---

![Build](https://github.com/petiaccja/DSPBB/workflows/Build/badge.svg)
[![codecov](https://codecov.io/gh/petiaccja/DSPBB/branch/master/graph/badge.svg?token=8A2I59KJ5H)](https://codecov.io/gh/petiaccja/DSPBB)

I couldn't really find a **low-level digital signal processing** (DSP) library that is 
- **free** (BSD, MIT)
- **modern** (C\++14 - C\++20)
- **maintained**
- **nice to use**

... so I decided to make one.

At this point, the interfaces change frequently, features come and go, so I don't suggest you use this library for anything other than your own entertainment. However, if the features below appeal to you, stay tuned!




## Features

More like notes to myself as to what I want to have.

### Design

- ❌️ No dependencies
  - ❌️ XSimd is optional
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
- ❔️ Generators
  - ✔️ Wave
    - ✔️ Sine
    - ✔️ Square
    - ✔️ Sawtooth (fw, bw, tri, any)
    - ✔️ PWM
  - ✔️ Chirp/Sweep
  - ❔️ Space
    - ✔️ Linspace
    - ❌️ Logspace
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
- ❔️ Statistics
  - ✔️ Sum
  - ✔️ SumSquare
  - ✔️ RootMeanSquare
  - ✔️ Mean
  - ✔️ Stddev
  - ❌️ Skewness
  - ❌️ Kurtosis
  - ❌️ Covariance
  - ❌️ Correlation
- ❔️ Filtering
  - ✔️ Convolution
  - ✔️ FFT
  - ❔️ FIR filtering
    - ✔️ Window method
    - ❌️ Least-squares method
    - ❌️ Parks-McClellan method
    - ✔️ Low-pass
    - ❌️ High-pass
    - ❌️ Band-pass
    - ✔️ Arbitrary response
    - ❌️ Quality assessment tools
  - ❌️ IIR filtering
    - ❌️ Butterworth
    - ❌️ Chebyshev I
    - ❌️ Chebyshev II
    - ❌️ Elliptic
    - ❌️ Quality assessment tools
  - ✔️ Polyphase FIR decomposition
  - ✔️ Resampling
    - ✔️ Decimation (every n-th)
    - ✔️ Expansion (zero-fill)
    - ✔️ Interpolation (polyphase)
    - ✔️ Arbitrary resampling (polyphase)
  - ❌️ Hilbert filter
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
      - ❌️ Triangular
      - ✔️ Hamming
      - ❌️ Blackman
      - ❌️ Blackman-Harris
      - ✔️ Flat top
      - ❔️ Kaiser (c++17 only)
      - ❌️ Gaussian
      - ❌️ Dolph-Chebyshev
      - ❌️ Lanczos
- ❌️ Literals
  - ❌️ Hz, kHz, MHz, GHz, THz
  - ❌️ dB, dB10, dB20
