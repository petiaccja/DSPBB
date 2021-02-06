Digital Signal Processing Building Blocks
---

![Build](https://github.com/petiaccja/DSPBB/workflows/Build/badge.svg)

Yeah, I made this library out of frustration that I couldn't find something that does low-level DSP stuff, for free, is maintained, and is modern. Well I guess this library is half-baked right now, so please don't use it, unless you wanna put up with API changes. Maybe it will turn into something useful one day.


## Features

More like notes to myself as to what I want to have.

- ✓ Primitives:
  - ✓ Signal
    - ✓ Real
    - ✓ Complex
    - ✓ 4 operators
  - ✓ SignalView
- Generators
  - ✓ Sine
  - Square
  - Sawtooth
  - Triangle
  - PWM
- Common math functions
  - ✓ Log
  - ✓ Abs
  - ✓ Real
  - ✓ Imag
  - Exp
  - Angle
- Vector math basics
  - ✓ Dot product
  - ✓ Norm
- Statistics
  - ✓ Sum
  - ✓ SumSquare
  - ✓ RootMeanSquare
  - ✓ Mean
  - ✓ Stddev
  - Skewness
  - Covariance
  - Correlation
- Filtering
  - ✓ Convolution
  - ✓ FFT
  - FIR
    - ✓ Window method
    - Parks-McClellan method
    - ✓ Low-pass
    - High-pass
    - Band-pass
    - ✓ Arbitrary response
    - Quality assessment tools
  - IIR
    - Butterworth
    - Chebyshev I
    - Chebyshev II
    - Elliptic
    - Quality assessment tools
  - ✓ Polyphase FIR
  - Resampling
    - Decimation
    - Integer downsampling
    - Integer upsampling
    - ✓ Arbitrary rational resampling
  - Hilbert
  - Wavelets
    - DWT
    - CWT
    - Wavelet design
  - Windowing
    - Derived properties
      - Gain
      - Energy
    - Functions
      - Rectangular
      - Triangular
      - ✓ Hamming
      - Blackman
      - Blackman-Harris
      - Flat top
      - ✓ Kaiser
      - Gaussian
      - Dolph-Chebyshev
      - Lanczos
- Literals
  - Hz, kHz, MHz, GHz, THz
  - dB, dB10, dB20
