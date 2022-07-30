Digital Signal Processing Building Blocks
---

![Language](https://img.shields.io/badge/Language-C++17-blue)
[![License](https://img.shields.io/badge/License-MIT-blue)](#license)
[![Build & test](https://github.com/petiaccja/DSPBB/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/petiaccja/DSPBB/actions/workflows/build_and_test.yml)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=petiaccja_DSPBB&metric=alert_status)](https://sonarcloud.io/dashboard?id=petiaccja_DSPBB)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=petiaccja_DSPBB&metric=coverage)](https://sonarcloud.io/dashboard?id=petiaccja_DSPBB)

DSPBB is a modern C++17/20 library for digital signal processing. DSPBB is not geared towards any specific signal processing domain such as audio processing or seismology, it only aims to provide the building blocks that can help build any signal processing application.

### Key features
- Signal and system classes
- Arithmetic and math on signals
- Waveform generation
- Statistics
- Filtering:
  - FFT
  - IIR filter design
  - FIR flter design
  - Running signals through filters
  - Resampling
  - Windowing

### Development status

The library now implements most of the features that I was planning and is largely covered by unit tests. It's currently in a consolidation phase where I'm adding smaller features, ironing out problems, and improving ergonomy as I start using it in larger projects. As such, I'm expecting moderate, but breaking changes to the interface.

<a name="user_guide"></a>
## User guide
- [Introduction](docs/introduction.md)
- [Installation](docs/installation.md)
- [Features](docs/features.md)
- [Performance](docs/performance.md)
- [Numerical accuracy](docs/accuracy.md)

<a name="license"></a>
## License

DSPBB is distributed under the MIT license, therefore can be used in commercial and non-commercial projects alike with very few restrictions. The dependencies of DSPBB are also distributed under premissive MIT or BSD licenses.