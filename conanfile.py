from conans import ConanFile, CMake, tools


class DSPBBConan(ConanFile):
    name = "dspbb"
    license = "MIT"
    author = "Péter Kardos"
    url = "https://github.com/petiaccja/DSPBB"
    description = "A library for digital signal processing fundamentals. Very WIP."
    topics = ("dsp", "digital-signal-processing", "fir-filter", "fft", "convolution")
    exports_sources = "include/*"

    def package(self):
        self.copy("*.hpp", dst="include", src="include")
        self.copy("*.h", dst="include", src="include")