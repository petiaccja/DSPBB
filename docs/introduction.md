# Introduction

## Why another library?

There are many libraries out there for DSP:
- [Aquila](https://github.com/zsiciarz/aquila)
- [KFR](https://github.com/kfrlib/kfr)
- [Q](https://github.com/cycfi/Q)
- [Intel MKL](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html)
- and several others.

However, they are either paid for commercial use, are not maintained anymore, lack some features or are tailored towards a specific field. I needed something different, but I encourage you to check out the libraries listed above as well.

What DSPBB bring to the table:
- **Template-based**: You can specify `float` or `double` anywhere. Some things will work for integrals and custom types (i.e. fixed point, arbitrary precision) as well, but as `std::complex` is not specialized for these, support is limited.
- **Low-level, generic design**: DSPBB does not focus on any specific field. No crossfeed, no equalizer, just the boring math needed to implement these quickly.
- **Stateless API**: Almost all functionality is provided by free functions, simplistic classes are only used to organize data. This minimizes memory allocations to make embedded development easier
- **Composability**: The API does not model continuous streams or processing pipelines. You can use the free function API as building blocks for your own stream and pipeline implementations, or use them as is.
- **Built-in filter design**: You can design your filters straight from C++ in real-time, no need to import them from Matlab. Of course, you can import filters designed offline, if you have to.
- **Extendibility**: If a feature, like a filter design method, is missing, chances are it's not so difficult to implement in user code using the rest of the library as a basis.

## A little example

Let's look at a simplified crossfeed implementation:
```c++
#include <dspbb/Primitives/Signal.hpp>
#include <dspbb/Primitives/SignalView.hpp>
#include <dspbb/Filtering/IIR.hpp>

using namespace dspbb;

auto Crossfeed(SignalView<const float> left, SignalView<const float> right)
  -> std::pair<Signal<float>, Signal<float>>
{
    // Design filter.
    const auto desc = Lowpass(BUTTERWORTH).Cutoff(0.16f);
    const size_t order = 2;
    const auto zpk = IirFilter<float>(order, desc);
    const auto filter = CascadedBiquad{ zpk };

    // Apply filter.
    CascadedForm<float> state{ filter.Order() };
    const auto lpLeft = Filter(left, filter, state);
    state.Reset();
    const auto lpRight = Filter(left, filter, state);

    // Combine results.
    const float gain = 0.25f;
    const auto outLeft = left + gain*lpRight;
    const auto outRight = right + gain*lpLeft;
    return { outLeft, outRight };
}
```

### 1. 
`Signal` is the most important class of the library, it stores the discrete samples of a signal. `SignalView` is a utility you can use to call DSPBB functions on any memory area, including a subset of a `Signal`.

### 2.
Generally, you want to have some filter to apply to your signal. For this, you can design IIR filters (like in this example) or FIR filters with DSPBB, or you can design them elsewhere (i.e. in Matlab) and copy the coefficients yourself.

Here, we chose to design a simple second-order Butterworth lowpass filter as a discrete zero-pole-gain system, which we promptly convert to a cascaded biquad representation.

### 3.
For IIR filters, which are recursive, we also need to store some state. The state corresponds to the realizations such as direct form I or direct form II. In this case, we will use a special state that is only applicable to biquad cascades.

Once you have the input signal, the filter, and the state, you can just call `Filter` to apply the filter to the input. If you want to feed some more signal into the filter and continue where you left off, just reuse the same state object, otherwise, you should reset it.

### 4.
Just some basic arithmetic to mix the filtered signals with the originals. Do I need to explain anything here?

## Conclusion

So you want to give it a go? Head to the [installation section](installation.md).

Want to know more? Check out the [features](features.md) or take a look at the [examples](../examples).