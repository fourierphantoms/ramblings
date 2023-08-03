Signal processing. C++, Python.


To use C++ from Python, we need to use a binding mechanism. There are several options available, such as ctypes, Boost.Python, SWIG, and pybind11. In this case, we'll use pybind11, which is a lightweight header-only library that exposes C++ types in Python and vice versa, mainly to create Python bindings for a C++ code.

To start, we need to install pybind11.

Install pybind11 using pip:

pip install pybind11


Next, let's write our C++ code with the LOWESS (LOcally WEighted Scatterplot Smoothing) implementation and the Python bindings.

The following is a simple C++ implementation of the LOWESS algorithm with pybind11 bindings:

c++
Copy code
#include <vector>
#include <algorithm>
#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

double tricube_weight(double x) {
    double absx = std::abs(x);
    if (absx > 1) return 0;
    double p = 1 - absx * absx * absx;
    return p * p * p;
}

std::vector<double> lowess(const std::vector<double>& x, const std::vector<double>& y, double bandwidth, int iterations) {
    std::vector<double> y_smoothed(y.size());
    std::vector<double> weights(x.size());
    for (int iteration = 0; iteration < iterations; ++iteration) {
        for (std::size_t i = 0; i < x.size(); ++i) {
            double sum_weight = 0, sum_x = 0, sum_y = 0, sum_xy = 0, sum_xx = 0;
            for (std::size_t j = 0; j < x.size(); ++j) {
                weights[j] = iteration == 0 ? 1 : tricube_weight((y[j] - y_smoothed[i]) / bandwidth);
                double weight = weights[j];
                double xj = x[j] - x[i];
                sum_weight += weight;
                sum_x += weight * xj;
                sum_y += weight * y[j];
                sum_xy += weight * xj * y[j];
                sum_xx += weight * xj * xj;
            }
            double slope = (sum_weight * sum_xy - sum_x * sum_y) / (sum_weight * sum_xx - sum_x * sum_x);
            double intercept = (sum_y - slope * sum_x) / sum_weight;
            y_smoothed[i] = slope * (x[i] - x[i]) + intercept;
        }
    }
    return y_smoothed;
}

PYBIND11_MODULE(lowess, m) {
    m.def("lowess", &lowess, py::arg("x"), py::arg("y"), py::arg("bandwidth") = 0.25, py::arg("iterations") = 3);
}

In the above code, we have a function lowess that accepts two vectors x and y (the data to be smoothed), a bandwidth parameter for controlling the smoothing, and an iterations parameter for controlling the number of smoothing iterations. The function calculates the smoothed y-values and returns them as a vector. The PYBIND11_MODULE macro creates a module named lowess that exposes the lowess function to Python.

The C++ code should be compiled to a shared library that can be imported in Python. For this, you will need a C++ compiler that is compatible with the Python interpreter. The exact command will depend on your system, but it might look something like this:

c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` lowess.cpp -o lowess`python3-config --extension-suffix`


This command compiles the lowess.cpp file to a Python extension named lowess. The -O3 -Wall options enable optimization and warnings. The -shared -std=c++11 -fPIC options are required for creating a shared library that can be loaded by Python. The `python3 -m pybind11 --includes` part adds the necessary include paths for pybind11 and Python. The -o lowess part specifies the output file name. Finally, the `python3-config --extension-suffix` part adds the appropriate file extension for a Python extension (e.g., .so or .pyd).

After compiling the C++ code, you can use it in Python like this here:

#include <vector>
#include <algorithm>
#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

double tricube_weight(double x) {
    double absx = std::abs(x);
    if (absx > 1) return 0;
    double p = 1 - absx * absx * absx;
    return p * p * p;
}

std::vector<double> lowess_rolling(const std::vector<double>& x, const std::vector<double>& y, int window_size, double bandwidth, int iterations) {
    std::vector<double> y_smoothed(y.size(), 0.0);
    for (std::size_t i = 0; i < x.size(); ++i) {
        std::vector<double> weights(window_size, 1.0);
        std::size_t window_start = i < window_size / 2 ? 0 : i - window_size / 2;
        std::size_t window_end = std::min(window_start + window_size, x.size());
        for (int iteration = 0; iteration < iterations; ++iteration) {
            double sum_weight = 0, sum_x = 0, sum_y = 0, sum_xy = 0, sum_xx = 0;
            for (std::size_t j = window_start; j < window_end; ++j) {
                double xj = x[j] - x[i];  // xj is now centered around 0
                if (iteration > 0) {
                    double yj_diff = y[j] - y_smoothed[j];
                    weights[j - window_start] = tricube_weight(yj_diff / (6.0 * bandwidth));
                }
                double weight = weights[j - window_start];
                sum_weight += weight;
                sum_x += weight * xj;
                sum_y += weight * y[j];
                sum_xy += weight * xj * y[j];
                sum_xx += weight * xj * xj;
            }
            double denominator = sum_weight * sum_xx - sum_x * sum_x;
            if (std::abs(denominator) < 1e-6) {
                y_smoothed[i] = sum_y / sum_weight;
            } else {
                double slope = (sum_weight * sum_xy - sum_x * sum_y) / denominator;
                double intercept = (sum_y - slope * sum_x) / sum_weight;
                y_smoothed[i] = intercept;
            }
        }
    }
    return y_smoothed;
}

PYBIND11_MODULE(lowess, m) {
    m.def("lowess_rolling", &lowess_rolling, py::arg("x"), py::arg("y"), py::arg("window_size"), py::arg("bandwidth") = 0.25, py::arg("iterations") = 3);
}

##### 

#include <vector>
#include <algorithm>
#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

double tricube_weight(double x) {
    double absx = std::abs(x);
    if (absx > 1) return 0;
    double p = 1 - absx * absx * absx;
    return p * p * p;
}

std::vector<double> lowess_rolling(const std::vector<double>& x, const std::vector<double>& y, int window_size, double bandwidth, int iterations) {
    std::vector<double> y_smoothed(y.size(), 0.0);
    for (std::size_t i = 0; i < x.size(); ++i) {
        std::vector<double> weights(window_size, 1.0);
        std::size_t window_start = i < window_size / 2 ? 0 : i - window_size / 2;
        std::size_t window_end = std::min(window_start + window_size, x.size());
        for (int iteration = 0; iteration < iterations; ++iteration) {
            double sum_weight = 0, sum_x = 0, sum_y = 0, sum_xy = 0, sum_xx = 0;
            for (std::size_t j = window_start; j < window_end; ++j) {
                double xj = x[j] - x[i];  // xj is now centered around 0
                if (iteration > 0) {
                    double yj_diff = y[j] - y_smoothed[j];
                    weights[j - window_start] = tricube_weight(yj_diff / (6.0 * bandwidth));
                }
                double weight = weights[j - window_start];
                sum_weight += weight;
                sum_x += weight * xj;
                sum_y += weight * y[j];
                sum_xy += weight * xj * y[j];
                sum_xx += weight * xj * xj;
            }
            double denominator = sum_weight * sum_xx - sum_x * sum_x;
            if (std::abs(denominator) < 1e-6) {
                y_smoothed[i] = sum_y / sum_weight;
            } else {
                double slope = (sum_weight * sum_xy - sum_x * sum_y) / denominator;
                double intercept = (sum_y - slope * sum_x) / sum_weight;
                y_smoothed[i] = intercept;
            }
        }
    }
    return y_smoothed;
}

PYBIND11_MODULE(lowess, m) {
    m.def("lowess_rolling", &lowess_rolling, py::arg("x"), py::arg("y"), py::arg("window_size"), py::arg("bandwidth") = 0.25, py::arg("iterations") = 3);
}










In this Python code, we first generate some noisy sine wave data. Then we call the lowess function to smooth the data. Finally, we plot the original data and the smoothed data.

Please note that the LOWESS implementation in the C++ code is a simple version and may not be as robust or flexible as implementations in mature scientific computing libraries. The main goal of this example is to show how to use pybind11 to call C++ code from Python.
