#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "Protread.hpp"

namespace py = pybind11;

py::tuple wrap_calculate_map_sparse(py::array_t<double> input, double threshold)
{
    auto buf = input.request();
    int n = buf.shape[0];

    auto contacts = calculate_map_sparse(static_cast<double *>(buf.ptr), n, threshold);

    py::array_t<int> rows(contacts.size());
    py::array_t<int> cols(contacts.size());

    auto rows_ptr = static_cast<int *>(rows.request().ptr);
    auto cols_ptr = static_cast<int *>(cols.request().ptr);

    for (size_t k = 0; k < contacts.size(); ++k)
    {
        rows_ptr[k] = contacts[k].i;
        cols_ptr[k] = contacts[k].j;
    }

    return py::make_tuple(rows, cols, n);
}

PYBIND11_MODULE(Protread, m)
{
    m.doc() = "Protein Contact Map Generator";

    m.def("generate_contact_map_sparse", &wrap_calculate_map_sparse, "Calculate a sparse contact map from coords");

    m.attr("__version__") = "0.0.2";
}