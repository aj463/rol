// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PYROL_TEUCHOS_CUSTOM
#define PYROL_TEUCHOS_CUSTOM

#include <Teuchos_ParameterList.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>

namespace py = pybind11;

template<typename T>
Teuchos::Array< T > copyNumPyToTeuchosArray(pybind11::array_t<T> array) {

    auto np_array = array.template mutable_unchecked<1>();
    int size = array.shape(0);
    Teuchos::Array< T > av(size);
    for (int i=0; i < size; ++i)
      av[i] = np_array(i);
    return av;
}

template<typename T>
pybind11::array_t<T> copyTeuchosArrayToNumPy(Teuchos::Array< T > & tArray) {

    pybind11::array_t<T> array(tArray.size());
    auto data = array.template mutable_unchecked<1>();
    for (int i=0; i < tArray.size(); ++i)
      data(i) = tArray[i];
    return array;
}

// Implementation based on:
// https://github.com/trilinos/Trilinos/tree/master/packages/PyTrilinos/src/PyTrilinos_Teuchos_Util.cpp
bool setPythonParameter(Teuchos::RCP<Teuchos::ParameterList> plist,
			const std::string      & name,
			py::object             value);

inline void setPythonParameterOrThrow(Teuchos::RCP<Teuchos::ParameterList> plist,
                                      const std::string & name,
                                      py::object value)
{
  if (!setPythonParameter(plist, name, value))
    throw py::type_error("Cannot store a value of type '" +
                         py::type::of(value).attr("__name__").cast<std::string>() +
                         "' in ParameterList entry '" + name +
                         "'; supported types are bool, int, float, str and ParameterList");
}


template <typename T>
bool setPythonParameterArray(Teuchos::RCP<Teuchos::ParameterList> plist,
			const std::string      & name,
			pybind11::array_t< T >   value)
{
  auto tArray = copyNumPyToTeuchosArray(value);
  plist->set(name, tArray);
  return true;
}

// Implementation based on:
// https://github.com/trilinos/Trilinos/tree/master/packages/PyTrilinos/src/PyTrilinos_Teuchos_Util.cpp
py::object getPythonParameter(Teuchos::RCP<Teuchos::ParameterList> plist,
			      const std::string            & name);

template <typename T>
void def_ParameterList_member_functions(T cl) {
  cl.def("__setitem__", [](Teuchos::RCP<Teuchos::ParameterList> &m, const std::string &name, Teuchos::ParameterList value) { m->set(name,value);  });
  cl.def("set", [](Teuchos::RCP<Teuchos::ParameterList> &m, const std::string &name, Teuchos::ParameterList value) { m->set(name,value);  });
  cl.def("sublist", [](Teuchos::RCP<Teuchos::ParameterList> &m, const std::string &name) { if (!m->isSublist(name)) { throw py::key_error(name); } return pybind11::cast(sublist(m, name)); }, pybind11::return_value_policy::reference);
  cl.def("__setitem__", [](Teuchos::RCP<Teuchos::ParameterList> &m, const std::string &name, pybind11::object value) { setPythonParameterOrThrow(m,name,value);  });
  cl.def("__getitem__", [](Teuchos::RCP<Teuchos::ParameterList> &m, const std::string &name) {
    // Sublist
    if (m->isSublist(name))
      return py::cast(Teuchos::sublist(m, name));
    return getPythonParameter(m,name);
  });
  cl.def("set", [](Teuchos::RCP<Teuchos::ParameterList> &m, const std::string &name, pybind11::object value) { setPythonParameterOrThrow(m,name,value);  });
  cl.def("get", [](Teuchos::RCP<Teuchos::ParameterList> &m, const std::string &name) {
    // Sublist
    if (m->isSublist(name))
      return py::cast(Teuchos::sublist(m, name));
    return getPythonParameter(m,name);
  });
}
#endif // PYROL_TEUCHOS_CUSTOM
