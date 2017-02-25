#ifndef GEN_RENDER_H
#define GEN_RENDER_H

#ifdef PYTHON_RENDERING_SUPPORTED
	#include <Python.h>
#endif

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>

#include "resources.h"

namespace gen {
namespace render {

void drawMap(std::vector<char> &drawdata, std::string filename);

#ifdef PYTHON_RENDERING_SUPPORTED
	PyObject* _get_PyString_FromString(const char *v);
	void _checkPyObjectNotNull(PyObject *obj, std::string err);
	void _checkPySuccess(int ret, std::string err);
#endif

}
}

#endif