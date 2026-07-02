/* CRPropa3 SWIG interface (for Python) */

/* Content:
 *
 * 1. SWIG settings and workarounds
 * 2. SWIG and CRPropa headers
 * 3. Pretty print for Python
 *
 */


/* 1. SWIG settings and workarounds */
%include "1_swig.i"

/* 2: SWIG and CRPropa headers */
%include "2_headers.i"

/* 3. Pretty print for Python */

%define __REPR__( classname )

%pythoncode %{
globals()["classname"["classname".find('::')+2:]].__repr__ = globals()["classname"["classname".find('::')+2:]].getDescription
%}

%enddef

%pythoncode %{
Vector3d.__repr__ = Vector3d.getDescription
Vector3f.__repr__ = Vector3f.getDescription
%}

%include "3_repr.i"