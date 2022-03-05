#!/usr/bin/python
# -*- coding: utf-8 -*-

r"""
Python script to create the source for the ``database_cubic_hecke`` package
in ``py`` format in the given path. This utility should be used to switch
to a new version of the data files if the original files have changed.

.. NOTE::

    This function needs `SymPy <https://de.wikipedia.org/wiki/SymPy>`__

INPUT:

- ``path`` -- name of path where to store the generated  files
  (if not provided, the the subdirectory ``database_cubic_hecke``)
  under the current working directory is used)

EXAMPLES::

    python create_marin_data.py
    python create_marin_data.py ~
"""

##############################################################################
#       Copyright (C) 2021 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################

import sys
from os import environ
from os.path import join, exists
from sympy import var, eye
from sympy import  __version__ as sympy_version
from database_cubic_hecke import FileNames, marin_url
from database_cubic_hecke import __version__ as this_version

cmdline_args = sys.argv[1:]

path = None

if len(cmdline_args) > 0:
    path = cmdline_args[0]
    if not exists(path):
        raise ValueError("the path '%s' does not exist" % path)

if not path:
    pwd = environ['PWD']
    path = join(pwd, 'database_cubic_hecke')
    if not exists(path):
        path = pwd


# -----------------------------------------------------------------------------
# Templates
# -----------------------------------------------------------------------------
py_header_text = r"""# ------------------------------------------------------------
# This file was generated using create_marin_data.py
#
# on    : Version %s
# under : SymPy %s
# using : Python %s
# ------------------------------------------------------------

"""
function_name = '_read_%s'
template=r"""# generated data function for cubic Hecke algebra
def %s(%s):
    %s
    data = %s
    return data

"""

rem=r"""
    This code has been generated by ``create_marin_data.py`` (from
    the `database_cubic_hecke repository  <https://github.com/soehms/database_cubic_hecke>`__),
    please don't edit.
"""

doc=r"""{}
    Return precomputed basis of Ivan Marin
%s
    OUTPUT:

    A list of lists where each of it represents a braid pre image
    in Tietze form.

    EXAMPLES::

        >>> from database_cubic_hecke import %s
        >>> res = %s()
        >>> len(res)
        648
    {}
""".format('r"""', '"""')

docrep=r"""{}
    Return precomputed representation matrices of Ivan Marin
%s
    INPUT:

    ``%s`` -- values for the variables of the representation matrices

    OUTPUT:

    A triple ``dim_list, repr_list, repr_list_inv`` where each member
    is a list indexed by the number of the representation.

    ``dim_list`` is a list of integers representing the dimension of
    the corresponding representation
    ``repr_list`` is a triple of dictionaries each describing the
    representation matrix of the corresponding generator of the
    cubic Hecke algebra.
    ``repr_list_inv`` is a triple of dictionaries each describing the
    inverse of the correspondig matrix of ``repr_list``.

    EXAMPLES::

        >>> from database_cubic_hecke import %s
        >>> dim_list, repr_list, repr_list_inv = %s(%s)
        >>> dim_list
        %s
        >>> len(repr_list)
        %s
        >>> g1, g2, g3 = repr_list[0]
        >>> list(g1.items())[0]
        %s
    {}
""".format('r"""', '"""')


# -----------------------------------------------------------------------------
# Auxillary function to generate the function code
# -----------------------------------------------------------------------------
def _write_data_function(filename, data):
        r"""
        Generate code for functions to access the data.
        """
        from textwrap import fill
        def fl(s):
            return fill(str(s), subsequent_indent='        ')
        func_name = function_name %filename.cache()

        if filename == FileNames.basis:
            var = ''
        elif filename == FileNames.irred_split:
            var = 'a, b, c, j'
            varexp = '3, 5, 7, 11'
            res1 = '[1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 8, 8, 8, 9, 9]'
            res2 = '24'
            res3 = '((0, 0), 3)'
        else:
            var = 'u, v, w'
            varexp = '3, 5, 7'
            res1 = '[648]'
            res2 = '1'
            res3 = '((0, 27), -5)'

        if var:
            docstr  = docrep %(rem, var, func_name,  func_name, varexp, res1, res2, res3)
            func  = template %(func_name,  var, docstr,  fl(str(data)))
        else:
            docstr  = doc %(rem, func_name, func_name)
            func  = template %(func_name,  '', docstr,  fl(str(data)))


        header = py_header_text %(this_version.value, sympy_version, sys.version.split(' ')[0])

        import_filen = join(path, filename.py())
        with open(import_filen, 'w') as f:
            f.write(header)
            f.write(func)



# -----------------------------------------------------------------------------
# Auxillary function to load the data from Ivan Marins web page
# -----------------------------------------------------------------------------
def _load_data(filename):
    r"""
    Download a data file in maple syntax.
    """
    # import directly from the internet page
    from six.moves.urllib.request import urlopen
    try:
        from urllib.error import HTTPError, URLError
    except ImportError:
        from urllib2 import HTTPError

    try:
        download = join(marin_url, filename)
        download_data = urlopen(download).read().decode()
        preparsed_data =download_data.replace(':=', '=').replace(';', '').replace('^', '**')
        return preparsed_data
    except (HTTPError, URLError):
        raise IOError('Data import file %s not found! Internet connection needed!' %(filename))
    return None


# -----------------------------------------------------------------------------
# Functions to create the Python files
# -----------------------------------------------------------------------------
def create_marin_split_repr():
    r"""
    Create the data base for split irreducible representations of the
    cubic Hecke algebra according to the original data from Iwan Marin's
    home page.
    """
    # ------------------------------------------------------
    # Ivan Marin's data file uses a, b, c for the variables
    # corresponding to the eigenvalues of the cubic equation.
    # ------------------------------------------------------
    fname = FileNames.irred_split
    preparsed_data = _load_data(fname.maple())

    global a, b, c, j          # set in exec
    a, b, c, j = var('a, b, c, j')
    exec(preparsed_data, globals())

    u, v, w = (-a-b-c, a*b+b*c+a*c, -a*b*c)

    def invert(matr, dim):
       """
       Return inverse matrix for generators
       """
       matri = v/w*eye(dim)
       matri += u/w*matr
       matri += 1/w*matr**2
       matri.simplify()
       return -matri

    # ----------------------------------------------------------------------
    # Restoring the split irreducibles from Iwan Marin's homepage
    # ----------------------------------------------------------------------
    from sympy import Matrix

    anz_reps = len(reps)

    dim_list = []
    repr_list = []
    repr_list_inv = []
    for i in range(anz_reps):
        repi = reps[i]
        if len(repi) != 3:
           raise RuntimeError( 'Error at position %d: three generators expected, got: %d' %( i, len(repi)))
        mt = []
        mtI = []
        for j in range(3):
            mat = Matrix(repi[j])
            d, e = mat.shape
            matI = invert(mat, d)
            mt.append(mat.todok())
            mtI.append(matI.todok())

        dim_list.append(d)
        repr_list.append( mt )
        repr_list_inv.append( mtI )

    _write_data_function(fname, (dim_list, repr_list, repr_list_inv))


def create_static_db_marin_regular(right=False):
    r"""
    Create the data base for regular representations of the cubic
    Hecke algebra according to the original data from Iwan Marin's home page.
    """

    if right == False:
        fname = FileNames.regular_left
    else:
        fname = FileNames.regular_right
    preparsed_data = _load_data(fname.maple())

    global baseH4, Matrix, u, v, w, mm1, mm2, mm3, mm1I, mm2I, mm3I, reps   # set in exec
    u, v, w = var('u, v, w')

    dim_list = [None]
    from sympy import Matrix as sMatrix
    def Matrix(d, e, data):
        dim_list[0]=d
        return sMatrix(data)

    exec(preparsed_data, globals())

    repr_list  = [[mm1.todok(), mm2.todok(), mm3.todok()]]
    repr_list_inv = [[mm1I.todok(), mm2I.todok(), mm3I.todok()]]

    if not right:
        _write_data_function(FileNames.basis, baseH4)
    _write_data_function(fname, (dim_list, repr_list, repr_list_inv))


# -------------------------------------------------------------------------------------
# Create all files
# -------------------------------------------------------------------------------------

create_marin_split_repr()
create_static_db_marin_regular()
create_static_db_marin_regular(right=True)
