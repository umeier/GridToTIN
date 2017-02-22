from distutils.core import setup
import numpy as np
from Cython.Build import cythonize

setup(
    name='GridToTIN',
    version='0.1',
    packages=['grid2tin'],
    url='',
    license='MIT',
    author='Ulrich Meier',
    author_email='ulrich.meier@ldbv.bayern.de',
    description='',
    ext_modules = cythonize("grid2tin/*.pyx"),
    include_path = [np.get_include()]
)
