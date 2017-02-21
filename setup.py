from distutils.core import setup

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
    ext_modules = cythonize("grid2tin/*.pyx")
)
