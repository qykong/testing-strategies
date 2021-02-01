from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize(["pyscripts/listdict.pyx",
                           "pyscripts/simulations.py",
                           "pyscripts/tracer.py"],  language_level="3")
)
