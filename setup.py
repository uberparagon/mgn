# nodoctest

from distutils.core import setup

try:
    from sage.env import SAGE_SRC, SAGE_VERSION
except ImportError:
    import warnings
    warnings.warn("This package only works with SageMath!")

version = "1.1.1"

print("This is version {}".format(version))

with open("description.rst") as f:
    long_description = f.read()
    
setup(
    name='mgn',
    version=version,
    packages=['strataalgebra', 'topintersections'],
    url='https://github.com/uberparagon/mgn',
    license='MIT',
    author='Drew Johnson',
    author_email='werd2.718@gmail.com',
    description='A Sage program for computing products, Faber-Zagier relations, and top intersections on the moduli space of stable curves.',
    long_description = long_description
)
