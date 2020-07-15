from setuptools import setup

# get the version here
pkg_vars  = {}

with open("version.py") as fp:
    exec(fp.read(), pkg_vars)

setup(
    name='sn_stackers',
    version= pkg_vars['__version__'],
    description='Stackers for supernovae',
    url='http://github.com/lsstdesc/sn_stackers',
    author='Philippe Gris',
    author_email='philippe.gris@clermont.in2p3.fr',
    license='BSD',
    packages=['sn_stackers'],
    python_requires='>=3.5',
    zip_safe=False
)
