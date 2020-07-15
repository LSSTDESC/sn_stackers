from setuptools import setup


setup(
    name='sn_stackers',
    version=__version__,
    description='Stackers for supernovae',
    url='http://github.com/lsstdesc/sn_stackers',
    author='Philippe Gris',
    author_email='philippe.gris@clermont.in2p3.fr',
    license='BSD',
    packages=['sn_stackers'],
    python_requires='>=3.5',
    zip_safe=False
)
