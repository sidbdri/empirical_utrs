from setuptools import setup

import empirical_utrs

setup(
    name='empirical_utrs',
    version=empirical_utrs.__version__,
    url='https://github.com/sidbdri/empirical_utrs',
    license='MIT License',
    author='Owen Dando',
    author_email='owen.dando@ed.ac.uk',
    packages=['empirical_utrs'],
    install_requires=[
        'docopt==0.6.2',
        'numpy',
        'pandas==0.13.0',
        'pysam==0.14.1',
        'python-dateutil==2.4.2',
        'pytz==2016.3',
        'schema==0.3.1',
        'six==1.10.0',
    ],
    scripts=[
        'bin/get_empirical_utrs',
        'bin/get_utr_sequences',
    ]
)
