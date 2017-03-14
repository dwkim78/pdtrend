from setuptools import find_packages, setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name='pdtrend',
    version='0.2.3',
    description='Photometric Detrending Algorithm',
    long_description=readme(),
    platforms=['any'],
    packages=find_packages(),
    include_package_data=True,
    url='https://github.com/dwkim78/pdtrend',
    license='MIT',
    author='Dae-Won Kim',
    author_email='dwkim78@gmail.com',
    install_requires=['numpy>=1.10', 'scikit-learn>=0.17', 'scipy>=0.17',
                      'matplotlib>=1.5'],
    keywords=['astronomy', 'light curves', 'time series',
              'machine learning', 'trend', 'detrend'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Topic :: Scientific/Engineering :: Astronomy'
    ]
)
