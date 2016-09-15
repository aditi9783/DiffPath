from distutils.core import setup

# setup.py syntax - https://docs.python.org/2/distutils/setupscript.html


# Running setup.py from PyCharm - https://www.jetbrains.com/help/pycharm/2016.1/creating-and-running-setup-py.html
# Package setup - https://docs.python.org/2/distutils/examples.html
## Whether to use specific modules or packages - http://programmers.stackexchange.com/questions/243044/single-python-file-distribution-module-or-package

# Installing packages - https://docs.python.org/2/install/

# For uploading to PyPl - https://docs.python.org/2/distutils/packageindex.html

# Specifying project dependencies:
## http://python-packaging.readthedocs.io/en/latest/dependencies.html
## https://packaging.python.org/requirements/

setup(
    name='DiffPath',
    version='0.0.1',
    py_modules=['DiffPath', 'fdrcorrection', 'geneexp', 'graphs', 'pathway'],
    url='https://github.com/aditi9783/DiffPath',
    license='MIT',
    author='Aditi Gupta',
    author_email='aditi9783@gmail.com',
    description='Differential Pathways',
    install_requires=[
        'matplotlib',
        'networkx',
        'numpy',
        'scipy'
    ]
)
