import setuptools

setuptools.setup(
    name='ccpbsa',
    description="A structure based method quantitive estimation of mutational \
    free energy",
    packages=setuptools.find_packages(),
    author='Linkai Zhang',
    author_email='linkai.zhang1@googlemail.com',
    entry_points={'console_scripts':['ccpbsa = ccpbsa:main']},
    include_package_data=True
)
