import setuptools
import subprocess
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name='ccpbsa',
    version='0.1',
    description="A structure based method quantitive estimation of mutational \
    free energy",
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=setuptools.find_packages(),
    author='Linkai Zhang',
    author_email='linkai.zhang1@googlemail.com',
    scripts=['ccpbsa/ccpbsa', 'ccpbsa/ccpbsa-setup'],
    include_package_data=True,
    install_requires=["numpy", "pandas", "pymol", "tqdm"],
    zip_safe = False
)

subprocess.run(['ccpbsa-setup'])
