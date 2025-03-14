from setuptools import setup

setup(
    name='nemus_magpie',
    version='0.0.5',
    packages=['magpie'],
    url='https://github.com/nemus-Project/magpie-python',
    license='GNU General Public License v3.0',
    author='NEMUS',
    author_email='michele.ducceschi@unibo.it',
    description='An open source framework for working with thin plates under generalized elastic boundary conditions',
    install_requires=['numpy', 'scipy','matplotlib']
)
