
try:
    from setuptools import setup
except ImportError:
    pass

setup(
    name='halflife',
    version='0.1dev',
    packages=['halflife',],
    license='MIT',
    long_description=open('README.md').read(),
)
