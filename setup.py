from setuptools import setup,find_packages
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "tempowork",
    version = "0.0.23",
    author = "Ali Jazayeri",
    author_email = "ali.jazayeri@drexel.edu",
    description = ("A library for mining temporal networks where time is represented as a continuous dimension!"),
    license = "MIT",
    keywords = "subgraph mining; temporal networks",
    url = "https://github.com/alijazayeri/tempowork",
    py_modules=['tempowork'],
    packages=find_packages(),
    long_description=read('README.rst'),
long_description_content_type='text/x-rst'
)