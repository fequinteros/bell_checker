import setuptools

with open("README.md", "r", encoding="utf-8") as readme:
    description = readme.read()

setuptools.setup(
    name='bell_checker',
    version='0.0.1',
    author='FeQuinteros',
    description='A Qiskit module to construct and evaluate bell inequalities',
    long_description=description,
    long_description_content_type="text/markdown",
    url='https://github.com/fequinteros/bell_checker',
    packages=['bell_checker'],
    install_requires=['qiskit', 'numpy', 'scipy', 'matplotlib'],
    include_package_data=True)
