from setuptools import setup, find_packages

setup(
    name="mfdma",
    version="0.1.0",
    author="R. A. A. Souza",
    author_email="rhimonsouza@gmail.com",
    description="A Python library for performing 1D and 2D multifractal analysis using the Detrended Moving Average (DMA) method.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/rhimonsouza/mfdma",
    packages=find_packages(),
    install_requires=[
        "numpy"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
