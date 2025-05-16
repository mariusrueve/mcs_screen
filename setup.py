from setuptools import setup, find_packages

setup(
    name="mcs_screen",
    version="1.0.0",
    author="Marius Rueve",
    author_email="marius.rueve@live.de",
    description="MCS Screening of two sets of molecules",
    packages=find_packages(),
    entry_points={"console_scripts": ["mcs_screen = mcs_screen.main:main"]},
)
