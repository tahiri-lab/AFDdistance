from setuptools import find_packages, setup

setup(
    name="afd",
    version="0.0.0",
    description="Ancestor Frequency Difference distance for mutation trees",
    author="Elio Torquet",
    packages=find_packages(),
    install_requires=["treelib==1.8.0", "scipy==1.16.0"],
)
