import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="rna-seeker",
    version="0.0.1",
    author="Shane Drabing",
    author_email="shane.drabing@gmail.com",
    packages=setuptools.find_packages(),
    url="https://github.com/shanedrabing/rna-seeker",
    description="Exploratory RNA-Seq data collection and k-means clustering analysis.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Microsoft :: Windows",
    ],
    data_files=[
        ("", ["LICENSE.txt"])
    ],
    install_requires=[
        "matplotlib", "requests",
    ]
)
