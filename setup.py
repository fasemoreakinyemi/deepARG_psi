import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="deepARG_psi", # Replace with your own username
    version="0.0.1",
    author="Akinyemi Mandela Fasemore",
    author_email="fasemoremandela@gmail.com",
    description="A python library to classify protein sequences of Antibiotic resitiance nature into AB groups based on psi-blast",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
