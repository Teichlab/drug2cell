from setuptools import setup, find_packages

setup(
    name='drug2cell',
    version='0.1.1',
    description='Gene group activity utility functions for scanpy',
    url='https://github.com/Teichlab/drug2cell',
    packages=find_packages(exclude=['docs', 'notebooks']),
    install_requires=[
        'anndata',
        'pandas',
        'numpy',
        'statsmodels',
        'scipy',
        'blitzgsea'
    ],
    package_data={
        "drug2cell": ["*.pkl"]
    },
    author='Krzysztof Polanski, Kazumasa Kanemaru',
    author_email='kp9@sanger.ac.uk',
    license='non-commercial license'
)
