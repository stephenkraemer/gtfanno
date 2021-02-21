from setuptools import setup

setup(
        name='gtfanno',
        version='0.1',
        author='Stephen Kraemer',
        author_email='stephenkraemer@gmail.com',
        license='MIT',
        packages = ['gtfanno'],
        install_requires=[
            'pandas',
            'numpy',
            'pybedtools',
        ],
        python_requires='>=3.8',
)



