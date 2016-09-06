from setuptools import setup, find_packages

setup(
    name='Auriferous',
    version='0.1.0',
    description='Variant filtering tools',
    license='MIT',
    author='Dominic Fitzgerald',
    author_email='dominicfitzgerald11@gmail.com',
    url='https://github.com/djf604/auriferous',
    download_url='https://github.com/djf604/chunky-pipes/tarball/0.2.3',
    packages=find_packages(),
    entry_points={
        'console_scripts': ['auri = auriferous:execute_from_command_line']
    },
    install_requires=['PyVCF'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries',
        'License :: OSI Approved :: MIT License'
    ]
)
