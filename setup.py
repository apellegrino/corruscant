from setuptools import setup, find_packages, Extension

def make_extension(lib, src_files):
    return Extension(
                'corruscant/bin/%s' % lib,
                sources=["src/%s" % f for f in src_files],
                include_dirs=["src"],
                extra_compile_args=[
                    '-Wall', '-Wextra', '-Werror', '-O3', '-march=native',
                    '-mtune=native',
                    ],
                )

setup(
    name="corruscant",
    version="0.7.0",
    description="A Python package for computing two-point correlation "
                "functions of astronomical data",
    # long_description = 
    url="https://github.com/apellegrino/corruscant",
    author="Andrew Pellegrino",
    author_email="ap994@drexel.edu",
    license="GPLv3",

    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: C",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
                 ],

    keywords="astronomy clustering correlation autocorrelation "
             "auto-correlation cross-correlation",

    ext_modules=[
        make_extension("libkdtree", ["kdbuild.c", "kdquery.c"]),
        make_extension("libcoords", ['coordmath.c']),
    ],

    packages=["corruscant", "corruscant.twopoint"],

    install_requires=[
                "numpy>=1.8",
                "healpy",
                "pymetis",
                ],
)
    
