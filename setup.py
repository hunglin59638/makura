# -*- coding: utf-8 -*-
from pathlib import Path
from setuptools import setup

packages = ["makura"]

package_data = {"": ["*"]}

install_requires = [
    "ete3>=3.1.2,<4.0.0",
    "pandas>=1.5.1,<2.0.0",
    "requests>=2.28.1,<3.0.0",
    "tqdm>=4.64.1,<5.0.0",
]


def get_readme():
    return Path("README.md").read_text()


entry_points = {"console_scripts": ["makura = makura.assembly:main"]}

setup_kwargs = {
    "name": "makura",
    "version": "1.0.2",
    "description": "Makura: NCBI Genome downloader",
    "long_description": get_readme(),
    "author": "Hung-Lin, Chen",
    "author_email": "hunglin59638@gmail.com",
    "maintainer": "Hung-Lin, Chen",
    "maintainer_email": "hunglin59638@gmail.com",
    "url": "https://github.com/hunglin59638/makura",
    "packages": packages,
    "package_data": package_data,
    "install_requires": install_requires,
    "entry_points": entry_points,
    "python_requires": ">=3.8,<4",
}


setup(**setup_kwargs)
