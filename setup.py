# -*- coding: utf-8 -*-
import toml
from pathlib import Path
from setuptools import setup


def read_pyproject_toml():
    pyproject_toml = Path(__file__).parent / "pyproject.toml"
    pyproject_dct = toml.loads(pyproject_toml.read_text())
    return pyproject_dct


pyproject_dct = read_pyproject_toml()

package_data = {"": ["*"]}
exclude_package_data = {"": ["*.db", "*.gz", "test"]}
python_requires = pyproject_dct["tool"]["poetry"]["dependencies"].pop("python")

install_requires = [
    f'{depend}{version.replace("^", ">=")}'
    for depend, version in pyproject_dct["tool"]["poetry"]["dependencies"].items()
]

packages = [pyproject_dct["tool"]["poetry"]["name"]]
package_data = {"": ["*"]}
entry_points = {"console_scripts": ["makura = makura.assembly:main"]}

author, author_email = [
    row.strip(" <>") for row in pyproject_dct["tool"]["poetry"]["authors"][0].split("<")
]


def get_readme():
    return Path(pyproject_dct["tool"]["poetry"]["readme"]).read_text()


setup_kwargs = {
    "name": pyproject_dct["tool"]["poetry"]["name"],
    "version": pyproject_dct["tool"]["poetry"]["version"],
    "description": pyproject_dct["tool"]["poetry"]["description"],
    "long_description": get_readme(),
    "author": author,
    "author_email": author_email,
    "maintainer": author,
    "maintainer_email": author_email,
    "url": pyproject_dct["tool"]["poetry"]["repository"],
    "packages": packages,
    "package_data": package_data,
    "install_requires": install_requires,
    "entry_points": entry_points,
    "python_requires": python_requires,
}

setup(**setup_kwargs)
