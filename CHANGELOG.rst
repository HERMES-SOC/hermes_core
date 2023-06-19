This project uses `semantic versioning <https://semver.org>`_.

Latest
======
* Added data class to hold measurements and to save to CDF files

0.2.0 (2023-03-22)
========
This release includes improvements tested in the second dataflow test. Since the last release, the improvements are as follows:

- Uses fstrings instead of the format syntax
- Adds a log message on import to show the version number of the package
- Documentation content and styling improvements
- Switches from using `setup.py` to `pyproject.toml` for package
- Bug fixes for supporting `pathlib`'s Path objects, and a permissions bug to the devcontainer

* Using f-strings and log on import by @ehsteve in https://github.com/HERMES-SOC/hermes_core/pull/28
* Added latest versions of python to testing by @ehsteve in https://github.com/HERMES-SOC/hermes_core/pull/29
* Add to the Documentation the location for where config files should be store (dynamically) by @dbarrous in https://github.com/HERMES-SOC/hermes_core/pull/35
* Fix devcontainer config to use new vscode user by @dbarrous in https://github.com/HERMES-SOC/hermes_core/pull/32
* Update to docs, added logo, updated theme colors and favicon by @ehsteve in https://github.com/HERMES-SOC/hermes_core/pull/37
* Fix to version number and move to pyproject.toml usage by @ehsteve in https://github.com/HERMES-SOC/hermes_core/pull/40
* Bug fix by @ehsteve in https://github.com/HERMES-SOC/hermes_core/pull/41


0.1.0 (2022-10-05)
==================
This version release was tested in the first HERMES Ground System data flow test.

* First draft of python packaging including sphinx documentation based on the sunpy package template
* First draft of the documentation including coding standards for the HERMES ecosystem
* Automated testing and coverage using GitHub actions
* Logging support
* Configuration support
* Utilities parsing compliant filenames for level 0 binary files and creating and parsing higher level filenames