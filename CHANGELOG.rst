This project uses `semantic versioning <https://semver.org>`_.

Latest
======
* Fixed automatic version numbering using `setuptools-scm <https://pypi.org/project/setuptools-scm/>`_.
* Removed use of towncrier for CHANGELOG generation
* Moved to using pyproject.toml file


0.1.0 (2022-10-05)
==================
This version release was tested in the first HERMES Ground System data flow test.

* First draft of python packaging including sphinx documentation based on the sunpy package template
* First draft of the documentation including coding standards for the HERMES ecosystem
* Automated testing and coverage using GitHub actions
* Logging support
* Configuration support
* Utilities parsing compliant filenames for level 0 binary files and creating and parsing higher level filenames