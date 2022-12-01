.. _customization:

**************************************
Customization and Global Configuration
**************************************

The :file:`configrc` file
=========================

This package uses a :file:`configrc` configuration file to customize
certain properties. You can control a number of key features of such as
where your data will download to. sunpy looks for this configuration file
in a platform specific directory, which you can see the path for by running::

  >>> import hermes_core
  >>> hermes_core.print_config()  # doctest: +SKIP

Using your own :file:`configrc` file
=====================================
To maintain your own customizations, you must place your customized :file:`configrc` inside the appropriate configuration folder (Which is based off the operating system you are working from). The `AppDirs module <https://github.com/sunpy/sunpy/blob/main/sunpy/extern/appdirs.py>`_  provided by the `sunpy` package is used to figure out where to look for your configuration file. 

.. warning::
    Do not edit the configrc file directly in the Python package as it will get overwritten  every time you re-install or update the package.

For example, on a Linux system, the configuration folder is :file:`~/.config/hermes_core/`. You can then place your customized :file:`configrc` file inside this folder.

More specifically if you work in our devcontainer environment you can place your configuration file in this directory:

.. code-block:: bash

  /home/vscode/.config/hermes_core/

You can also run the following code-block to see where to place it on your specific machine as well:

.. doctest::

  >>> from hermes_core import util
  >>> print(util.config._get_user_configdir())

  /home/vscode/.config/hermes_core/


.. note:: 
  To get more information on where to place your configuration file depending on your operating system, you can refer to the `AppDirs module docstrings <https://github.com/sunpy/sunpy/blob/1459206e11dc0c7bfeeeec6aede701ca60a8630c/sunpy/extern/appdirs.py#L165>`_. 


To learn more about how to set-up your development environment see :ref:`dev_env`.

See below (:ref:`configrc-sample`) for an example configuration file.

.. _customizing-with-dynamic-settings:

Dynamic settings
================

You can also dynamically change the default settings in a Python script or
interactively from the python shell. All of the settings are stored in a
Python ConfigParser instance called ``sunpy.config``, which is global to
the package. Settings can be modified directly, for example::

    import hermes_core
    hermes_core.config.set('downloads', 'download_dir', '/home/user/Downloads')


.. _configrc-sample:

A sample configrc file
--------------------------------------------------------------------

.. only:: html

    `(download) <../_static/configrc>`__

.. literalinclude:: ../../hermes_core/data/configrc