.. _customization:

**************************************
Customization and Global Configuration
**************************************

The :file:`configrc` file
=========================

This package uses a :file:`configrc` configuration file to customize
certain properties. You can control a number of key features of such as
where your data will download to. HERMES packages look for this configuration file
in a platform specific directory, which you can see the path for by running::

  >>> import hermes_core
  >>> hermes_core.print_config()  # doctest: +SKIP

Using your own :file:`configrc` file
=====================================
To maintain your own customizations, you must place your customized :file:`configrc` inside the appropriate configuration folder (which is based off the operating system you are working on). The `AppDirs module <https://github.com/sunpy/sunpy/blob/main/sunpy/extern/appdirs.py>`_  provided by the `sunpy` package is used to figure out where to look for your configuration file. 

.. warning::
    Do not edit the configrc file directly in the Python package as it will get overwritten  every time you re-install or update the package.

You can copy the file below, customize it, and then place your customized :file:`configrc` file inside your config folder.

If you work in our developer environment you can place your configuration file in this directory:

.. code-block:: bash

  /home/vscode/.config/hermes_core/

If you do not use our developer environment, you can run the following code to see where to place it on your specific machine as well:

.. doctest::

  >>> from hermes_core import util
  >>> print(util.config._get_user_configdir())
  /home/vscode/.config/hermes_core


.. note:: 
  For more information on where to place your configuration file depending on your operating system, you can refer to the `AppDirs module docstrings <https://github.com/sunpy/sunpy/blob/1459206e11dc0c7bfeeeec6aede701ca60a8630c/sunpy/extern/appdirs.py#L165>`_. 


To learn more about how to set-up your development environment see :ref:`dev_env`.

See below (:ref:`configrc-sample`) for an example configuration file.

.. _customizing-with-dynamic-settings:

Dynamic settings
================

You can also dynamically change most of the default settings. One setting that cannot be changed is the location of the log file which is set on import. All settings are stored in a Python ConfigParser instance called ``hermes_core.config``, which is global to the package. Settings can be modified directly, for example::

    import hermes_core
    hermes_core.config.set('downloads', 'download_dir', '/home/user/Downloads')


.. _configrc-sample:

A sample configrc file
--------------------------------------------------------------------

.. only:: html

    `(download) <../_static/configrc>`__

.. literalinclude:: ../../hermes_core/data/configrc