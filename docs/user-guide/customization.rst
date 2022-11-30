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
To maintain your own customizations place a your customized :file:`configrc` in the directory of your project. The `AppDirs module <https://github.com/sunpy/sunpy/blob/main/sunpy/extern/appdirs.py>`_  provided by the `sunpy` package is used to figure out what file to be used. 

If you work in our devcontainer environment you can place your configuration file in this directory:

.. exec_code::

   # Location to place configuration file:

   # --- hide: start ---
   from hermes_core import util
   print(util.config._get_user_configdir())
   #hide:toggle

.. Note:: 

    The above directory is provided by the :code:`_get_user_configdir` function in the hermes-core util module within the :file:`config.py` file. You can find this file `here <https://github.com/HERMES-SOC/hermes_core/blob/main/hermes_core/util/config.py>`_.
  
Do not edit the default file directly as every time you install or update, this file will be overwritten.

To learn more about how to set-up your development environment to be able to use your own configuation file see :ref:`dev_env`.

See below (:ref:`configrc-sample`) for the example config file.

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

A sample sunpyrc file
--------------------------------------------------------------------

.. only:: html

    `(download) <../_static/configrc>`__

.. literalinclude:: ../../hermes_core/data/configrc