.. _customization:

**************************************
Customization and Global Configuration
**************************************

HERMES & SWxSOC Configuration
=============================

The `hermes_core` package uses an inherited configuration file from the `~swxsoc.util.config` module.
The SWxSOC configuration is made accessible through the `~hermes_core.config` property that loads the appropriate parent configuration for the HERMES mission. 

You can view the default configuration in the following way:

  >>> import swxsoc
  >>> import hermes_core
  >>> swxsoc.print_config(hermes_core.config) # doctest:+ELLIPSIS
  FILES USED:
  ...
  CONFIGURATION:
  ...


Using A Custom Configuration File
=================================

You can also use a custom configuration file if needed that can be overridden when loading the parent SWxSOC configuration which is shared with this HERMES package.
For more information on how to override the SWxSOC configuration please see the package documentation here `Using your own config.yml file <https://swxsoc.readthedocs.io/en/latest/user-guide/customization.html#using-your-own-config-yml-file>`_ 
