.. cdf_format_guide:

*******************************
HERMES CDF Format Guide
*******************************

===============
1. Introduction
===============

The :py:class:`~hermes_core.util.schema.HERMESDataSchema` class provides an interface to
examine the HERMES CDF Format Guide.

---------------------
1.1 Purpose and Scope
---------------------

This document is provided as a reference for construction of HERMES standard CDF
files. It is intended to complement information available from the Space Physics Data
Facility (listed in Sec. 1.2). It lays down REQUIREMENTS and
RECOMMENDATIONS for Level 2 (and above) CDF files that are intended for public
access, and should be taken as RECOMMENDATIONs for all other mission CDFs.
This document is based on discussions within the HERMES Science Data Working
Group (HSDWG) and personnel at NASA's Space Physics Data Facility (SPDF). It is
intended to provide sufficient reference material to understand CDF files and
the requirements for creating HERMES CDF files, and to understand the structure and
contents of the resulting CDF files.

--------------
1.2 References
--------------

Relevant documents that provide background material and support details provided in
this guide are listed below:

- `SPDF CDF User's Guide <http://cdf.gsfc.nasa.gov/>`_
- `SKTEditor <http://spdf.gsfc.nasa.gov/sp_use_of_cdf.html>`_
- `ISTP Guidelines <http://spdf.gsfc.nasa.gov/istp_guide/istp_guide.html>`_
- `ISTP/IACG Global Attributes <http://spdf.gsfc.nasa.gov/istp_guide/gattributes.html>`_

================================
2. HERMES Science Investigations
================================

The HERMES Instrument Suite will make high-time resolution measurements of plasmas
(ions and electrons) and magnetic fields. The HERMES Instrument Suite consists of the
following complement of instruments:

- Electron Electrostatic Analyzer (EEA): The EEA provides measurements of
    low-energy electrons in the solar wind and in Earth’s deep magnetotail by
    measuring electron flux as functions of energy and direction.
- Miniaturized Electron pRoton Telescope (MERIT): The MERiT instrument
    measures the flux of high-energy electrons and ions with two telescopes
    pointing in opposite directions and nominally spanning the forward and
    reverse Parker Spiral.
- Noise Eliminating Magnetometer In a Small Integrated System (NEMISIS):
    NEMISIS is comprised of a fluxgate magnetometer (M0) at the end of a
    deployable boom and two inductive magnetometers (M1, M2) mounted on
    the HERMES platform. Each sensor measures the vector magnetic field at
    its location. Measurements from the 3 sensors are combined to reduce the
    contribution to the local field due to Gateway.
- Solar Probe Analyzer for Ions (SPAN-I): The SPAN-i ion sensor measures
    Interplanetary and Magnetotail ion flux as functions of direction and
    energy/charge from several eV/q to 20 keV/q. A time-of-flight section
    enables it to sort particles by their mass/charge ratio, permitting
    differentiation of ion species.

HERMES Instrument Team Facilities (ITFs) are the principal institutions associated with
each of the HERMES science investigations. These facilities and their personnel provide
support to the operation of their instruments and the overall data processing and
distribution effort for HERMES science data products. The institutions listed in Table
2-1 have responsibility for each of the investigations and their corresponding instruments.

.. table:: Table 2-1 HERMES ITF Summary
   :widths: 20 30 30

   +-------------------------------------------------------------------------+---------------------------------------------------------------------------+---------------------------------+
   | HERMES Investigation                                                    | Managing Institution                                                      | Principal Investigator          |
   +=========================================================================+===========================================================================+=================================+
   | Electron Electrostatic Analyzer (EEA)                                   | Goddard Space Flight Center (GSFC)                                        | D. Gershman                     |
   +-------------------------------------------------------------------------+---------------------------------------------------------------------------+---------------------------------+
   | Miniaturized Electron pRoton Telescope (MERIT)                          | Goddard Space Flight Center (GSFC)                                        | S. Kanekal                      |
   +-------------------------------------------------------------------------+---------------------------------------------------------------------------+---------------------------------+
   | Noise Eliminating Magnetometer In a Small Integrated System (NEMISIS)   | Goddard Space Flight Center (GSFC), University of Michigan                | E. Zesta, M. Moldwin (Co-I),    |
   +-------------------------------------------------------------------------+---------------------------------------------------------------------------+---------------------------------+
   | Solar Probe Analyzer for Ions (SPAN-I)                                  | University of California, Berkeley (UCB), Space Sciences Laboratory (SSL) | R. Livi                         |
   +-------------------------------------------------------------------------+---------------------------------------------------------------------------+---------------------------------+
