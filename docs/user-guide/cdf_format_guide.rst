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

==============
3. Conventions
==============

All HERMES scientific data products that will be shared between HERMES entities (e.g.
ITFs, IDS groups) or made available to the general research community will be stored as
CDF data files and are expected to be compatible with CDF version 3.5. Data that will
not be shared beyond an individual team may be stored in any format that is convenient
for that team.

--------------------------------------
3.1 Science Product Naming Conventions
--------------------------------------

The HERMES data products will be produced with the following filename format where
the individual identifying components are described in Table 3-1. Additionally, to ensure
software compatibility between disparate systems, filenames will consist of all lowercase
characters. Filenames are used as a system identifier for a logical grouping of data and
are also stored in the `Logical_file_id` global attribute field (see Section 4.1.8). It is
expected that filenames will be created dynamically from the attributes identified in
Section 4 of this document.

**Filename Format**
    `scId_instrumentId_mode_dataLevel_optionalDataProductDescriptor_startTime_vX.Y.Z.ext`

.. list-table:: Table 3-1: Filename Component Description
   :widths: 25 50 25
   :header-rows: 1

   * - Short Name
     - Description
     - Valid Options
   * - scID
     - Spacecraft ID
     - `hermes`
   * - instrumentId
     - Instrument or investigation identifier shortened to three letter acronym.
     - `eea`, `mrt`, `nms`, `spn`
   * - mode
     - *TBS*
     - *TBS*
   * - dataLevel
     - The level to which the data product has been processed
     - `l0`, `l1`, `ql`, `l2`, `l3`, `l4`
   * - optionalDataProductDescriptor
     - This is an optional field that may not be needed for all products. Where it is used, identifier should be short (e.q. 3-8 characters) descriptors that are helpful to end-users. If a descriptor contains multiple components, underscores are used to separate those components.
     - An optional time span may be specified as "2s" to represent a data file that spans two seconds. In this case, "10s" and "5m" are other expected values that correspond with ten seconds and 5 minutes respectively.
   * - startTime
     - The start time of the contained data given in "YYYYMMDD_hhmmss"
     - `20220601_101520`
   * - vX.Y.Z
     - The 3-part version number of the data product. Full description of this identifier is provided in Section 3.1.1 of this document.
     - `v0.0.0`, `v<#.#.#>`
   * - .ext
     - The required file extension, where CDF is required.
     - `.cdf`

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
3.1.1 Version Numbering Guidelines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The three-part version number contains the interface number, quality number, and bug
fix/revision number. The initial release of CDF data that is suitable for scientific
publication should begin with “v1.Y.Z”. Each component of the version number is
incremented in integer steps, as needed, and Table 3-2 describes the instances in which
the value should be incremented. Release “v0.Y.Z” may be used for early development
purposes.

====================
4. Global Attributes
====================

.. csv-table:: HERMES Global Metadata Schema
   :file: global_attributes.csv
   :widths: 30, 70, 30, 30, 30, 30, 30
   :header-rows: 1