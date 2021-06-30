Rubik's Cube
============
This repo contains the C-library that I wrote during my ADACS internship in 2021.
It implements a small Rubik's cube model that is used for the `PUMAS muon transport engine`_.
PUMAS is used for mapping out the flux of muons detected in the Stawell gold mine in Victoria, Australia.

.. _PUMAS muon transport engine: https://github.com/niess/pumas

Installation
------------
This library can be easily built with the ``Makefile`` by simply using ``make`` in the root directory.

Model data
----------
The model requires a data file with 5 columns (``X``, ``Y``, ``Z``, ``density``, ``rock_id``) separated by either commas or spaces.
The data is assumed to be part of a regular 3D grid with no holes or missing cubes.
The coordinates in the data file themselves are assumed to change the fastest in the positive X-direction, then in the positive Y-direction and then in the negative Z-direction (e.g., for every Z-value, all possible Y-values are given and for every of these Y-values, all X-values are given).

Usage
-----
A model data file that satisfies the conditions given above can be loaded into the model using ``rubiks_cube_create``, which takes a ``struct rubiks_cube **`` and the name of the file to load.
Afterward, using ``rubiks_cube_find_cube``, one can find any cube within the loaded Rubik's cube model by providing the function with the coordinates and direction of the particle within the model for which the cube needs to be known.