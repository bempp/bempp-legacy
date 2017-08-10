"""This module gives access to various standard shapes."""

from bempp.api.shapes.shapes import \
    sphere, ellipsoid, cube, almond, regular_sphere
from bempp.api.shapes.shapes import reentrant_cube
from bempp.api.shapes.shapes import rectangle_with_hole
from bempp.api.shapes.shapes import union
from bempp.api.shapes.letters import letters_bempp as bempp
from bempp.api.shapes.fractals import sierpinski_pyramid, menger_sponge
from bempp.api.shapes.fractals import \
    sierpinski_triangle, sierpinski_carpet, koch_snowflake


__all__ = ['sphere', 'ellipsoid', 'cube',
           'almond', 'regular_sphere', 'reentrant_cube',
           'bempp', 'sierpinski_pyramid', 'menger_sponge',
           'sierpinski_triangle', 'sierpinski_carpet', 'koch_snowflake',
           'rectangle_with_hole','union']
