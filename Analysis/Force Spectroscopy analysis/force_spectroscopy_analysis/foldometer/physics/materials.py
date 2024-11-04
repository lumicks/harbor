#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import namedtuple

from tweezer.lux import SellmeierCoefficients, dispersion

MaterialDescription = namedtuple("MaterialDescription", ['n', 'density', 'SellmeierCoefficients'])


class Material:

    def __init__(self, name):
        self.name = name


Units = namedtuple('Units', ['density'])
units = Units('g/cm^3')

#silica = Material()
#polystyrene = Material()
#water = Material()

# Glasses ----
# BK7 - borosilicate crown glass
_SellmeierBK7 = SellmeierCoefficients(1.03961212, 0.231792344, 1.01046945, 0,
                                      6.00069867E-3, 2.00179144E-2, 1.03560653E2, 0)

_IndexOfRefractionBK7At1064 = dispersion(1.064, _SellmeierBK7)

BK7 = MaterialDescription(_IndexOfRefractionBK7At1064, 1.05E-21, _SellmeierBK7)

# Silica, Silicon dioxide, quartz
_SellmeierSilica = SellmeierCoefficients(0.28604141, 1.07044083, 1.10202242, 0,
                                         0, 1.00585997e-2, 100, 0)

_IndexOfRefractionSilicaAt1064 = dispersion(1.064, _SellmeierSilica)

Silica = MaterialDescription(_IndexOfRefractionSilicaAt1064, 3, _SellmeierSilica)

# Polymer materials ----
# Polystyrene
_SellmeierPolystyrene = SellmeierCoefficients(1.4435, 0, 0, 0,
                                              0.020216, 0, 0, 0)

_IndexOfRefractionPolystyreneAt1064 = dispersion(1.064, _SellmeierPolystyrene)

Polystyrene = MaterialDescription(_IndexOfRefractionPolystyreneAt1064, 1.05E-21, _SellmeierPolystyrene)


# Liquids
_SellmeierWater = SellmeierCoefficients(5.666959820E-1, 1.731900098E-1, 2.095951857E-2, 1.125228406E-1,
                                        5.084151894E-3, 1.818488474E-2, 2.625439472E-2, 1.073842352E1)

_IndexOfRefractionWaterAt1064 = dispersion(1.064, _SellmeierWater)

Water = MaterialDescription(_IndexOfRefractionWaterAt1064, 1.05E-21, _SellmeierWater)
