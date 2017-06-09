# -*- coding: utf-8 -*-
r"""
Monkey patches to make the MacLane code work in standalone mode, i.e., without
modifying the sage source code at build time.
"""
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Fix doctests so they work in standalone mode (when invoked with sage -t, they run within the mac_lane/ directory)
import sys, os
if hasattr(sys.modules['__main__'], 'DC') and 'standalone' in sys.modules['__main__'].DC.options.optional:
    sys.path.append(os.getcwd())
    sys.path.append(os.path.dirname(os.getcwd()))

import valuation_space
from valuation_space import DiscretePseudoValuationSpace
import trivial_valuation
from trivial_valuation import TrivialValuation, TrivialPseudoValuation
import padic_valuation
from padic_valuation import pAdicValuation
import gauss_valuation
from gauss_valuation import GaussValuation
import value_group
from value_group import DiscreteValuationCodomain, DiscreteValueGroup, DiscreteValueSemigroup
import function_field_valuation
from function_field_valuation import FunctionFieldValuation
import augmented_valuation
from augmented_valuation import AugmentedValuation
import scaled_valuation
from scaled_valuation import ScaledValuation

# fix unpickling and type checks of classes (otherwise, the instances of the
# local file and the instances that come from the mac_lane import define
# different types)
from trivial_valuation import TrivialDiscreteValuation, TrivialDiscretePseudoValuation
from function_field_valuation import FunctionFieldValuation_base, DiscreteFunctionFieldValuation_base, RationalFunctionFieldValuation_base, InducedFunctionFieldValuation_base, ClassicalFunctionFieldValuation_base, FunctionFieldFromLimitValuation, InfiniteRationalFunctionFieldValuation, FiniteRationalFunctionFieldValuation, NonClassicalRationalFunctionFieldValuation, InfiniteRationalFunctionFieldValuation, FunctionFieldMappedValuation_base, FunctionFieldExtensionMappedValuation, RationalFunctionFieldMappedValuation
from limit_valuation import LimitValuation, MacLaneLimitValuation, LimitValuation_generic
from mapped_valuation import MappedValuation_base, FiniteExtensionFromLimitValuation, FiniteExtensionFromInfiniteValuation, MappedValuation_base
from augmented_valuation import FiniteAugmentedValuation, InfiniteAugmentedValuation
from gauss_valuation import GaussValuation_generic
from valuation import DiscretePseudoValuation, DiscreteValuation, InfiniteDiscretePseudoValuation
from padic_valuation import pAdicValuation_base, pAdicValuation_int, pAdicValuation_padic, pAdicFromLimitValuation
from developing_valuation import DevelopingValuation
from augmented_valuation import AugmentedValuation_base, FinalAugmentedValuation, NonFinalAugmentedValuation, FinalFiniteAugmentedValuation, NonFinalFiniteAugmentedValuation
from inductive_valuation import FiniteInductiveValuation, FinalInductiveValuation, InfiniteInductiveValuation, NonFinalInductiveValuation
from scaled_valuation import ScaledValuation_generic

# =============================
# FORMER MONKEY PATCHES TO SAGE
# =============================
r"""
Check that :trac:`23166` has been resolved::

sage: K.<x> = FunctionField(QQ)
sage: K(x) in K._ring # indirect doctest
True

Check that :trac:`23166` has been resolved::

sage: K.<x> = FunctionField(QQ)
sage: K(1) in QQ # indirect doctest
True

Check that :trac:`23167` has been resolved::

sage: R.<x> = QQ[]
sage: K.<x> = FunctionField(QQ)
sage: R.fraction_field().is_subring(K)
True

Check that :trac:`23185` has been resolved::

sage: R.<x> = QQ[]
sage: K.<x> = FunctionField(QQ)
sage: R.is_subring(K)
True
sage: R.is_subring(R.fraction_field())
True

Check that :trac:`23167` has been resolved::

sage: R.fraction_field().is_subring(K)
True

Check that :trac:`21879` has been resolved::

sage: QQ.coerce_map_from(ZZ).is_injective()
True

Check that :trac:`21879` has been resolved::

sage: Hom(ZZ,QQ['x']).natural_map().is_injective()
True

Check that :trac:`21879` has been resolved::

sage: R.<xbar> = R.quo(x^2+x+1)
sage: Hom(ZZ,R).natural_map().is_injective()
True

Check that :trac:`21879` has been resolved::

sage: R.<x> = QQbar[]
sage: R.coerce_map_from(QQbar).is_injective()
True

Check that rings embed into polynomial rings::

sage: G = GaussianIntegers()
sage: R.<x> = G[]
sage: G.hom(R).is_injective()
True

Check that :trac:`23186` has been resolved::

sage: QQ.coerce_map_from(ZZ).is_surjective()
False

Check that :trac:`23187` has been resolved::

sage: QQ['x'].coerce_map_from(ZZ['x']).is_injective()
True
sage: GF(2)['x'].coerce_map_from(ZZ['x']).is_surjective()
True

Check that :trac:`23188` has been resolved::

sage: R.<a> = ZqCA(9)
sage: R['x'].is_subring(R.fraction_field()['x'])
True

Check that :trac:`21879` has been resolved::

sage: GaussianIntegers().fraction_field().coerce_map_from(QQ).is_injective()
True
sage: GaussianIntegers().fraction_field().coerce_map_from(QQ).is_injective()
True

Check that :trac:`23189` has been resolved::

sage: CyclotomicField(5).maximal_order().coerce_map_from(ZZ).is_injective()

Check that :trac:`23190` has been resolved::

sage: R.<x> = ZZ[]
sage: S.<x> = QQ[]
sage: S.quo(x^2 + 1).coerce_map_from(R.quo(x^2 + 1)).is_injective()
True

Check that :trac:`23191` has been resolved::

sage: R.<x> = ZZ[]
sage: S = R.quo(x^2+x+1)
sage: S(1).inverse_of_unit()
1

"""


# ===========================================
# MONKEY PATCH BUGS/REQUIRED FEATURES IN SAGE
# ===========================================
import sage

# ===============================
# MONKEY PATCH FEATURES INTO SAGE
# ===============================

# Implement Qp/Zp.valuation
sage.rings.padics.padic_generic.pAdicGeneric.valuation = lambda self: pAdicValuation(self)

# Implement principal_part() and sides() for newton polygons
# TODO: inline this when these methods get called, instead of modifying Sage
r"""
TESTS::

    sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
    sage: NP = sage.geometry.newton_polygon.NewtonPolygon([(0,1),(1,0),(2,1)])
    sage: NP.principal_part()
    Infinite Newton polygon with 2 vertices: (0, 1), (1, 0) ending by an infinite line of slope 0

"""
import sage.geometry.newton_polygon
sage.geometry.newton_polygon.NewtonPolygon_element.principal_part = lambda self: sage.geometry.newton_polygon.NewtonPolygon(self.vertices(), last_slope=0)
sage.geometry.newton_polygon.NewtonPolygon_element.sides = lambda self: zip(self.vertices(), self.vertices()[1:])

######################################
# MONKEY PATCH FOR NON-TRIVIAL TESTING
######################################

# make some_elements() non-trivial for number fields
def some_elements(self):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: K = GaussianIntegers().fraction_field()
        sage: list(K.some_elements())
        [I, 0, 1, 1/2, 2*I, -I, -2, 0, 0]

    """
    for element in self.polynomial_ring().some_elements():
        yield element(self.gen())
sage.rings.number_field.number_field.NumberField_generic.some_elements = some_elements
del(some_elements)

# make some_elements() deterministic for function fields
def some_elements(self):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: list(K.some_elements()) == list(K.some_elements())
        True

    """
    for num in self._ring.some_elements():
        for den in self._ring.some_elements():
            if den != 0:
                yield self(num) / self(den)
sage.rings.function_field.function_field.RationalFunctionField.some_elements = some_elements
del(some_elements)

def some_elements(self):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: list(L.some_elements()) == list(L.some_elements())
        True

    """
    for element in self._ring.some_elements():
        yield self(element)
sage.rings.function_field.function_field.FunctionField_polymod.some_elements = some_elements
del(some_elements)

# make some_elements() non-trivial for fraction fields
def some_elements(self):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: K = R.fraction_field()
        sage: len(list(K.some_elements()))
        72

    """
    for num in self.ring().some_elements():
        for den in self.ring().some_elements():
            if den != 0:
                yield self(num) / self(den)
sage.rings.fraction_field.FractionField_generic.some_elements = some_elements

# make some_elements() non-trivial for orders in number fields
def some_elements(self):
    r"""
    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: R = GaussianIntegers()
        sage: list(R.some_elements())
        [I, 0, 1, 2*I, -I, -2, 0, 0]

    """
    for element in self.fraction_field().some_elements():
        if element in self:
            yield self(element)
sage.rings.number_field.order.Order.some_elements = some_elements
del(some_elements)

###########################################################
# PATCHING ONLY NECESSARY BECAUSE WE DO NOT RUN WITHIN SAGE
###########################################################

# register modules at some standard places so imports work as exepcted
r"""
sage: from sage.rings.valuation.gauss_valuation import GaussValuation
"""
import imp, sys
sage.rings.valuation = sys.modules['sage.rings.valuation'] = imp.new_module('sage.rings.valuation')
sage.rings.valuation.gauss_valuation = sys.modules['sage.rings.valuation.gauss_valuation'] = gauss_valuation
sage.rings.valuation.valuation = sys.modules['sage.rings.valuation.valuation'] = valuation
sage.rings.valuation.valuation_space = sys.modules['sage.rings.valuation.valuation_space'] = valuation_space
sage.rings.valuation.augmented_valuation = sys.modules['sage.rings.valuation.augmented_valuation'] = augmented_valuation
sage.rings.function_field.function_field_valuation = sys.modules['sage.rings.function_field.function_field_valuation'] = function_field_valuation

# fix unpickling of factories
from sage.structure.factory import register_factory_unpickle
register_factory_unpickle("pAdicValuation", pAdicValuation)
register_factory_unpickle("GaussValuation", GaussValuation)
register_factory_unpickle("TrivialValuation", TrivialValuation)
register_factory_unpickle("TrivialPseudoValuation", TrivialPseudoValuation)
register_factory_unpickle("FunctionFieldValuation", FunctionFieldValuation)
register_factory_unpickle("AugmentedValuation", AugmentedValuation)
register_factory_unpickle("LimitValuation", LimitValuation)
register_factory_unpickle("ScaledValuation", ScaledValuation)
