# -*- coding: utf-8 -*-
r"""
Inductive valuations on polynomial rings

This module provides functionality for inductive valuations, i.e., finite
chains of :class:`AugmentedValuation`s on top of a :class:`GaussValuation`.

AUTHORS:

- Julian Rüth (2016-11-01): initial version

REFERENCES:

.. [ML1936] Mac Lane, S. (1936). A construction for prime ideals as absolute
values of an algebraic field. Duke Mathematical Journal, 2(3), 492-510.

.. [ML1936'] MacLane, S. (1936). A construction for absolute values in
polynomial rings. Transactions of the American Mathematical Society, 40(3),
363-395.

"""
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from valuation import DiscreteValuation, InfiniteDiscretePseudoValuation
from developing_valuation import DevelopingValuation

from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method

class InductiveValuation(DevelopingValuation):
    r"""
    Abstract base class for iterated :class:`AugmentedValuation` on top of a
    :class:`GaussValuation`.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, pAdicValuation(QQ, 5))

    TESTS::

        sage: TestSuite(v).run() # long time

    """
    def is_equivalence_unit(self, f, valuations=None):
        r"""
        Return whether ``f`` is an equivalence unit, i.e., an element of
        :meth:`effective_degree` zero (see [ML1936'] p.497.)

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_equivalence_unit(x)
            False
            sage: v.is_equivalence_unit(S.zero())
            False
            sage: v.is_equivalence_unit(2*x + 1)
            True

        """
        f = self.domain().coerce(f)

        if f.is_zero():
            return False
        return self.effective_degree(f, valuations=valuations) == 0

    def equivalence_reciprocal(self, f, coefficients=None, valuations=None, check=True):
        r"""
        Return an equivalence reciprocal of ``f``.

        An equivalence reciprocal of `f` is a polynomial `h` such that `f\cdot
        h` is equivalent to 1 modulo this valuation (see [ML1936'] p.497.)

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation which is an
          :meth:`equivalence_unit`

        - ``coefficients`` -- the coefficients of ``f`` in the :meth:`phi`-adic
          expansion if known (default: ``None``)

        - ``valuations`` -- the valuations of ``coefficients`` if known
          (default: ``None``)

        - ``check`` -- whether or not to check the validity of ``f`` (default:
          ``True``)

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R = Zp(3,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = 3*x + 2
            sage: h = v.equivalence_reciprocal(f); h # optional: integrated (needs xgcd for polynomials with p-adic coefficients)
            2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)
            sage: v.is_equivalent(f*h, 1) # optional: integrated
            True

        In an extended valuation over an extension field::

            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.augmentation(x^2 + x + u, 1)
            sage: f = 2*x + u
            sage: h = v.equivalence_reciprocal(f); h
            (u + 1) + O(2^5)
            sage: v.is_equivalent(f*h, 1)
            True

        Extending the valuation once more::

            sage: v = v.augmentation((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
            sage: h = v.equivalence_reciprocal(f); h
            (u + 1) + O(2^5)
            sage: v.is_equivalent(f*h, 1)
            True

        TESTS:

        A case that caused problems at some point::

            sage: K = Qp(2, 4)
            sage: R.<x> = K[]
            sage: L.<a> = K.extension(x^4 + 4*x^3 + 6*x^2 + 4*x + 2)
            sage: R.<t> = L[]
            sage: v = GaussValuation(R)
            sage: w = v.augmentation(t + 1, 5/16)
            sage: w = w.augmentation(t^4 + (a^8 + a^12 + a^14 + a^16 + a^17 + a^19 + a^20 + a^23)*t^3 + (a^6 + a^9 + a^13 + a^15 + a^18 + a^19 + a^21)*t^2 + a^10*t + 1 + a^4 + a^5 + a^8 + a^13 + a^14 + a^15, 17/8)
            sage: f = a^-15*t^2 + (a^-11 + a^-9 + a^-6 + a^-5 + a^-3 + a^-2)*t + a^-15
            sage: f_ = w.equivalence_reciprocal(f)
            sage: w.reduce(f*f_)
            1
            sage: f = f.parent()([f[0], f[1].add_bigoh(1), f[2]])
            sage: f_ = w.equivalence_reciprocal(f)
            sage: w.reduce(f*f_)
            1

        """
        f = self.domain().coerce(f)

        from sage.categories.fields import Fields
        if not self.domain().base_ring() in Fields():
            # the xgcd does in general not work, i.e., return 1, unless over a field
            raise NotImplementedError("only implemented for polynomial rings over fields")

        if check:
            if coefficients is None:
                coefficients = list(self.coefficients(f))
            if valuations is None:
                valuations = list(self.valuations(f, coefficients=coefficients))
            if not self.is_equivalence_unit(f, valuations=valuations):
                raise ValueError("f must be an equivalence unit but %r is not"%(f,))

        if coefficients is None:
            e0 = self.coefficients(f).next()
        else:
            e0 = coefficients[0]

        # f is an equivalence unit, its valuation is given by the constant coefficient
        if valuations is None:
            vf = self(e0)
        else:
            vf = valuations[0]

        e0 = self.simplify(e0, error=vf)
        s_ = self.equivalence_unit(-vf)
        residue = self.reduce(e0 * s_)
        if not isinstance(self, FinalInductiveValuation):
            assert residue.is_constant()
            residue = residue[0]
        h = self.lift(~residue) * s_

        h = self.simplify(h, -vf)

        # it might be the case that f*h has non-zero valuation because h has
        # insufficient precision, so we must not assert that here but only
        # until we lifted to higher precision

        # We do not actually need g*phi + h*e0 = 1, it is only important that
        # the RHS is 1 in reduction.
        # This allows us to do two things:
        # - we may lift h to arbitrary precision
        # - we can add anything which times e0 has positive valuation, e.g., we
        # may drop coefficients of positive valuation
        h = h.map_coefficients(lambda c:_lift_to_maximal_precision(c))

        return h

    @cached_method
    def mu(self):
        r"""
        Return the valuation of :meth:`phi`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: v.mu()
            0

        """
        return self(self.phi())

    @abstract_method
    def equivalence_unit(self, s, reciprocal=False):
        """
        Return an equivalence unit of valuation ``s``.

        INPUT:

        - ``s`` -- an element of the :meth:`value_group`

        - ``reciprocal`` -- a boolean (default: ``False``); whether or not to
          return the equivalence unit as the :meth:`equivalence_reciprocal` of
          the equivalence unit of valuation ``-s``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: S.<x> = Qp(3,5)[]
            sage: v = GaussValuation(S)
            sage: v.equivalence_unit(2)
            (3^2 + O(3^7))
            sage: v.equivalence_unit(-2)
            (3^-2 + O(3^3))

        Note that this might fail for negative ``s`` if the domain is not
        defined over a field::

            sage: v = pAdicValuation(ZZ, 2)
            sage: R.<x> = ZZ[]
            sage: w = GaussValuation(R, v)
            sage: w.equivalence_unit(1)
            2
            sage: w.equivalence_unit(-1)
            Traceback (most recent call last):
            ...
            ValueError: s must be in the value semigroup of this valuation but -1 is not in Additive Abelian Semigroup generated by 1

        """

    @abstract_method
    def augmentation_chain(self):
        r"""
        Return a list with the chain of augmentations down to the underlying
        :class:`GaussValuation`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.augmentation_chain()
            [Gauss valuation induced by 2-adic valuation]

        """

    @abstract_method
    def is_gauss_valuation(self):
        r"""
        Return whether this valuation is a Gauss valuation over the domain.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_gauss_valuation()
            True

        """

    @abstract_method
    def E(self):
        """
        Return the ramification index of this valuation over its underlying
        Gauss valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.E()
            1

        """

    @abstract_method
    def F(self):
        """
        Return the residual degree of this valuation over its Gauss extension.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.F()
            1

        """

    @abstract_method
    def monic_integral_model(self, G):
        r"""
        Return a monic integral irreducible polynomial which defines the same
        extension of the base ring of the domain as the irreducible polynomial
        ``G`` together with maps between the old and the new polynomial.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: v.monic_integral_model(5*x^2 + 1/2*x + 1/4)
            (Ring endomorphism of Univariate Polynomial Ring in x over Rational Field
               Defn: x |--> 1/2*x,
             Ring endomorphism of Univariate Polynomial Ring in x over Rational Field
               Defn: x |--> 2*x,
            x^2 + 1/5*x + 1/5)

        """

    @abstract_method
    def element_with_valuation(self, s):
        r"""
        Return a polynomial of minimal degree with valuation ``s``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: v.element_with_valuation(-2)
            1/4

        Depending on the base ring, an element of valuation ``s`` might not
        exist::

            sage: R.<x> = ZZ[]
            sage: v = GaussValuation(R, pAdicValuation(ZZ, 2))
            sage: v.element_with_valuation(-2)
            Traceback (most recent call last):
            ...
            ValueError: s must be in the value semigroup of this valuation but -2 is not in Additive Abelian Semigroup generated by 1

        """

    def _test_element_with_valuation_inductive_valuation(self, **options):
        r"""
        Test the correctness of :meth:`element_with_valuation`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: v._test_element_with_valuation_inductive_valuation()

        """
        tester = self._tester(**options)
        chain = self.augmentation_chain()
        for s in tester.some_elements(self.value_group().some_elements()):
            try:
                R = self.element_with_valuation(s)
            except (ValueError, NotImplementedError):
                # this is often not possible unless the underlying ring of
                # constants is a field
                from sage.categories.fields import Fields
                if self.domain().base() not in Fields():
                    continue
                raise
            tester.assertEqual(self(R), s)
            if chain != [self]:
                base = chain[1]
                if s in base.value_group():
                    S = base.element_with_valuation(s)
                    tester.assertEqual(self(S), s)
                    tester.assertGreaterEqual(S.degree(), R.degree())

    def _test_EF(self, **options):
        r"""
        Test the correctness of :meth:`E` and :meth:`F`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v._test_EF()

        """
        tester = self._tester(**options)
        chain = self.augmentation_chain()
        for w,v in zip(chain, chain[1:]):
            from sage.rings.all import infinity, ZZ
            if w(w.phi()) is infinity:
                tester.assertEqual(w.E(), v.E())
            tester.assertIn(w.E(), ZZ)
            tester.assertIn(w.F(), ZZ)
            tester.assertGreaterEqual(w.E(), v.E())
            tester.assertGreaterEqual(w.F(), v.F())

    def _test_augmentation_chain(self, **options):
        r"""
        Test the correctness of :meth:`augmentation_chain`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_augmentation_chain()

        """
        tester = self._tester(**options)
        chain = self.augmentation_chain()
        tester.assertIs(chain[0], self)
        tester.assertTrue(chain[-1].is_gauss_valuation())
        for w,v in zip(chain, chain[1:]):
            tester.assertGreaterEqual(w, v)

    def _test_equivalence_unit(self, **options):
        r"""
        Test the correctness of :meth:`lift_to_key`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_equivalence_unit()

        """
        tester = self._tester(**options)

        if self.is_gauss_valuation():
            value_group = self.value_group()
        else:
            value_group = self.augmentation_chain()[1].value_group()

        for s in tester.some_elements(value_group.some_elements()):
            try:
                R = self.equivalence_unit(s)
            except (ValueError, NotImplementedError):
                # this is often not possible unless the underlying ring of
                # constants is a field
                from sage.categories.fields import Fields
                if self.domain().base() not in Fields():
                    continue
                raise
            tester.assertIs(R.parent(), self.domain())
            tester.assertEqual(self(R), s)
            tester.assertTrue(self.is_equivalence_unit(R))

    def _test_is_equivalence_unit(self, **options):
        r"""
        Test the correctness of :meth:`is_equivalence_unit`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_is_equivalence_unit()

        """
        tester = self._tester(**options)
        tester.assertFalse(self.is_equivalence_unit(self.phi()))

    def _test_equivalence_reciprocal(self, **options):
        r"""
        Test the correctness of :meth:`equivalence_reciprocal`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_equivalence_reciprocal()

        """
        tester = self._tester(**options)
        S = tester.some_elements(self.domain().some_elements())
        for f in S:
            if self.is_equivalence_unit(f):
                try:
                    g = self.equivalence_reciprocal(f)
                except (ValueError, NotImplementedError):
                    # this is often not possible unless the underlying ring of
                    # constants is a field
                    from sage.categories.fields import Fields
                    if self.domain().base() not in Fields():
                        continue
                    raise
                tester.assertEqual(self.reduce(f*g), 1)

    def _test_inductive_valuation_inheritance(self, **options):
        r"""
        Test that every instance that is a :class:`InductiveValuation` is
        either a :class:`FiniteInductiveValuation` or a
        :class:`InfiniteInductiveValuation`. Same for
        :class:`FinalInductiveValuation` and
        :class:`NonFinalInductiveValuation`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_inductive_valuation_inheritance()

        """
        tester = self._tester(**options)
        tester.assertTrue(isinstance(self, InfiniteInductiveValuation) != isinstance(self, FiniteInductiveValuation))
        tester.assertTrue(isinstance(self, FinalInductiveValuation) != isinstance(self, NonFinalInductiveValuation))


class FiniteInductiveValuation(InductiveValuation, DiscreteValuation):
    r"""
    Abstract base class for iterated :class:`AugmentedValuation` on top of a
    :class:`GaussValuation` which is a discrete valuation, i.e., the last key
    polynomial has finite valuation.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, TrivialValuation(QQ))

    """
    def __init__(self, parent, phi):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: isinstance(v, FiniteInductiveValuation)
            True

        """
        InductiveValuation.__init__(self, parent, phi)
        DiscreteValuation.__init__(self, parent)

    def extensions(self, other):
        r"""
        Return the extensions of this valuation to ``other``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = ZZ[]
            sage: v = GaussValuation(R, TrivialValuation(ZZ))
            sage: K.<x> = FunctionField(QQ)
            sage: v.extensions(K)
            [Trivial valuation on Rational Field]

        """
        from sage.categories.function_fields import FunctionFields
        if other in FunctionFields() and other.ngens() == 1:
            # extend to K[x] and from there to K(x)
            v = self.extension(self.domain().change_ring(self.domain().base().fraction_field()))
            from function_field_valuation import FunctionFieldValuation
            return [FunctionFieldValuation(other, v)]
        return super(FiniteInductiveValuation, self).extensions(other)


class NonFinalInductiveValuation(FiniteInductiveValuation, DiscreteValuation):
    r"""
    Abstract base class for iterated :class:`AugmentedValuation` on top of a
    :class:`GaussValuation` which can be extended further through
    :meth:`augmentation`.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: R.<u> = Qq(4,5)
        sage: S.<x> = R[]
        sage: v = GaussValuation(S)
        sage: v = v.augmentation(x^2 + x + u, 1)

    """
    def __init__(self, parent, phi):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.augmentation(x^2 + x + u, 1)
            sage: isinstance(v, NonFinalInductiveValuation)
            True

        """
        FiniteInductiveValuation.__init__(self, parent, phi)
        DiscreteValuation.__init__(self, parent)

    def augmentation(self, phi, mu, check=True):
        r"""
        Return the inductive valuation which extends this valuation by mapping
        ``phi`` to ``mu``.

        INPUT:

        - ``phi`` -- a polynomial in the domain of this valuation; this must be
          a key polynomial, see :meth:`is_key` for properties of key
          polynomials.

        - ``mu`` -- a rational number or infinity, the valuation of ``phi`` in
          the extended valuation

        - ``check`` -- a boolean (default: ``True``), whether or not to check
          the correctness of the parameters

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.augmentation(x^2 + x + u, 1)
            sage: v = v.augmentation((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
            sage: v
            [ Gauss valuation induced by 2-adic valuation,
              v((1 + O(2^5))*x^2 + (1 + O(2^5))*x + u + O(2^5)) = 1,
              v((1 + O(2^5))*x^4 + (2^2 + O(2^6))*x^3 + (1 + (u + 1)*2 + O(2^5))*x^2 + ((u + 1)*2^2 + O(2^6))*x + (u + 1) + (u + 1)*2 + (u + 1)*2^2 + (u + 1)*2^3 + (u + 1)*2^4 + O(2^5)) = 3 ]

        TESTS:

        Make sure that we do not make the assumption that the degrees of the
        key polynomials are strictly increasing::

            sage: v_K = pAdicValuation(QQ,3)
            sage: A.<t> = QQ[]
            sage: v0 = GaussValuation(A,v_K)

            sage: v1 = v0.augmentation(t, 1/12)
            sage: v2 = v1.augmentation(t^12 + 3, 7/6)
            sage: v3 = v2.augmentation(t^12 + 3*t^2 + 3, 9/4)
            sage: v4 = v1.augmentation(t^12 + 3*t^2 + 3, 9/4)
            sage: v3 <= v4 and v3 >= v4
            True

        .. SEEALSO::

            :meth:`AugmentedValuation`

        """
        from augmented_valuation import AugmentedValuation
        return AugmentedValuation(self, phi, mu, check)

    def mac_lane_step(self, G, principal_part_bound=None, assume_squarefree=False, assume_equivalence_irreducible=False, report_degree_bounds_and_caches=False, coefficients=None, valuations=None, check=True):
        r"""
        Perform an approximation step towards the squarefree monic non-constant
        integral polynomial ``G`` which is not an :meth:`equivalence_unit`.

        This performs the individual steps that are used in
        :meth:`mac_lane_approximants`.

        INPUT:

        - ``G`` -- a sqaurefree monic non-constant integral polynomial ``G``
          which is not an :meth:`equivalence_unit`

        - ``principal_part_bound`` -- an integer or ``None`` (default:
          ``None``), a bound on the length of the principal part, i.e., the
          section of negative slope, of the Newton polygon of ``G``

        - ``assume_squarefree`` -- whether or not to assume that ``G`` is
          squarefree (default: ``False``)

        - ``assume_equivalence_irreducible`` -- whether or not to assume that
          ``G`` is equivalence irreducible (default: ``False``)

        - ``report_degree_bounds_and_caches`` -- whether or not to include internal state with the returned value (used by :meth:`mac_lane_approximants` to speed up sequential calls)

         - ``coefficients`` -- the coefficients of ``G`` in the
           :meth:`phi`-adic expansion if known (default: ``None``)

        - ``valauations`` -- the valuations of ``coefficients`` if known
          (default: ``None``)

        - ``check`` -- whether to check that ``G`` is a squarefree monic
          non-constant  integral polynomial and not an :meth:`equivalence_unit`
          (default: ``True``)

        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: K.<x> = FunctionField(QQ)
            sage: S.<y> = K[]
            sage: F = y^2 - x^2 - x^3 - 3
            sage: v0 = GaussValuation(K._ring,pAdicValuation(QQ, 3))
            sage: v1 = v0.augmentation(K._ring.gen(), 1/3)
            sage: mu0 = FunctionFieldValuation(K, v1)
            sage: eta0 = GaussValuation(S, mu0)
            sage: eta1 = eta0.mac_lane_step(F)[0]
            sage: eta2 = eta1.mac_lane_step(F)[0]
            sage: eta2
            [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by 3-adic valuation, v(x) = 1/3 ], v(y + x) = 2/3 ]

        """
        G = self.domain().coerce(G)

        if G.is_constant():
            raise ValueError("G must not be constant")

        from itertools import islice
        from sage.misc.misc import verbose
        verbose("Augmenting %s towards %s"%(self, G), level=10)

        if not G.is_monic():
            raise ValueError("G must be monic")

        if coefficients is None:
            coefficients = self.coefficients(G)
            if principal_part_bound:
                coefficients = islice(coefficients, 0, principal_part_bound + 1, 1)
            coefficients = list(coefficients)
        if valuations is None:
            valuations = self.valuations(G, coefficients=coefficients)
            if principal_part_bound:
                valuations = islice(valuations, 0, principal_part_bound + 1, 1)
            valuations = list(valuations)

        if check and min(valuations) < 0:
            raise ValueError("G must be integral")

        if check and self.is_equivalence_unit(G, valuations=valuations):
            raise ValueError("G must not be an equivalence-unit")

        if check and not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.all import infinity
        assert self(G) is not infinity # this is a valuation and G is non-zero

        ret = []

        F = self.equivalence_decomposition(G, assume_not_equivalence_unit=True, coefficients=coefficients, valuations=valuations, compute_unit=False, degree_bound=principal_part_bound)
        assert len(F), "%s equivalence-decomposes as an equivalence-unit %s"%(G, F)
        if len(F) == 1 and F[0][1] == 1 and F[0][0].degree() == G.degree():
            assert self.is_key(G, assume_equivalence_irreducible=assume_equivalence_irreducible)
            ret.append((self.augmentation(G, infinity, check=False), G.degree(), principal_part_bound, None, None))
        else:
            for phi,e in F:
                if G == phi:
                    # Something strange happened here:
                    # G is not a key (we checked that before) but phi==G is; so phi must have less precision than G
                    # this can happen if not all coefficients of G have the same precision
                    # if we drop some precision of G then it will be a key (but is
                    # that really what we should do?)
                    assert not G.base_ring().is_exact()
                    prec = min([c.precision_absolute() for c in phi.list()])
                    g = G.map_coefficients(lambda c:c.add_bigoh(prec))
                    assert self.is_key(g)
                    ret.append((self.augmentation(g, infinity, check=False), g.degree(), principal_part_bound, None, None))
                    assert len(F) == 1
                    break

                if phi == self.phi():
                    # a factor phi in the equivalence decomposition means that we
                    # found an actual factor of G, i.e., we can set
                    # v(phi)=infinity
                    # However, this should already have happened in the last step
                    # (when this polynomial had -infinite slope in the Newton
                    # polygon.)
                    if self.is_gauss_valuation(): # unless in the first step
                        pass
                    else:
                        continue

                verbose("Determining the augmentation of %s for %s"%(self, phi), level=11)
                old_mu = self(phi)
                w = self.augmentation(phi, old_mu, check=False)

                # we made some experiments here: instead of computing the
                # coefficients again from scratch, update the coefficients when
                # phi - self.phi() is a constant.
                # It turned out to be slightly slower than just recomputing the
                # coefficients. The main issue with the approach was that we
                # needed to keep track of all the coefficients and not just of
                # the coefficients up to principal_part_bound.

                w_coefficients = w.coefficients(G)
                if principal_part_bound:
                    w_coefficients = islice(w_coefficients, 0, principal_part_bound + 1, 1)
                w_coefficients = list(w_coefficients)

                w_valuations = w.valuations(G, coefficients=w_coefficients)
                if principal_part_bound:
                    w_valuations = islice(w_valuations, 0, principal_part_bound + 1, 1)
                w_valuations = list(w_valuations)

                NP = w.newton_polygon(G, valuations=w_valuations).principal_part()

                verbose("Newton-Polygon for v(phi)=%s : %s"%(self(phi), NP), level=11)
                slopes = NP.slopes(repetition=True)
                multiplicities = {slope : len([s for s in slopes if s == slope]) for slope in slopes}
                slopes = multiplicities.keys()
                if NP.vertices()[0][0] != 0:
                    slopes = [-infinity] + slopes
                    multiplicities[-infinity] = 1

                if not slopes:
                    q,r = G.quo_rem(phi)
                    assert not r.is_zero()
                    phi = phi.coefficients(sparse=False)
                    for i,c in enumerate(r.coefficients(sparse=False)):
                        if not c.is_zero():
                            v = w(c)
                            # for a correct result we need to add O(pi^v) in degree i
                            # we try to find the coefficient of phi where such an
                            # error can be introduced without losing much absolute
                            # precision on phi
                            best = i
                            for j in range(i):
                                if w(q[j]) < w(q[best]):
                                    best = j
                            # now add the right O() to phi in degree i - best
                            phi[i-best] = phi[i-best].add_bigoh(w(c)-w(q[best]))

                    phi = G.parent()(phi)
                    w = self._base_valuation.augmentation(phi, infinity, check=False)
                    ret.append((w, phi.degree(), principal_part_bound, None, None))
                    continue

                for i, slope in enumerate(slopes):
                    slope = slopes[i]
                    verbose("Slope = %s"%slope, level=12)
                    new_mu = old_mu - slope
                    new_valuations = [val - (j*slope if slope is not -infinity else (0 if j == 0 else -infinity)) for j,val in enumerate(w_valuations)]
                    base = self
                    if phi.degree() == base.phi().degree():
                        assert new_mu > self(phi)
                        if not base.is_gauss_valuation():
                            base = base._base_valuation
                    w = base.augmentation(phi, new_mu, check=False)
                    assert slope is -infinity or 0 in w.newton_polygon(G).slopes(repetition=False)

                    from sage.rings.all import ZZ
                    assert (phi.degree() / self.phi().degree()) in ZZ
                    degree_bound = multiplicities[slope] * phi.degree()
                    assert degree_bound <= G.degree()
                    assert degree_bound >= phi.degree()
                    ret.append((w, degree_bound, multiplicities[slope], w_coefficients, new_valuations))

        assert ret
        if not report_degree_bounds_and_caches:
            ret = [v for v,_,_,_,_ in ret]
        return ret

    def is_key(self, phi, explain=False, assume_equivalence_irreducible=False):
        r"""
        Return whether ``phi`` is a key polynomial for this valuation, i.e.,
        whether it is monic, whether it :meth:`is_equivalence_irreducible`, and
        whether it is :meth:`is_minimal`.

        INPUT:

        - ``phi`` -- a polynomial in the domain of this valuation

        - ``explain`` -- a boolean (default: ``False``), if ``True``, return a
          string explaining why ``phi`` is not a key polynomial

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_key(x)
            True
            sage: v.is_key(2*x, explain = True)
            (False, 'phi must be monic')
            sage: v.is_key(x^2, explain = True)
            (False, 'phi must be equivalence irreducible')

            sage: w = v.augmentation(x, 1)
            sage: w.is_key(x + 1, explain = True)
            (False, 'phi must be minimal')

        """
        phi = self.domain().coerce(phi)

        reason = None

        if not phi.is_monic():
            reason = "phi must be monic"
        elif not assume_equivalence_irreducible and not self.is_equivalence_irreducible(phi):
            reason = "phi must be equivalence irreducible"
        elif not self.is_minimal(phi, assume_equivalence_irreducible=True):
            reason = "phi must be minimal"

        if explain:
            return reason is None, reason
        else:
            return reason is None

    def is_minimal(self, f, assume_equivalence_irreducible=False):
        r"""
        Return whether the polynomial ``f`` is minimal with respect to this
        valuation, i.e., whether ``f`` is not constant any non-constant
        polynomial `h` has at least the degree of ``f`` or ``f`` is not
        divisible by `h` with respect to this valuation, i.e., there is no `c`
        such that `c h` :meth:`is_equivalent` to `f`.

        ALGORITHM:

        Based on Theorem 9.4 of [ML1936'].

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_minimal(x + 1)
            True
            sage: w = v.augmentation(x, 1)
            sage: w.is_minimal(x + 1)
            False

        TESTS::

            sage: K = Qp(2, 10)
            sage: R.<x> = K[]
            sage: vp = pAdicValuation(K)
            sage: v0 = GaussValuation(R, vp)
            sage: v1 = v0.augmentation(x, 1/4)
            sage: v2 = v1.augmentation(x^4 + 2, 5/4)
            sage: v2.is_minimal(x^5 + x^4 + 2)
            False

        Polynomials which are equivalent to the key polynomial are minimal if
        and only if they have the same degree as the key polynomial::

            sage: v2.is_minimal(x^4 + 2)
            True
            sage: v2.is_minimal(x^4 + 4)
            False

        """
        f = self.domain().coerce(f)

        if f.is_constant():
            return False

        if not assume_equivalence_irreducible and not self.is_equivalence_irreducible(f):
            # any factor divides f with respect to this valuation
            return False

        if not f.is_monic():
            # divide out the leading factor, it does not change minimality
            v = self
            if not self.domain().base_ring().is_field():
                domain = self.domain().change_ring(self.domain().base_ring().fraction_field())
                v = self.extension(domain)
                f = domain(f)
            return v.is_minimal(f / f.leading_coefficient())

        if self.is_gauss_valuation():
            if self(f) == 0:
                F = self.reduce(f, check=False)
                assert not F.is_constant()
                return F.is_irreducible()
            else:
                assert(self(f) <= 0) # f is monic
                # f is not minimal:
                # Let g be f stripped of its leading term, i.e., g = f - x^n.
                # Then g and f are equivalent with respect to this valuation
                # and in particular g divides f with respect to this valuation
                return False

        if self.is_equivalent(self.phi(), f):
            assert f.degree() >= self.phi().degree()
            # If an h divides f with respect to this valuation, then it also divides phi:
            # v(f - c*h) > v(f) = v(c*h) => v(phi - c*h) = v((phi - f) + (f - c*h)) > v(phi) = v(c*h)
            # So if f were not minimal then phi would not be minimal but it is.
            return f.degree() == self.phi().degree()

        else:
            tau = self.value_group().index(self._base_valuation.value_group())
            # see Theorem 9.4 of [ML1936']
            return list(self.valuations(f))[-1] == self(f) and \
                   list(self.coefficients(f))[-1].is_constant() and \
                   list(self.valuations(f))[0] == self(f) and \
                   tau.divides(len(list(self.coefficients(f))) - 1)

    def _equivalence_reduction(self, f, coefficients=None, valuations=None, degree_bound=None):
        r"""
        Helper method for :meth:`is_equivalence_irreducible` and
        :meth:`equivalence_decomposition` which essentially returns the
        reduction of ``f`` after multiplication with an ``R`` which
        :meth:`is_equivalence_unit`.

        This only works when ``f`` is not divisible by :meth:`phi` with respect
        to this valuation. Therefore, we also return the number of times that
        we took out :meth:`phi` of ``f`` before we computed the reduction.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: v._equivalence_reduction(2*x^6 + 4*x^5 + 2*x^4 + 8)
            (1, 4, x^2 + 1)

        """
        f = self.domain().coerce(f)

        # base change from R[x] to K[x], so divisions work and sufficient
        # elements of negative valuation exist
        if not self.domain().base_ring().is_field():
            domain = self.domain().change_ring(self.domain().base_ring().fraction_field())
            v = self.extension(domain)
            assert self.residue_ring() is v.residue_ring()
            return v._equivalence_reduction(f)

        if coefficients is None:
            coefficients = list(self.coefficients(f))
        if valuations is None:
            valuations = list(self.valuations(f, coefficients=coefficients))
        valuation = min(valuations)
        for phi_divides in range(len(valuations)):
            # count how many times phi divides f
            if valuations[phi_divides] <= valuation:
                break

        if phi_divides:
            from sage.rings.all import PolynomialRing
            R = PolynomialRing(f.parent(), 'phi')
            f = R(coefficients[phi_divides:])(self.phi())
        valuations = [v-self.mu()*phi_divides for v in valuations[phi_divides:]]
        coefficients = coefficients[phi_divides:]
        valuation = min(valuations)

        R = self.equivalence_unit(-valuation)
        R = self.coefficients(R).next()
        fR_valuations = [v-valuation for v in valuations]
        from sage.rings.all import infinity
        fR_coefficients = [self.coefficients(c*R).next() if v is not infinity and v == 0 else 0 for c,v in zip(coefficients,fR_valuations)]

        return valuation, phi_divides, self.reduce(f*R, check=False, degree_bound=degree_bound, coefficients=fR_coefficients, valuations=fR_valuations)

    def is_equivalence_irreducible(self, f, coefficients=None, valuations=None):
        r"""
        Return whether the polynomial ``f`` is equivalence-irreducible, i.e.,
        whether its :meth:`equivalence_decomposition` is trivial.

        ALGORITHM:

        We use the same algorithm as in :meth:`equivalence_decomposition` we
        just do not lift the result to key polynomials.

        INPUT:

        - ``f`` -- a non-constant polynomial in the domain of this valuation

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_equivalence_irreducible(x)
            True
            sage: v.is_equivalence_irreducible(x^2)
            False
            sage: v.is_equivalence_irreducible(x^2 + 2)
            False

        """
        f = self.domain().coerce(f)

        if not self.domain().base_ring().is_field():
            domain = self.domain().change_ring(self.domain().base_ring().fraction_field())
            v = self.extension(domain)
            return v.is_equivalence_irreducible(v.domain()(f))

        if f.is_constant():
            raise ValueError("f must not be constant")

        _, phi_divides, F = self._equivalence_reduction(f, coefficients=coefficients, valuations=valuations)
        if phi_divides == 0:
            return F.is_constant() or F.is_irreducible()
        if phi_divides == 1:
            return F.is_constant()
        if phi_divides > 1:
            return False

    def equivalence_decomposition(self, f, assume_not_equivalence_unit=False, coefficients=None, valuations=None, compute_unit=True, degree_bound=None):
        r"""
        Return an equivalence decomposition of ``f``, i.e., a polynomial
        `g(x)=e(x)\prod_i \phi_i(x)` with `e(x)` an equivalence unit (see
        :meth:`is_equivalence_unit()`) and the `\phi_i` key polynomials (see
        :meth:`is_key`) such that ``f`` :meth:`is_equivalent` to `g`.

        INPUT:

        - ``f`` -- a non-zero polynomial in the domain of this valuation

        - ``assume_not_equivalence_unit`` -- whether or not to assume that
          ``f`` is not an :meth:`equivalence_unit` (default: ``False``)

        - ``coefficients`` -- the coefficients of ``f`` in the :meth:`phi`-adic
          expansion if known (default: ``None``)

        - ``valuations`` -- the valuations of ``coefficients`` if known
          (default: ``None``)

        - ``compute_unit`` -- whether or not to compute the unit part of the
          decomposition (default: ``True``)

        - ``degree_bound`` -- a bound on the degree of the
          :meth:`_equivalence_reduction` of ``f`` (default: ``None``)

        ALGORITHM:

        We use the algorithm described in Theorem 4.4 of [ML1936']. After
        removing all factors `\phi` from a polynomial `f`, there is an
        equivalence unit `R` such that `Rf` has valuation zero. Now `Rf` can be
        factored as `\prod_i \alpha_i` over the :meth:`residue_field`. Lifting
        all `\alpha_i` to key polynomials `\phi_i` gives `Rf=\prod_i R_i f_i`
        for suitable equivalence units `R_i` (see :meth:`lift_to_key`). Taking
        `R'` an :meth:`equivalence_reciprocal` of `R`, we have `f` equivalent
        to `(R'\prod_i R_i)\prod_i\phi_i`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.equivalence_decomposition(S.zero())
            Traceback (most recent call last):
            ...
            ValueError: equivalence decomposition of zero is not defined
            sage: v.equivalence_decomposition(S.one())
            1 + O(2^10)
            sage: v.equivalence_decomposition(x^2+2)
            ((1 + O(2^10))*x)^2
            sage: v.equivalence_decomposition(x^2+1)
            ((1 + O(2^10))*x + 1 + O(2^10))^2

        A polynomial that is an equivalence unit, is returned as the unit part
        of a :class:`sage.structure.factorization.Factorization`, leading to a unit
        non-minimal degree::

            sage: w = v.augmentation(x, 1)
            sage: F = w.equivalence_decomposition(x^2+1); F
            (1 + O(2^10))*x^2 + 1 + O(2^10)
            sage: F.unit()
            (1 + O(2^10))*x^2 + 1 + O(2^10)

        However, if the polynomial has a non-unit factor, then the unit might
        be replaced by a factor of lower degree::

            sage: f = x * (x^2 + 1)
            sage: F = w.equivalence_decomposition(f); F
            (1 + O(2^10))*x
            sage: F.unit()
            1 + O(2^10)

        Examples over an iterated unramified extension:

            sage: v = v.augmentation(x^2 + x + u, 1)
            sage: v = v.augmentation((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)

            sage: v.equivalence_decomposition(x)
            (1 + O(2^10))*x
            sage: F = v.equivalence_decomposition( v.phi() )
            sage: len(F)
            1
            sage: F = v.equivalence_decomposition( v.phi() * (x^4 + 4*x^3 + (7 + 2*u)*x^2 + (8 + 4*u)*x + 1023 + 3*u) )
            sage: len(F)
            2

        TESTS::

            sage: R.<x> = QQ[]
            sage: K1.<pi>=NumberField(x^3-2)
            sage: K.<alpha>=K1.galois_closure()
            sage: R.<x>=K[]
            sage: vp=pAdicValuation(QQ,2)
            sage: vp=vp.extension(K)
            sage: v0=GaussValuation(R,vp)
            sage: G=x^36 + 36*x^35 + 630*x^34 + 7144*x^33 + 59055*x^32 + 379688*x^31 +1978792*x^30 + 8604440*x^29 + 31895428*x^28 + 102487784*x^27 + 289310720*x^26 + 725361352*x^25 + 1629938380*x^24 + 3307417800*x^23 + 6098786184*x^22+10273444280*x^21 + 15878121214*x^20 + 22596599536*x^19 + 29695703772*x^18 +36117601976*x^17 + 40722105266*x^16 + 42608585080*x^15 + 41395961848*x^14 +37344435656*x^13 + 31267160756*x^12 + 24271543640*x^11 + 17439809008*x^10 + 11571651608*x^9 + 7066815164*x^8 + 3953912472*x^7 + 2013737432*x^6 + 925014888*x^5 + 378067657*x^4 + 134716588*x^3 + 40441790*x^2 + 9532544*x + 1584151
            sage: v1=v0.mac_lane_step(G)[0]
            sage: V=v1.mac_lane_step(G)
            sage: v2=V[0]
            sage: v2.equivalence_decomposition(G)
            (1/387420489) * (x^4 + 2*x^2 + alpha^4 + alpha^3 + 1)^3 * (x^4 + 2*x^2 + 1/2*alpha^4 + alpha^3 + 5*alpha + 1)^3 * (x^4 + 2*x^2 + 3/2*alpha^4 + alpha^3 + 5*alpha + 1)^3

        REFERENCES:

        .. [ML1936'] MacLane, S. (1936). A construction for absolute values in
        polynomial rings. Transactions of the American Mathematical Society, 40(3),
        363-395.

        """
        f = self.domain().coerce(f)

        if f.is_zero():
            raise ValueError("equivalence decomposition of zero is not defined")

        from sage.structure.factorization import Factorization
        if not assume_not_equivalence_unit and self.is_equivalence_unit(f):
            return Factorization([], unit=f, sort=False)

        if not self.domain().base_ring().is_field():
            domain = self.domain().change_ring(self.domain().base_ring().fraction_field())
            v = self.extension(domain)
            ret = v.equivalence_decomposition(v.domain()(f))
            return Factorization([(g.change_ring(self.domain().base_ring()),e) for g,e in ret], unit=ret.unit().change_ring(self.domain().base_ring()), sort=False)

        valuation, phi_divides, F = self._equivalence_reduction(f, coefficients=coefficients, valuations=valuations, degree_bound=degree_bound)

        # fixing broken factorization for polynomials over fraction fields
        # original code:
        # F = F.factor()
        # the following code has been added by S.Wewers on Sept. 10, 2017
        # if the coefficients of f lie in a quotient ring of a polynomial ring
        # we change into elements of a rational function field, factor F
        # and reverse the procedure for the factors
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        from sage.rings.function_field.constructor import FunctionField
        K = F.base_ring()
        if hasattr(K, "base") and not hasattr(K, "modulus") and is_PolynomialRing(K.base()):
            K1 = FunctionField(K.base_ring(), K.variable_name())
            factors = F.change_ring(K1).factor()
            to_K = lambda c: K(c.numerator())/K(c.denominator())
            F = Factorization([(g.map_coefficients(to_K, K), m) for g, m in factors if not g.is_constant()], to_K(K1(factors.unit())), sort=False)
            # the if clause is strange but necessary
        else:
            F = F.factor()

        from sage.misc.misc import verbose
        verbose("%s factors as %s = %s in reduction"%(f, F.prod(), F), level=20)

        unit = self.domain().one()
        if compute_unit:
            R_ = self.equivalence_unit(valuation, reciprocal=True)
            unit = self.lift(self.residue_ring()(F.unit())) * R_
        F = list(F)

        if compute_unit:
            from sage.misc.all import prod
            unit *= self.lift(self.residue_ring()(prod([ psi.leading_coefficient()**e for psi,e in F ])))

        # A potential speedup that we tried to implement here:
        # When F factors as T^n - a, then instead of using any lift of T^n - a
        # we tried to take a lift that approximates well an n-th root of the
        # constant coefficient of f[0]. Doing so saved a few invocations of
        # mac_lane_step but in the end made hardly any difference.

        F = [(self.lift_to_key(psi/psi.leading_coefficient()),e) for psi,e in F]

        if compute_unit:
            for g,e in F:
                v_g = self(g)
                unit *= self._pow(self.equivalence_unit(-v_g, reciprocal=True), e, error=-v_g*e, effective_degree=0)
            unit = self.simplify(unit)

        if phi_divides:
            for i,(g,e) in enumerate(F):
                if g == self.phi():
                    F[i] = (self.phi(),e+phi_divides)
                    break
            else:
                F.append((self.phi(),phi_divides))

        ret = Factorization(F, unit=unit, sort=False)

        if compute_unit:
            assert self.is_equivalent(ret.prod(), f) # this might fail because of leading zeros in inexact rings
            assert self.is_equivalence_unit(ret.unit())

        return ret

    def minimal_representative(self, f):
        r"""
        Return a minimal representative for ``f``, i.e., a pair `e, a` such
        that ``f`` :meth:`is_equivalent` to `e a`, `e` is an
        :meth:`equivalence_unit`, and `a` :meth:`is_minimal` and monic.

        INPUT:

        - ``f`` -- a non-zero polynomial which is not an equivalence unit

        OUTPUT:

        A factorization which has `e` as its unit and `a` as its unique factor.

        ALGORITHM:

        We use the algorithm described in the proof of Lemma 4.1 of [ML1936'].
        In the expansion `f=\sum_i f_i\phi^i` take `e=f_i` for the largest `i`
        with `f_i\phi^i` minimal (see :meth:`effective_degree`).
        Let `h` be the :meth:`equivalence_reciprocal` of `e` and take `a` given
        by the terms of minimal valuation in the expansion of `e f`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.minimal_representative(x + 2)
            (1 + O(2^10))*x

            sage: v = v.augmentation(x, 1)
            sage: v.minimal_representative(x + 2)
            (1 + O(2^10))*x + 2 + O(2^11)
            sage: f = x^3 + 6*x + 4
            sage: F = v.minimal_representative(f); F
            (2 + 2^2 + O(2^11)) * ((1 + O(2^10))*x + 2 + O(2^11))
            sage: v.is_minimal(F[0][0])
            True
            sage: v.is_equivalent(F.prod(), f)
            True

        """
        f = self.domain().coerce(f)

        from sage.categories.fields import Fields
        if not self.domain().base_ring() in Fields():
            raise NotImplementedError("only implemented for polynomial rings over fields")

        if f.is_zero():
            raise ValueError("zero has no minimal representative")

        degree = self.effective_degree(f)
        if degree == 0:
            raise ValueError("equivalence units can not have a minimal representative")

        e = list(self.coefficients(f))[degree]
        h = self.equivalence_reciprocal(e).map_coefficients(lambda c:_lift_to_maximal_precision(c))
        g = h*f
        vg = self(g)

        coeffs = [c if v == vg else c.parent().zero() for v,c in zip(self.valuations(g), self.coefficients(g))]
        coeffs[degree] = self.domain().base_ring().one()
        ret = sum([c*self._phi**i for i,c in enumerate(coeffs)])

        assert self.effective_degree(ret) == degree
        assert ret.is_monic()
        assert self.is_minimal(ret)

        from sage.structure.factorization import Factorization
        ret = Factorization([(ret, 1)], unit=e, sort=False)

        assert self.is_equivalent(ret.prod(), f) # this might fail because of leading zeros
        return ret

    @abstract_method
    def lift_to_key(self, F):
        """
        Lift the irreducible polynomial ``F`` from the :meth:`residue_ring` to
        a key polynomial over this valuation.

        INPUT:

        - ``F`` -- an irreducible non-constant monic polynomial in
          :meth:`residue_ring` of this valuation

        OUTPUT:

        A polynomial `f` in the domain of this valuation which is a key
        polynomial for this valuation and which is such that an
        :meth:`augmentation` with this polynomial adjoins a root of ``F`` to
        the resulting :meth:`residue_ring`.

        More specifically, if ``F`` is not the generator of the residue ring,
        then multiplying ``f`` with the :meth:`equivalence_reciprocal` of the
        :meth:`equivalence_unit` of the valuation of ``f``, produces a unit
        which reduces to ``F``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: y = v.residue_ring().gen()
            sage: u0 = v.residue_ring().base_ring().gen()
            sage: f = v.lift_to_key(y^2 + y + u0); f
            (1 + O(2^10))*x^2 + (1 + O(2^10))*x + u + O(2^10)

        """

    def _test_lift_to_key(self, **options):
        r"""
        Test the correctness of :meth:`lift_to_key`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_lift_to_key()

        """
        tester = self._tester(**options)

        try:
            k = self.residue_ring()
        except NotImplementedError:
            from sage.categories.fields import Fields
            if self.domain().base() in Fields():
                raise
            return

        S = tester.some_elements(self.residue_ring().some_elements())
        for F in S:
            if F.is_monic() and not F.is_constant() and F.is_irreducible():
                try:
                    f = self.lift_to_key(F)
                except NotImplementedError:
                    from sage.categories.fields import Fields
                    if self.domain().base() in Fields():
                        raise
                    continue
                tester.assertIs(f.parent(), self.domain())
                tester.assertTrue(self.is_key(f))

                # check that augmentation produces a valuation with roots of F
                # in the residue ring
                from sage.rings.all import infinity
                w = self.augmentation(f, infinity)
                F = F.change_ring(w.residue_ring())
                roots = F.roots(multiplicities=False)
                tester.assertGreaterEqual(len(roots), 1)

                # check that f has the right reduction
                if F == F.parent().gen():
                    tester.assertTrue(self.is_equivalent(f, self.phi()))
                else:
                    tester.assertEqual(self.reduce(f * self.equivalence_reciprocal(self.equivalence_unit(self(f)))), F)


    def _test_is_equivalence_irreducible(self, **options):
        r"""
        Test the correctness of :meth:`is_equivalence_irreducible`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, TrivialValuation(QQ))
            sage: v._test_is_equivalence_irreducible()

        """
        tester = self._tester(**options)
        S = tester.some_elements(self.domain().some_elements())
        for f in S:
            if f.is_constant(): continue
            is_equivalence_irreducible = self.is_equivalence_irreducible(f)
            F = self.equivalence_decomposition(f)
            tester.assertEqual(is_equivalence_irreducible, len(F)==0 or (len(F)==1 and F[0][1]==1))
            if self.is_equivalence_unit(f):
                tester.assertTrue(f.is_constant() or self.is_equivalence_irreducible(f))

        tester.assertTrue(self.is_equivalence_irreducible(self.phi()))
        tester.assertTrue(self.is_equivalence_irreducible(-self.phi()))
        tester.assertFalse(self.is_equivalence_irreducible(self.phi() ** 2))


class FinalInductiveValuation(InductiveValuation):
    r"""
    Abstract base class for an inductive valuation which can not be augmented further.

    TESTS::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, TrivialValuation(QQ))
        sage: w = v.augmentation(x^2 + x + 1, infinity)
        sage: isinstance(w, FinalInductiveValuation)
        True

    """


class InfiniteInductiveValuation(FinalInductiveValuation, InfiniteDiscretePseudoValuation):
    r"""
    Abstract base class for an inductive valuation which is not discrete, i.e.,
    which assigns infinite valuation to its last key polynomial.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
        sage: w = v.augmentation(x^2 + x + 1, infinity)

    """
    def __init__(self, parent, base_valuation):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x^2 + x + 1, infinity)
            sage: isinstance(w, InfiniteInductiveValuation)
            True

        """
        FinalInductiveValuation.__init__(self, parent, base_valuation)
        InfiniteDiscretePseudoValuation.__init__(self, parent)

    def change_domain(self, ring):
        r"""
        Return this valuation over ``ring``.

        EXAMPLES:

        We can turn an infinite valuation into a valuation on the quotient::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, pAdicValuation(QQ, 2))
            sage: w = v.augmentation(x^2 + x + 1, infinity)
            sage: w.change_domain(R.quo(x^2 + x + 1))
            2-adic valuation

        """
        from sage.rings.polynomial.polynomial_quotient_ring import is_PolynomialQuotientRing
        if is_PolynomialQuotientRing(ring) and ring.base() is self.domain() and ring.modulus() == self.phi():
            return self.restriction(self.domain().base())._extensions_to_quotient(ring, approximants=[self])[0]
        return super(InfiniteInductiveValuation, self).change_domain(ring)


def _lift_to_maximal_precision(c):
    r"""
    Lift ``c`` to maximal precision if the parent is not exact.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: from mac_lane.inductive_valuation import _lift_to_maximal_precision # optional: standalone
        sage: R = Zp(2,5)
        sage: x = R(1,2); x
        1 + O(2^2)
        sage: _lift_to_maximal_precision(x)
        1 + O(2^5)

        sage: x = 1
        sage: _lift_to_maximal_precision(x)
        1

    """
    return c if c.parent().is_exact() else c.lift_to_precision()
