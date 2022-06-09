from .shapes import polished_tetrahedra_shapes
from ..sage_helper import _within_sage, sage_method
from .polished_reps import polished_holonomy
from . import nsagetools, interval_reps, slice_obs_HKL, ManifoldNT
from .character_varieties import  character_variety, character_variety_ideal

if _within_sage:
    from .find_field import ListOfApproximateAlgebraicNumbers

@sage_method
def tetrahedra_field_gens(manifold):
    """
    The shapes of the tetrahedra as ApproximateAlgebraicNumbers. Can be
    used to compute the tetrahedra field, where the first two parameters
    are bits of precision and maximum degree of the field::

        sage: M = Manifold('m015')
        sage: tets = M.tetrahedra_field_gens()
        sage: tets.find_field(100, 10, optimize=True)    # doctest: +NORMALIZE_WHITESPACE +NUMERIC9
        (Number Field in z with defining polynomial x^3 - x - 1
         with z = -0.6623589786223730? - 0.5622795120623013?*I,
        <ApproxAN: -0.662358978622 - 0.562279512062*I>, [-z, -z, -z])
    """
    if manifold.is_orientable():
        def func(prec):
            return polished_tetrahedra_shapes(manifold, bits_prec=prec)
    else:
        double_cover = manifold.orientation_cover()
        def func(prec):
            return polished_tetrahedra_shapes(double_cover, bits_prec=prec)[::2]
    return ListOfApproximateAlgebraicNumbers(func)

@sage_method
def trace_field_gens(manifold, fundamental_group_args = []):
    """
    The generators of the trace field as ApproximateAlgebraicNumbers. Can be
    used to compute the tetrahedra field, where the first two parameters
    are bits of precision and maximum degree of the field::

        sage: M = Manifold('m125')
        sage: traces = M.trace_field_gens()
        sage: traces.find_field(100, 10, optimize=True)    # doctest: +NORMALIZE_WHITESPACE
        (Number Field in z with defining polynomial x^2 + 1
         with z = -1*I,
        <ApproxAN: -1.0*I>, [z + 1, z, z + 1])
    """
    def func(prec):
        return polished_holonomy(manifold, prec,
                                 fundamental_group_args).trace_field_generators()
    return ListOfApproximateAlgebraicNumbers(func)

@sage_method
def invariant_trace_field_gens(manifold, fundamental_group_args = []):
    """
    The generators of the trace field as ApproximateAlgebraicNumbers. Can be
    used to compute the tetrahedra field, where the first two parameters
    are bits of precision and maximum degree of the field::

        sage: M = Manifold('m007(3,1)')
        sage: K = M.invariant_trace_field_gens().find_field(100, 10, optimize=True)[0]
        sage: L = M.trace_field_gens().find_field(100, 10, optimize=True)[0]
        sage: K.polynomial(), L.polynomial()
        (x^2 - x + 1, x^4 - 2*x^3 + x^2 + 6*x + 3)
    """
    def func(prec):
        return polished_holonomy(manifold, prec,
                                 fundamental_group_args).invariant_trace_field_generators()
    return ListOfApproximateAlgebraicNumbers(func)

@sage_method
def holonomy_matrix_entries(manifold,
                            fundamental_group_args = [],
                            match_kernel = True):

    """
    The entries of the matrices of the holonomy as list of ApproximateAlgebraicNumbers
    (four consecutive numbers per matrix). The numbers are guaranteed to lie in the
    trace field only if match_kernel = False::

        sage: M = Manifold("m004")
        sage: mat_entries = M.holonomy_matrix_entries(match_kernel = False) # doctest: +NORMALIZE_WHITESPACE +NUMERIC9
        sage: mat_entries
        <SetOfAAN: [0.5 + 0.8660254037844386*I, 0.5 - 0.8660254037844386*I, 0.5 + 0.8660254037844386*I, 1.0 - 1.7320508075688772*I, 1.0 - 3.4641016151377544*I, -2.0 + 1.7320508075688772*I, -1.0 - 1.7320508075688772*I, 1.7320508075688772*I]>
        sage: K = mat_entries.find_field(100, 10, optimize = True)[0]
        sage: K.polynomial()
        x^2 - x + 1
    """

    def func(prec):
        G = polished_holonomy(manifold,
                              prec,
                              fundamental_group_args = fundamental_group_args,
                              match_kernel = match_kernel)
        return sum( [G.SL2C(g).list() for g in G.generators()], [])
    return ListOfApproximateAlgebraicNumbers(func)


def add_methods(mfld_class, hyperbolic=True):
    mfld_class.alexander_polynomial = nsagetools.alexander_polynomial
    mfld_class.homological_longitude = nsagetools.homological_longitude
    mfld_class.slice_obstruction_HKL = slice_obs_HKL.slice_obstruction_HKL
    if hyperbolic:
        mfld_class.polished_holonomy = polished_holonomy
        mfld_class.tetrahedra_field_gens = tetrahedra_field_gens
        mfld_class.trace_field_gens = trace_field_gens
        mfld_class.invariant_trace_field_gens = invariant_trace_field_gens
        mfld_class.holonomy_matrix_entries = holonomy_matrix_entries
        mfld_class.hyperbolic_torsion = nsagetools.hyperbolic_torsion
        mfld_class.hyperbolic_adjoint_torsion = nsagetools.hyperbolic_adjoint_torsion
        mfld_class.hyperbolic_SLN_torsion = nsagetools.hyperbolic_SLN_torsion

@sage_method
def trace_field(manifold, prec=None, degree=None):
    return ManifoldNT.ManifoldNT(spec=None, snappy_mfld=manifold).trace_field(prec, degree)

@sage_method
def invariant_trace_field(manifold, prec=None, degree=None):
    return ManifoldNT.ManifoldNT(spec=None, snappy_mfld=manifold).invariant_trace_field(prec, degree)

@sage_method
def quaternion_algebra(manifold, prec=None, degree=None):
    mfld_nt = ManifoldNT.ManifoldNT(spec=None, snappy_mfld=manifold)
    mfld_nt.trace_field(prec, degree)
    return mfld_nt.quaternion_algebra(prec)

@sage_method
def invariant_quaternion_algebra(manifold, prec=None, degree=None):
    mfld_nt = ManifoldNT.ManifoldNT(spec=None, snappy_mfld=manifold)
    mfld_nt.invariant_trace_field(prec, degree)
    return mfld_nt.invariant_quaternion_algebra(prec)

@sage_method
def denominators(manifold, prec=None, degree=None):
    mfld_nt = ManifoldNT.ManifoldNT(spec=None, snappy_mfld=manifold)
    mfld_nt.trace_field(prec, degree)
    return mfld_nt.denominators()

@sage_method
def denominator_residue_characteristics(manifold, prec=None, degree=None):
    mfld_nt = ManifoldNT.ManifoldNT(spec=None, snappy_mfld=manifold)
    mfld_nt.trace_field(prec, degree)
    return mfld_nt.denominator_residue_characteristics()

@sage_method
def is_arithmetic(manifold, prec=None, degree=None):
    mfld_nt = ManifoldNT.ManifoldNT(spec=None, snappy_mfld=manifold)
    mfld_nt.compute_arithmetic_invariants(prec, degree)
    return mfld_nt.is_arithmetic()

@sage_method
def p_arith(manifold, prec=None, degree=None):
    mfld_nt = ManifoldNT.ManifoldNT(spec=None, snappy_mfld=manifold)
    mfld_nt.compute_arithmetic_invariants(prec, degree)
    mfld_nt.p_arith()

@sage_method
def compare_arithmetic_invariants(manifold1, manifold2, prec=None, degree=None):
    mfld1_nt = ManifoldNT.ManifoldNT(spec=None, snappy_mfld=manifold1)
    mfld2_nt = ManifoldNT.ManifoldNT(spec=None, snappy_mfld=manifold2)
    mfld1_nt.compute_arithmetic_invariants(prec, degree)
    mfld2_nt.compute_arithmetic_invariants(prec, degree)
    return mfld1_nt.compare_arithmetic_invariants(mfld2_nt)

@sage_method
def has_same_arithmetic_invariants(manifold1, manifold2, prec=None, degree=None):
    mfld1_nt = ManifoldNT.ManifoldNT(spec=None, snappy_mfld=manifold1)
    mfld2_nt = ManifoldNT.ManifoldNT(spec=None, snappy_mfld=manifold2)
    mfld1_nt.compute_arithmetic_invariants(prec, degree)
    mfld2_nt.compute_arithmetic_invariants(prec, degree)
    return mfld1_nt.has_same_arithmetic_invariants(mfld2_nt)


def add_arith_methods(mfld_class):
    mfld_class.trace_field = trace_field
    mfld_class.invariant_trace_field = invariant_trace_field
    mfld_class.quaternion_algebra = quaternion_algebra
    mfld_class.invariant_quaternion_algebra = invariant_quaternion_algebra
    mfld_class.denominators = denominators
    mfld_class.denominator_residue_characteristics = denominator_residue_characteristics
    mfld_class.is_modtwo_homology_sphere = ManifoldNT.ManifoldNT.is_modtwo_homology_sphere
    mfld_class.is_arithmetic = is_arithmetic
    mfld_class.p_arith = p_arith
    mfld_class.compare_arithmetic_invariants = compare_arithmetic_invariants
    mfld_class.has_same_arithmetic_invariants = has_same_arithmetic_invariants