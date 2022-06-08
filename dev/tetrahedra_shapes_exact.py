
"""
This function adds to features to find_field. 

First, it will later be 
"""
from sage.all import ComplexIntervalField

from snappy.snap.find_field import *



# this is a helper function that tells us which equations
# correspond to edges, cusp equations, and filling equations
# it returns a vector e, u, f
# where e is the number of edges
# u is the number of unfilled equations
# f is the number of filled equations
def types_of_gluing_equations(M):
    """
    >>> M = Manifold('m125')
    >>> M.dehn_fill((4,5),0)
    >>> types_of_gluing_equations(M)
    [4,2,1]

    """
    u = 0
    f = 0
    for ci in M.cusp_info():
        if ci.is_complete:
            u = u+2
        else:
            f = f+1
    return [M.num_tetrahedra(), u, f]

def fast_pi_approx(prec):
    if prec < 20:
        npi = RealInterval(3.14159265358979323846, 3.14159265358979323847)
    else:
        npi = RealInterval(pi.n(prec)-.1^prec, pi.n(10)+.1^prec)
    return npi

def full_shapes_for_embedding(shapes, zGen):
    # zGen should be an interval or a exact value
    interval_shapes = [shape.polynomial()(zGen) for shape in shapes]
    interval_z_primes = [1/(1-z) for z in interval_shapes]
    interval_z_douple_primes = [(z-1)/z for z in interval_shapes]

    return zip(interval_shapes, interval_z_primes, interval_z_double_primes)




def expected_angle_sums(M):
    """
    Given a triangulated manifold with n tetrahedra, u unfilled cusps, and f filled cusps
    returns a list of length n+u+f which gives the integral multiple of pi expected 
    for each of the gluing equations for M

    >>> M = Manifold('m125') 
    >>> M.dehn_fill((1,10),1)
    >>> expected_angle_sums(M)
    [2,2,2,2,0,0,2] 
    """
    angle_sums = [2 for i in range(M.num_tetrahedra())]

    for ci in M.cusp_info():
        if ci.is_complete:
            angle_sums.append(0)
            angle_sums.append(0)
        else:
            # todo: handle orbifolds/cone manifold here
            # 
            angle_sums.append(2)
    return angle_sums


def isolated_root_interval_for_embedding(K, root, prec): 
    """
    Given an embedding of a number field, this function returns an 
    complex interval with exactly one root in it 
    
    """
    # add more exposition about why this is enough
    p = K.defining_polynomial()
    CIF = ComplexIntervalField(prec)
    z = CIF(root.real(), root.imag())
    # add a check that z is actually contains a root
    if 0<abs(p.derivative()(z)):
        return z
    else:
        return None # add more here
    
    
# to do:
# add more to caching
def tetrahedra_shapes_exact(M, prec, deg, verify=False, optimize=False, verbosity = False):

    """
    Computes the tetrahedral shapes as algebraic numbers using find_field and the LLL algorithm. 
    Then checks that the gluing equations are satisfied algebraically.
    The options "optimize" and "verbosity"  are described in more detail `here <snap.html>`_::
      
      sage: M = Manifold('m016')
      sage: tetrahedra_shapes_exact(M, 100, 30)
      sage: K, zApprox, shapes = tetrahedra_shapes_exact(M, 100, 20)
      sage: K.degree()
      3
      sage: K.discriminant()
      -23
    
    If verified=True, the tetrahedral shapes are computed as algebraic numbers
    using find_field and the LLL algorithm.       
    Then the gluing equations are checked algebraically the additional check that the argument 
    angle sum conditions are satisfied. Here an interval which isolates the root of the polynomial
    defined in the embedding of the tetrahedral field is returned instead of the approximate roots. 


      sage: K, zVerifiedInterval, shapes = tetrahedra_shapes_exact(M, 100, 20, verify=True)
      sage: K.degree()
      3
      sage: K.discriminant()
      -23
       
    """
    tet_gens = M.tetrahedra_field_gens()
    K = tet_gens.find_field(prec, deg, optimize, verbosity)
    if K == None:
        return None


    # this function always checks that the gluing equations are satisfied
    
    shape_field, approx_gen, shapes = K

    
    # todo: be more careful about selecting equations
    #       only select on peripheral equation of an unfilled cusp
    #       and only the necessary n-c equations

    # this for loop goes through each gluing equation
    # checking that the equation evaluated at
    # the solutions for the tetrahedral shapes is 1
    # in order to avoid inverting 1-z in z'= 1/(1-z) and
    # z in z''=(z-1)/z
    # it bulids a left-hand side and a right-hand side of an
    # equation.
    # For each positive power of z that arises the left-hand side
    # multiplied by z and the right-hand side is left alone 
    # for each positive power of z' that comes up the left-hand side is left alone
    # and the right-hand side gets an extra factor of 1-z
    # finally for each z'' the left-hand side gets a factor of
    # z-1 and the right-hand side gets a factor of z
    for exponents in M.gluing_equations():
        left_side = 1
        right_side = 1
        for i, exponent in enumerate(exponents):

            if exponent == 0:
                continue
            
            tet_num = i//3
            cur_shape = shapes[tet_num]

            if i % 3 == 0:
                 if exponent > 0:
                    left_side *= cur_shape**exponent)
                elif exponent < 0:
                    right_side *= cur_shape**-exponent
                    
            elif i % 3 == 1: # deal with z'=1/(1-z)
                if exponent > 0:
                    right_side *= (1-cur_shape)**exponent
                else:
                    left_side *= (1-cur_shape)**-exponent
                    
            elif exponent != 0: # deal with z''= (z-1)/z
                if exponent > 0:
                    left_side *= (cur_shape-1)**exponent
                    right_side *= (shapes[int((i-2)/3)])**exponents[i])
                else:
                    right_side *= (cur_shape-1)**-exponent
                    left_side *= (cur_shape)**-exponent


        if left_side - right_side != 0:

            raise ValueError("Equation: ", exponents, "and shapes", shapes, " is not 1. Left side:", left_side, " right side:", right_side, " difference:", left_side - right_side, " aprrox difference:", (left_side - right_side).n())
            return None

    if not verify:
        return K # consider also adding shapes to cache at this stage
    else:
        t, u, f = types_of_gluing_equations(M)
        if f ==0: # all cusps are unfilled
            # we just check that all tet shapes have positive imaginary part
            # Theorem: Given an algebraic solution to the gluing equations
            #          where all tetrahedral shapes have positive imaginary part
            #          the logarithmic gluing equations are satisfied.
            #
            # Proof (check Benedetti and Petronio or recall nathan's argument

            zInt = isolated_root_interval_for_embedding(shape_field, shape_field.gen().n(), prec)

            for shape in shapes:
                if 0 >= (shape.polynomial()(zInt)).imag().lower():
                    raiseValueError("Could not verify shapes are positive")
                    #return None
            return [K[0], zInt, K[2]]
        else:
            print(" to do")

            interval_shapes = [shape.polynomial()(zInt) for shape in shapes]
            interval_z_primes = [1/(1-z) for z in interval_shapes]
            interval_z_douple_primes = [(z-1)/z for z in interval_shapes]
            
            full_set_of_interval_shapes =zip(interval_shapes, interval_z_primes, interval_z_double_primes)

            for int_shape in interval_shapes:
                if 0 >= (int_shape.imag()).lower():
                    raiseValueError("Could not verify shapes")
                    return None
            expected_arguments = expected_angle_sums(M)

            pi_approx = fast_pi_approx(10) # this should be sufficent

            for i in range(len(M.gluing_equations())):
                cur_equation = M.gluing_equations()[i]
                cur_argument = 0
                for j in range(len(cur_equation)):
                    if cur_equation !=0:
                        cur_argument = cur_equation[j]*(full_set_of_interval_shapes[j]).arg()
                if cur_argument.lower() < ((expected_arguments-2)*pi_approx).upper() or cur_argument.upper() < ((expected_arguments+2)*pi_approx).lower():
                    raiseValueError("Could not verify arguments for ", M.gluing_equations()[i], " with shapes ", full_set_of_interval_shapes[j])
                    return None
            
                   
                                  
            return [K[0], zInt, K[2]]
           
            # check shapes and arguments
