import poly_arithmetic as pa
import random 

def poly_mod_reduction(f, h):
    """
        Helper function to reduce polynomial f modulo
        Polynomial h. Uses LD from poly arithmetic.
        Is used for finite field addition, multiplication and subtraction
    """
    r = pa.polynomial_LD(f, h)[1]
    if r is None:
        return pa.Polynomial([0], f.mod)
    return r

def finite_field_inversion(f: pa.Polynomial, h: pa.Polynomial):
    """
    Obtain the multiplicative inverse of a given polynomial f in the field Z/pZ/(h)
    This uses the extended Euclidean algorithm to find polynomials a, b such that:
        a*f + b*h = d = gcd(f, h).
    
    Args:
        f: polynomial to invert
        h: modulus polynomial in the given field; irreducible

    """

    # Sanity check for the mods of polynomials
    if f.mod != h.mod:
        return None
    
    p = f.mod

    f = poly_mod_reduction(f,h)

    # zero has no inverse
    if all((c % p) == 0 for c in f.coefficients):
        return None

    # Compute gcd(f,h) and a,b s.t a*f + b*h = d = gcd(f, h)
    a, b, d = pa.poly_extended_euclidean_algorithm(f, h)

    # invertible iff d == 1 (mod p)
    if not (len(d.coefficients) == 1 and (d.coefficients[0] % p) == 1):
        return None

    inv = poly_mod_reduction(a, h)
    return inv

def prime_factors(n):
    """
        Get unique prime factors of n.
        Trial and error approach
        Used for primality check
    """
    
    # NOTE: design decision to use set data structure, 
    # as sets provide uniqueness
    factors = set()
    
    # We start with the smallest prime
    d = 2
    # Only check divisors up to sqrt(n), because
    # NOTE: after sqrt(n) we get the same pairs of divisors, 
    # but in a swapped order, e.g not 3*21 but 21*3, so
    # it is quicker to check just up until sqrt(n)
    while d * d <= n:
        while n % d == 0:
            factors.add(d)
            n //= d
        d += 1
    if n > 1:
        factors.add(n)
    return factors


def finite_field_multiply(f, g, h):
    """
        Multiply f and g in the finite field Z_p[X]/(h).
        Reduces degree during multiplication to keep coefficients in check.
    """
    p = f.mod
    deg_h = h.degree()
    
    if f.degree() == -1 or g.degree() == -1:
        return pa.Polynomial([0], p)
    
    # Start with zero polynomial
    result = pa.Polynomial([0], p)
    
    # Multiply term by term and reduce modulo h frequently
    for i in range(len(g.coefficients)):
        # Multiply f by g[i] * X^i
        term_coeffs = [0] * i + [c * g.coefficients[i] % p for c in f.coefficients]
        term = pa.Polynomial(term_coeffs, p)
        
        # Add to result
        result = result + term
        
        # Reduce modulo h if degree gets too large
        # This prevents the result from having huge degree
        if result.degree() >= 2 * deg_h:
            result = poly_mod_reduction(result, h)
    
    # Final reduction
    result = poly_mod_reduction(result, h)
    return result

def finite_field_division(f: pa.Polynomial, g: pa.Polynomial, h: pa.Polynomial):
    """
    Divide f by g in the finite field Z/pZ/(h), by computing the product f * g^(-1)

    Args:
        f: Dividient polynomial in the given field
        g: Divisor polynomial in the given field
        h: modulus polynomial in the given field; irreducible
    """
    p = f.mod

    # Reduce the divisor by the appropriate modulus
    g = poly_mod_reduction(g, h)
    
    # Division by zero is not allowed
    if len(g.coefficients) == 1 and g.coefficients[0] % p == 0:
        return None

    g_inv = finite_field_inversion(g, h)

    # Divide by multipling dividient with the inverse of the divisor
    prod = finite_field_multiply(f, g_inv,h)
    prod = poly_mod_reduction(prod, h)

    return prod

def is_primitive(f, h, p):
    """
        Check if f is a primitive element in Z_p[X]/(h).
        f is primitive if its order is p^deg(h) - 1.
        Uses modular exponentiation with reduction mod h.
    """
    n = h.degree()
    order = p ** n - 1
    
    if f.degree() == -1:
        return False
    
    if n > 1 and f.degree() == 0:
        return False
    
    # Check if f^order = 1 using fast exponentiation with reduction
    result = power_mod(f, order, h)
    
    if not (len(result.coefficients) == 1 and result.coefficients[0] == 1):
        return False
    
    # Check that f^(order/q) != 1 for all prime divisors q of order
    factors = prime_factors(order)
    
    for q in factors:
        exp = order // q
        result = power_mod(f, exp, h)
        
        # If result equals 1, then f is NOT primitive
        if len(result.coefficients) == 1 and result.coefficients[0] == 1:
            return False
    
    return True


def power_mod(base, exp, h):
    """
        Helper function to compute base^exp mod h efficiently using RTL square and multiply method.
        Reduces modulo h at each step to keep polynomials small.
    """
    p = base.mod
    result = pa.Polynomial([1], p)
    base = pa.Polynomial(base.coefficients[:], p)
    
    # Ensure base is already reduced
    base = poly_mod_reduction(base, h)
    
    while exp > 0:
        # check the least significant bit
        if exp % 2 == 1:
            result = finite_field_multiply(result, base, h)
        # square the base, reduce mod h
        base = finite_field_multiply(base, base, h)
        # shift the "cursor"
        exp //= 2
    
    result = poly_mod_reduction(result, h)
    return result

def generate_polynomial_h(h: pa.Polynomial):
    """
    Helper function to generate a random polynomial in Z/pZ/(h) with deg < deg(h)
    Coefficients are sampled uniformly from the set {0, 1, ..., p-1}

    Args:
        h: modulus polynomial in the given field; irreducible; determines the maximal degree
    """
    p = h.mod
    n = h.degree()
    
    # Generate polynomial with degree < deg(h)
    coeffs = [random.randint(0, p - 1) for _ in range(n)]
    return pa.Polynomial(coeffs, p)

def primitive_generation(h: pa.Polynomial, p: int):
    """
    Compute a primitive polynomial in a given field Z/pZ/(h) by randomly sampling polynomials 
        until a primitive one is found.

    Args:
        h: modulus polynomial in the given field; irreducible
        p: Prime modulus of the coefficient field
    """
    p = h.mod
    
    while True:
        # If a primitive element is found, return the polynomial
        f = generate_polynomial_h(h)
        if is_primitive(f, h, p):
            return f