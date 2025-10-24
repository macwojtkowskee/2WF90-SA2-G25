import poly_arithmetic as pa

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
    
    # Check if f^order = 1 using fast exponentiation with reduction
    result = power_mod(f, order, h)
    
    if result.coefficients != [1]:
        return False
    
    # Check that f^(order/q) != 1 for all prime divisors q of order
    factors = prime_factors(order)
    
    for q in factors:
        exp = order // q
        result = power_mod(f, exp, h)
        
        if result.coefficients == [1]:
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
    
    return result