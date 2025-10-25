# Finally finished
import random

class Polynomial:
    """ 
    
    Represents a polynomial in Z[X]/pZ, for a prime p modulus.
    To avoid the inversion/reversion issues that we faced in the
    first assignment, all of the coefficients are stored in ascending order
    (since python iterates 0 -> n):
    
    [C_0 + C_1 * X + C_2 * X^2]
    
    All the modulo prime p operations are applied stepwise, to maintain
    the form-appropriate form.
    
    """
    def __init__ (self, coefficients, mod):
        """
        
        Simple initialization, same as SA1, but here we provide the logic
        for instantenous formalization into canonical form (so that any
        [logical] input can work).

        Args:
            coefficients (list[int]): List of coefficients, given in order
                                      described in class docstring.
            mod (int): Prime p (for some Z/pZ) 
        """
        self.mod = mod # It's just an integer this time :-)
        
        # Simply use the built-in modulo function to make sure moduli appropriate,
        # then trim leading zeroes (also for input robustness, didn't have time last time.)
        temp_coeffs = [c % mod for c in coefficients]
        
        while len(temp_coeffs) > 1 and temp_coeffs[-1] == 0:
            temp_coeffs.pop()
        if len(temp_coeffs) == 1 and temp_coeffs[0] == 0:
            self.coefficients = [0]
        else:
            self.coefficients = temp_coeffs
            
    def degree(self):
        """
            Returns the degree of the polynomial. 
            As provided in the assignment description, degree(-1) indicates
            the zero polynomial.
        """
        if self.coefficients == [0]:
            return -1
        return len(self.coefficients) - 1
    
    def get_coefficient(self, index):
        """
            Gets the coefficient of a polynomial variable, based on its index
            (its degree.) Required for any sort of basic operation to properly
            increment/decrement (easier than doing it per-function)
        """
        if index < len(self.coefficients):
            return self.coefficients[index]
        return 0
    
    # This time, to avoid the overtly long  naming conventions in the first submission,
    # we use dunders (double leading/trailing underscores) - we only used this for the
    # __init__ function in SA1, but these functions are also crucial to the behaviour
    # of the objects, and so to reduce overhead it makes sense to utilize them for all 
    # the methods below.
    
    def __add__(self, other_poly):
        """
        
            Performs polynomial addition: (polynomial + second polynomial) mod p.
            Args:
                self: base polynomial for addition
                other_poly: second polynomial to be used
        
        """
        # Note that 
        max_len = max(len(self.coefficients), len(other_poly.coefficients))
        new_coeffs = []
        
        # Incredibly simple incrementation of coefficients, since we're dealing with
        # integer operations. And we don't have to do any inversions...
        for i in range(max_len):
            new_coeff = (self.get_coefficient(i) + other_poly.get_coefficient(i)) % self.mod
            new_coeffs.append(new_coeff)
        
        # As I mentioned previously, the modulo operations are applied at the conclusion
        # of all basic operations to maintain conformity to p.
        return Polynomial(new_coeffs, self.mod)
    
    # Forgive the constant comparisons to the implementation of SA1 - a lot of content
    # here is similar or even recycled (for example, the maximal length setters), and
    # some feedback (mainly self-critique for now) has been incorporated into these 
    # methods.

    def __sub__(self, other_poly):
        """
            Performs polynomial subtraction: (polynomial - second polynomial) mod p.
            Analogous to '__add__'
            Args:
                self: base polynomial for addition
                other_poly: second polynomial to be used
        """
        max_len = max(len(self.coefficients), len(other_poly.coefficients))
        new_coeffs = []
        
        for i in range(max_len):
            new_coeff = (self.get_coefficient(i) - other_poly.get_coefficient(i)) % self.mod
            new_coeffs.append(new_coeff)
            
        return Polynomial(new_coeffs, self.mod)
    
    def __mul__(self, other):
        """
            Performs polynomial multiplication: (self * other) mod p.
            
        """
        deg_self = self.degree()
        deg_other = other.degree()
        
        if deg_self == -1 or deg_other == -1:
            return Polynomial([0], self.mod)

        new_deg = deg_self + deg_other
        new_coeffs = [0] * (new_deg + 1)
        
        # Whereas the logic above simply fetches the necessary data and increments
        # the appropriate lists to reflect the new polynomial, this handles the
        # incrementation of the coefficients.
        
        for i in range(len(self.coefficients)):
            for j in range(len(other.coefficients)):
                new_coeffs[i + j] = (new_coeffs[i + j] 
                                     + self.coefficients[i] 
                                     * other.coefficients[j]) % self.mod
                
        return Polynomial(new_coeffs, self.mod)
    
# The rest of the logic is handled as a regular function, though named in a manner
# less obnoxious than in the standard arithmetics implementation, and uses all of the
# (basic) functions given above.
# [NOTE]: This also removes the issue of having to call every single function on the 
# polynomial function we are working with itself, which was... a choice.

def polynomial_LD(f, g):
    """
    
    Performs long division f/g, using the script-provided algorithm.

    Args:
        f (Polynomial): The first polynomial.
        g (Polynomial): The second polynomial.
    """
    
    p = f.mod
    q = Polynomial([0], p)
    r = Polynomial(f.coefficients, p)
    
    # Calculate modular inverse of the leading coefficient of g
    g_lead_coeff = g.coefficients[-1]
    # Use pow(base, exponent, modulus) for modular inverse
    # exponent = -1 in Z_p is p-2 by Fermat's Little Theorem
    g_lead_coeff_inv = pow(g_lead_coeff, p - 2, p)
    
    while r.degree() >= g.degree() and r.degree() != -1:
        deg_r = r.degree()
        r_lead_coeff = r.coefficients[-1]
        
        # Calculate the coefficient and degree of the term to subtract
        term_coeff = (r_lead_coeff * g_lead_coeff_inv) % p
        term_degree = deg_r - g.degree()
        
        # Create the term as a polynomial: (term_coeff) * X^(term_degree)
        term_coeffs = [0] * (term_degree + 1)
        term_coeffs[term_degree] = term_coeff
        term = Polynomial(term_coeffs, p)
        
        # q = q + term
        # r = r - (term * g)
        q = q + term
        r = r - (term * g)
        
    return q, r

def poly_extended_euclidean_algorithm(f, g):
    """
    Performs the Extended Euclidean Algorithm for polynomials f and g.
    Returns (a, b, d) such that a*f + b*g = d = gcd(f, g).
    
    Adapts the outcome to be monic via comparison/multiplication by inverse.
    
    Args:
        f (Polynomial): The first polynomial.
        g (Polynomial): The second polynomial.
    """
    p = f.mod
    
    if g.degree() == -1:
        if f.degree() == -1: # Both are zero
            return Polynomial([0], p), Polynomial([0], p), Polynomial([0], p)
        
        # For now, 
        lead_coeff = f.coefficients[-1]
        lead_coeff_inv = pow(lead_coeff, p - 2, p)
        inv_poly = Polynomial([lead_coeff_inv], p)
        
        d = f * inv_poly
        a = inv_poly
        b = Polynomial([0], p)
        return a, b, d

    r_prev, r_curr = f, g
    a_prev, a_curr = Polynomial([1], p), Polynomial([0], p)
    b_prev, b_curr = Polynomial([0], p), Polynomial([1], p)

    while r_curr.degree() != -1:
        q, r_next = polynomial_LD(r_prev, r_curr)
        
        r_prev, r_curr = r_curr, r_next
        a_prev, a_curr = a_curr, a_prev - q * a_curr
        b_prev, b_curr = b_curr, b_prev - q * b_curr
        
    # r_prev is the gcd, but may not be monic
    d = r_prev
    a = a_prev
    b = b_prev
    
    # As the description notes - 'ensure the outcome is monic':
    # if the lead coefficient is already minimnal, we are OK;
    # otherwise, need to normalize by multiplying by inverse of
    # lead coefficient.
    
    lead_coeff = d.coefficients[-1]
    if lead_coeff == 0: # Should only be the zero polynomial
        return a, b, d 
        
    lead_coeff_inv = pow(lead_coeff, p - 2, p)
    inv_poly = Polynomial([lead_coeff_inv], p)
    
    # Normalize as described above.
    d_monic = d * inv_poly
    a_monic = a * inv_poly 
    b_monic = b * inv_poly
    
    return a_monic, b_monic, d_monic

def poly_irreducibility_check(f):
    """
    Checks if a polynomial f is irreducible in Z_p[X].
    Uses the property that f (deg n) is reducible iff it has a
    divisor of degree k, where 1 <= k <= n/2.
    
    Much like other alrgorithms, this is based on the algorithm 
    found in Chapter 7 of the script (p. 71-72.)
    """
    n = f.degree()
    p = f.mod
    
    if n == 1:
        return True
        
    max_divisor_deg = n // 2
    
    for deg_k in range(1, max_divisor_deg + 1):
        
        # Iterate through all monic polynomials of degree deg_k
        # There are p^deg_k such polynomials, and so
        # we iterate through all of the p^deg_k combinations for the
        # coefficients.
        
        num_polys = p ** deg_k
        for i in range(num_polys):
            coeffs = [0] * (deg_k + 1)
            # simply SET the degree to be monic
            coeffs[deg_k] = 1 
            
            temp_i = i
            for j in range(deg_k):
                coeffs[j] = temp_i % p
                temp_i //= p
            
            # Simple check for divisibility (just ver deg)
            divisor = Polynomial(coeffs, p)
            _,r = polynomial_LD(f, divisor)
            if r.degree() == -1:
                # f is divisible by the potential divisor, so it's reducible!
                return False
                
    # The polynomial is irreducible if no divisors are found up to degree n/2, so:
    return True

def poly_generate_irreducible(p, n):
    """
    Generates a random irreducible polynomial of degree n in Z/pZ.
    This is as simple as sampling a random polynomial in the given range,
    and using the previous
    """
    while True:
        # Generate a random polynomial of the given degree n:
        coeffs = [0] * (n + 1)
        for i in range(n):
            coeffs[i] = random.randint(0, p - 1)
        coeffs[n] = 1 # Again, make it monic
        
        f = Polynomial(coeffs, p)
        
        if poly_irreducibility_check(f):
            return f
