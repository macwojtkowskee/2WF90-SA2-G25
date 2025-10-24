# STILL IS MISSING SOME BEHAVIOUR: 
# I did the basic operations, but we still need some of the longer algorithms.

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
        self.mod = mod # It's just an integer :-)
        
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
            new_coeff = (self.get_coeff(i) + other_poly.get_coeff(i)) % self.modulus
            new_coeffs.append(new_coeff)
        
        # As I mentioned previously, the modulo operations are applied at the conclusion
        # of all basic operations to maintain conformity to p.
        return Polynomial(new_coeffs, self.modulus)
    
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
            new_coeff = (self.get_coeff(i) - other_poly.get_coeff(i)) % self.modulus
            new_coeffs.append(new_coeff)
            
        return Polynomial(new_coeffs, self.modulus)