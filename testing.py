import unittest
from poly_arithmetic import (
    Polynomial, 
    polynomial_LD, 
    poly_extended_euclidean_algorithm, 
    poly_irreducibility_check, 
    poly_generate_irreducible
)

class TestPolynomialArithmetic(unittest.TestCase):

    # --- Utility Functions for easy Polynomial creation ---

    def P(self, coeffs, mod=5):
        """Helper to create a Polynomial with a default modulus of 5."""
        return Polynomial(coeffs, mod)

    def A(self, coeffs):
        """Helper to create an expected list of coefficients."""
        # Use Polynomial to normalize coefficients, matching your __init__
        return Polynomial(coeffs, 5).coefficients 

    # --- Test Cases for Polynomial Class Initialization and Canonical Form ---
    
    def test_initialization_and_canonical_form(self):
        # Case 1: Standard polynomial, no modulus reduction needed
        p1 = self.P([1, 2, 3]) # 1 + 2X + 3X^2 mod 5
        self.assertEqual(p1.coefficients, [1, 2, 3])
        self.assertEqual(p1.degree(), 2)
        self.assertEqual(p1.mod, 5) # Check correct attribute name

        # Case 2: Zero polynomial
        p2 = self.P([0, 0, 0])
        self.assertEqual(p2.coefficients, [0])
        self.assertEqual(p2.degree(), -1)

        # Case 3: Leading zeros trimmed
        p3 = self.P([1, 2, 0, 0])
        self.assertEqual(p3.coefficients, [1, 2])
        self.assertEqual(p3.degree(), 1)
        
        # Case 4: Coefficients reduced modulo p
        p4 = self.P([6, 11, -3]) # mod 5: 1, 1, 2
        self.assertEqual(p4.coefficients, [1, 1, 2])
        self.assertEqual(p4.degree(), 2)
        
        # Case 5: Zero polynomial that needs reduction
        p5 = self.P([5, 10], mod=5)
        self.assertEqual(p5.coefficients, [0])
        self.assertEqual(p5.degree(), -1)

    # --- Test Cases for Basic Arithmetic Operations (+, -, *) ---

    def test_addition(self):
        # (1 + 3X^2) + (2 + 4X) mod 5 = 3 + 4X + 3X^2
        f = self.P([1, 0, 3])
        g = self.P([2, 4])
        h = f + g
        self.assertEqual(h.coefficients, self.A([3, 4, 3]))

        # Addition resulting in zero
        f = self.P([1, 2])
        g = self.P([4, 3]) # -1 -2X mod 5
        h = f + g
        self.assertEqual(h.coefficients, self.A([0]))
        
    def test_subtraction(self):
        # (1 + 3X^2) - (2 + 4X) mod 5 = (1-2) + (0-4)X + 3X^2 = 4 + X + 3X^2
        f = self.P([1, 0, 3])
        g = self.P([2, 4])
        h = f - g
        self.assertEqual(h.coefficients, self.A([4, 1, 3]))

        # Subtraction resulting in zero
        f = self.P([1, 2, 3])
        g = self.P([1, 2, 3])
        h = f - g
        self.assertEqual(h.coefficients, self.A([0]))
        
        # Subtraction resulting in coefficient reduction
        f = self.P([4, 0]) # 4
        g = self.P([1])   # 1
        h = f - g # 3
        self.assertEqual(h.coefficients, self.A([3]))

    def test_multiplication(self):
        # (2 + X) * (1 + 2X) mod 5
        # = 2*1 + 2*2X + X*1 + X*2X = 2 + 5X + 2X^2 = 2 + 0X + 2X^2 (mod 5)
        f = self.P([2, 1])
        g = self.P([1, 2])
        h = f * g
        self.assertEqual(h.coefficients, self.A([2, 0, 2])) 
        
        # Multiplication by zero
        f = self.P([1, 2, 3])
        g = self.P([0])
        h = f * g
        self.assertEqual(h.coefficients, self.A([0]))

        # Higher degree multiplication
        # (1 + X + 2X^2) * (1 + 3X^2) mod 5
        # = 1(1+3X^2) + X(1+3X^2) + 2X^2(1+3X^2)
        # = 1 + 3X^2 + X + 3X^3 + 2X^2 + 6X^4
        # = 1 + X + 5X^2 + 3X^3 + X^4
        # = 1 + X + 0X^2 + 3X^3 + X^4 (mod 5)
        f = self.P([1, 1, 2])
        g = self.P([1, 0, 3])
        h = f * g
        self.assertEqual(h.coefficients, self.A([1, 1, 0, 3, 1]))

    # --- Test Cases for Long Division ---

    def test_polynomial_LD(self):
        # Example 1: (X^3 + X + 1) / (X^2 + 1) in Z_2[X]
        # X^3 + X + 1 = X(X^2 + 1) + 1. Quotient = X, Remainder = 1
        p = 2
        f = Polynomial([1, 1, 0, 1], p) # X^3 + X + 1
        g = Polynomial([1, 0, 1], p)   # X^2 + 1
        q_expected = Polynomial([0, 1], p) # X
        r_expected = Polynomial([1], p)    # 1

        q, r = polynomial_LD(f, g)
        self.assertEqual(q.coefficients, q_expected.coefficients)
        self.assertEqual(r.coefficients, r_expected.coefficients)

        # Example 2: (X^4 + 3X^3 + 2X^2 + X + 1) / (X^2 + 4X + 1) in Z_5[X]
        # (X^4 + 3X^3 + 2X^2 + X + 1) / (X^2 + 4X + 1) = (X^2 + 4X + 1) + (4X)
        p = 5
        f = Polynomial([1, 1, 2, 3, 1], p) # X^4 + 3X^3 + 2X^2 + X + 1
        g = Polynomial([1, 4, 1], p)       # X^2 + 4X + 1


        # q = X^2 + 4X
        # r = 2X + 1
        q_expected = Polynomial([0, 4, 1], p)
        r_expected = Polynomial([1, 2], p)

        q, r = polynomial_LD(f, g)
        self.assertEqual(q.coefficients, q_expected.coefficients)
        self.assertEqual(r.coefficients, r_expected.coefficients)

        # Example 3: Division where remainder is zero
        f = Polynomial([2, 0, 3], p) # 2 + 3X^2
        g = Polynomial([4, 1], p)    # 4 + X
        # (2 + 3X^2) / (4 + X) = (3X + 3)
        # (3X + 3)(X + 4) = 3X^2 + 12X + 3X + 12 = 3X^2 + 0X + 2 (mod 5)
        q_expected = Polynomial([3, 3], p)
        r_expected = Polynomial([0], p)
        
        q, r = polynomial_LD(f, g)
        self.assertEqual(q.coefficients, q_expected.coefficients)
        self.assertEqual(r.coefficients, r_expected.coefficients)


    # --- Test Cases for Extended Euclidean Algorithm (EEA) ---

    def test_poly_eea_coprime(self):
        # Example 1: EEA with coprime polynomials in Z_5[X]
        # f = X^3 + 2X^2 + 3X + 4
        # g = X^2 + 3
        p = 5
        f = Polynomial([4, 3, 2, 1], p)
        g = Polynomial([3, 0, 1], p)
        
        # Expected gcd is 1 (monic)
        d_expected = Polynomial([1], p)
        
        a, b, d = poly_extended_euclidean_algorithm(f, g)

        self.assertEqual(d.coefficients, d_expected.coefficients, "GCD must be 1 (coprime)")
        
        # Check Bêzout's identity: a*f + b*g = d
        lhs = a * f + b * g

        _q_lhs, r_lhs = polynomial_LD(lhs, d) # Check if lhs - d = 0
        
        self.assertEqual(lhs.coefficients, d.coefficients, "Bézout's identity failed: a*f + b*g != d")
        
    def test_poly_eea_common_factor(self):
        # Example 2: EEA with a common factor in Z_3[X]
        # f = X^3 + 2X^2 + 2X + 1 = (X+1)(X^2+X+1) mod 3
        # g = X^2 + 2X + 1 = (X+1)^2 mod 3
        p = 3
        f = Polynomial([1, 2, 2, 1], p) 
        g = Polynomial([1, 2, 1], p)
        
        # Expected gcd: X+1 (monic)
        d_expected = Polynomial([1, 1], p)

        a, b, d = poly_extended_euclidean_algorithm(f, g)

        self.assertEqual(d.coefficients, d_expected.coefficients, "GCD must be X+1")
        
        # Check Bêzout's identity: a*f + b*g = d
        lhs = a * f + b * g

        _q_lhs, r_lhs = polynomial_LD(lhs, d) 
        self.assertEqual(r_lhs.degree(), -1, "Bézout's identity failed: Remainder should be zero after division by GCD")
        self.assertEqual(lhs.coefficients, d.coefficients, "Bézout's identity failed: a*f + b*g != d")

    def test_poly_eea_non_monic_input(self):
        # Test non-monic inputs to ensure monic output for gcd
        # f = 2X^2 + X + 1 mod 3
        # g = 2X + 1 mod 3
        p = 3
        f = Polynomial([1, 1, 2], p)
        g = Polynomial([1, 2], p)

        # (2X^2 + X + 1) = X * (2X + 1) + 1.
        # (2X + 1) = (2X+1) * 1 + 0.
        # The last non-zero remainder is 1. The GCD is 1.
        d_expected = Polynomial([1], p) # GCD is 1

        a, b, d = poly_extended_euclidean_algorithm(f, g)

        self.assertEqual(d.coefficients, d_expected.coefficients, "GCD must be 1")
        
        # Check Bêzout's identity: a*f + b*g = d
        lhs = a * f + b * g
        self.assertEqual(lhs.coefficients, d.coefficients, "Bézout's identity failed: a*f + b*g != d")


    # --- Test Cases for Irreducibility Check and Generation ---

    # Renamed test and function calls
    def test_poly_irreducibility_check_small_degree(self):
        # Test Z_2[X]
        p = 2
        # f = X^2 + X + 1 (Irreducible)
        self.assertTrue(poly_irreducibility_check(Polynomial([1, 1, 1], p)))
        
        # f = X^2 + 1 = (X+1)^2 (Reducible)
        self.assertFalse(poly_irreducibility_check(Polynomial([1, 0, 1], p)))

        # f = X^3 + X + 1 (Irreducible)
        self.assertTrue(poly_irreducibility_check(Polynomial([1, 1, 0, 1], p)))

        # f = X^4 + X^3 + 1 (Irreducible)
        self.assertTrue(poly_irreducibility_check(Polynomial([1, 0, 0, 1, 1], p)))
        
        # Test Z_3[X]
        p = 3
        # f = X^2 + 1 (Irreducible since X^2+1 has no roots in Z_3: 0^2+1=1, 1^2+1=2, 2^2+1=5=2)
        self.assertTrue(poly_irreducibility_check(Polynomial([1, 0, 1], p)))

        # f = X^2 + X + 1 (Reducible since root is 1: 1^2+1+1=3=0)
        self.assertFalse(poly_irreducibility_check(Polynomial([1, 1, 1], p)))
        
    def test_poly_generate_irreducible_small_degree(self):
        p = 2
        n = 3
        # Generate an irreducible polynomial of degree 3 in Z_2[X]
        f = poly_generate_irreducible(p, n)
        self.assertEqual(f.mod, p)
        self.assertEqual(f.degree(), n)
        self.assertTrue(poly_irreducibility_check(f))
        self.assertEqual(f.coefficients[-1], 1, "Generated polynomial must be monic")
        
        p = 5
        n = 2
        # Generate an irreducible polynomial of degree 2 in Z_5[X]
        f = poly_generate_irreducible(p, n)
        self.assertEqual(f.mod, p)
        self.assertEqual(f.degree(), n)
        self.assertTrue(poly_irreducibility_check(f))
        self.assertEqual(f.coefficients[-1], 1, "Generated polynomial must be monic")


if __name__ == '__main__':
    # Add a note explaining how to run the tests
    print("--- Starting Polynomial Arithmetic Tests ---")
    print("To run tests from the command line, use: python -m unittest test_arithmetic.py")
    print("------------------------------------------")
    unittest.main(exit=False)