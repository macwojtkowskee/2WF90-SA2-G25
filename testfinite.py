import unittest
import finite_field_arithmetic as ffa
import poly_arithmetic as pa

class TestFiniteArithmetic(unittest.TestCase):
    def P(self, coeffs, mod=5):
        """Helper to create a Polynomial with a default modulus of 5."""
        return pa.Polynomial(coeffs, mod)

    def A(self, coeffs):
        """Helper to create an expected list of coefficients."""
        # Use Polynomial to normalize coefficients, matching your __init__
        return pa.Polynomial(coeffs, 5).coefficients 
    
    def h(self):
        return pa.Polynomial([1, 1, 0, 1], 5)

    def test_ff_multiplication(self):
        # Base multiplication
        f = pa.Polynomial([3, 2, 0, 1], 5)  
        g = pa.Polynomial([0, 1, 1], 5) 
        t = ffa.finite_field_multiply(f, g, self.h())
        self.assertEqual(t.coefficients, [4, 1, 3])

        # Check commutativity
        t2 = ffa.finite_field_multiply(g, f, self.h())
        self.assertEqual(t2.coefficients, t.coefficients)

        # Check identity
        f = pa.Polynomial([0, 1], 5)
        g = pa.Polynomial([1], 5)
        t = ffa.finite_field_multiply(f, g, self.h())
        self.assertEqual(t.coefficients, [0,1])

        # Check multiplication by 0
        t = ffa.finite_field_multiply(f, pa.Polynomial([0], 5), self.h())
        self.assertEqual(t.coefficients, [0])

        # Check reduction 
        f = pa.Polynomial([1, 1, 1, 1, 1], 5)
        g = pa.Polynomial([1, 1, 1, 1, 1], 5)
        t = ffa.finite_field_multiply(f, g, self.h())
        self.assertEqual(t.coefficients, [0, 0, 1])
    
    def test_ff_inversion(self):
        f = pa.Polynomial([0, 1], 5)
        inv_f = ffa.finite_field_inversion(f, self.h())

        # Test if inversion of X is not none
        self.assertIsNotNone(inv_f)

        # Test if f * inv(f) == 1
        prod = ffa.finite_field_multiply(f, inv_f, self.h())
        self.assertEqual(prod.coefficients, [1])

        # Zero has no inverse
        zero_poly = pa.Polynomial([0], 5)
        inv_zero = ffa.finite_field_inversion(zero_poly, self.h())
        self.assertIsNone(inv_zero)

        # Should be None if f.mod != h.mod
        f = pa.Polynomial([0, 1], 3)
        inv_f = ffa.finite_field_inversion(f, self.h())
        self.assertIsNone(inv_f)

    def test_ff_division(self):
        # Basic division
        f = pa.Polynomial([3, 2, 0, 1], 5)   
        g = pa.Polynomial([0, 1, 1], 5)
        t = ffa.finite_field_division(f, g, self.h()) 
        self.assertEqual(t.coefficients, [1, 1, 2])

        # Division by 0 should not be possible
        f = pa.Polynomial([1, 1], 5)
        g = pa.Polynomial([0], 5)
        div = ffa.finite_field_division(f, g, self.h())
        self.assertIsNone(div)

        # Division when g == h should not be possible
        g = self.h()
        div = ffa.finite_field_division(f, g, self.h())
        self.assertIsNone(div)

    def test_ff_primitivity(self):
        # Base case for testing
        f = pa.Polynomial([0, 0, 2], 5)
        self.assertTrue(ffa.is_primitive(f, self.h(), f.mod))
        
        # Check if zero polynomial is primitive
        zero_poly = pa.Polynomial([0], 5)
        self.assertFalse(ffa.is_primitive(zero_poly, self.h(), 5))

        # Constants should not be primitive 
        const_poly = pa.Polynomial([7], 5)
        self.assertFalse(ffa.is_primitive(const_poly, self.h(), 5))

    def test_ffa_primitive_generation(self):
        # Base case 
        f = ffa.primitive_generation(self.h(), 5)
        self.assertTrue(ffa.is_primitive(f, self.h(), 5))
        
        # Extensive check
        g = ffa.primitive_generation(self.h(), 5)
        self.assertLess(g.degree(), self.h().degree())
        self.assertTrue(all(0 <= c < 5 for c in g.coefficients))
        self.assertEqual(g.mod, 5)
        
if __name__ == '__main__':
    # Add a note explaining how to run the tests
    print("--- Starting Polynomial Arithmetic Tests ---")
    print("To run tests from the command line, use: python -m unittest test_arithmetic.py")
    print("------------------------------------------")
    unittest.main(exit=False)
