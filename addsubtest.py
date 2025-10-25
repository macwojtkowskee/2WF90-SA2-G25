import unittest
import poly_arithmetic as pa
import finite_field_arithmetic as ffa  # we only need poly_mod_reduction here

class TestFFAddSub(unittest.TestCase):
    def h(self):
        # GF(5)[x]/(x^3 + x + 1)
        return pa.Polynomial([1, 1, 0, 1], 5)

    # ---------- ADDITION ----------
    def test_ff_add_no_reduction(self):
        f = pa.Polynomial([2, 3], 5)
        g = pa.Polynomial([4, 1], 5)
        s = f + g
        s = ffa.poly_mod_reduction(s, self.h())
        self.assertEqual(s.coefficients, [1, 4])

    def test_ff_add_with_reduction(self):
        f = pa.Polynomial([1, 0, 0, 1], 5) 
        g = pa.Polynomial([4], 5)            
        s = f + g
        s = ffa.poly_mod_reduction(s, self.h())
        self.assertEqual(s.coefficients, [4, 4])

    def test_ff_add_identity(self):
        f = pa.Polynomial([3, 2, 1], 5)
        z = pa.Polynomial([0], 5)
        s = f + z
        s = ffa.poly_mod_reduction(s, self.h())
        self.assertEqual(s.coefficients, [3, 2, 1])

    def test_ff_add_commutativity(self):
        f = pa.Polynomial([4, 0, 1], 5)
        g = pa.Polynomial([3, 3], 5)
        s1 = ffa.poly_mod_reduction(f + g, self.h())
        s2 = ffa.poly_mod_reduction(g + f, self.h())
        self.assertEqual(s1.coefficients, s2.coefficients)

    # ---------- SUBTRACTION ----------
    def test_ff_sub_no_reduction(self):
        f = pa.Polynomial([0, 1], 5)  
        g = pa.Polynomial([3], 5)       
        d = f - g
        d = ffa.poly_mod_reduction(d, self.h())
        self.assertEqual(d.coefficients, [2, 1])

    def test_ff_sub_with_reduction(self):
        f = pa.Polynomial([0, 0, 0, 1], 5) 
        g = pa.Polynomial([0, 0, 1], 5)   
        d = f - g
        d = ffa.poly_mod_reduction(d, self.h())
        self.assertEqual(d.coefficients, [4, 4, 4])

    def test_ff_sub_identity(self):
        f = pa.Polynomial([1, 1, 1], 5)
        z = pa.Polynomial([0], 5)
        d = f - z
        d = ffa.poly_mod_reduction(d, self.h())
        self.assertEqual(d.coefficients, [1, 1, 1])

    def test_ff_sub_self_is_zero(self):
        f = pa.Polynomial([2, 4, 3], 5)
        d = f - f
        d = ffa.poly_mod_reduction(d, self.h())
        self.assertEqual(d.coefficients, [0])

    def test_ff_sub_anti_commutativity(self):
        f = pa.Polynomial([1, 2, 0, 1], 5)
        g = pa.Polynomial([4, 3], 5)
        d1 = ffa.poly_mod_reduction(f - g, self.h())
        d2 = ffa.poly_mod_reduction(g - f, self.h())
        neg_d2 = pa.Polynomial([(5 - c) % 5 for c in d2.coefficients], 5)
        self.assertEqual(d1.coefficients, neg_d2.coefficients)

if __name__ == "__main__":
    unittest.main()
