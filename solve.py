##
# 2WF90 Algebra for Security -- Software Assignment 2
# Polynomial and Finite Field Arithmetic
# solve.py
#
#
# Group number:
# group_number 
#
# Author names and student IDs:
# Maciej Wojtkowski (1984209)
# Sofiia Larina (1971018))
# author_name_3 (author_student_ID_3)
# author_name_4 (author_student_ID_4)
##

# Import built-in json library for handling input/output 
import json
import poly_arithmetic as pa
import finite_field_arithmetic as ffa

def solve_exercise(exercise_location : str, answer_location : str):
    """
    solves an exercise specified in the file located at exercise_location and
    writes the answer to a file at answer_location. Note: the file at
    answer_location might not exist yet and, hence, might still need to be created.
    """
    
    # Open file at exercise_location for reading.
    with open(exercise_location, "r") as exercise_file:
        # Deserialize JSON exercise data present in exercise_file to corresponding Python exercise data 
        exercise = json.load(exercise_file)
        

    ### Parse and solve ###

    answer = {}
    
    p = exercise["integer_modulus"]
    task = exercise["task"]
    ex_type = exercise["type"]
    # Check type of exercise
    try:
        if exercise["type"] == "polynomial_arithmetic":
            if task in ["addition", "subtraction", "multiplication", "long_division",                "extended_euclidean_algorithm"]:
                f = pa.Polynomial(exercise["f"], p)
                g = pa.Polynomial(exercise["g"], p)
                # Check what task within the polynomial arithmetic tasks we need to perform
                if exercise["task"] == "addition":
                    # Solve polynomial arithmetic addition exercise
                    result = f + g
                    answer["answer"] = result.coefficients
                    pass
                elif exercise["task"] == "subtraction":
                    # Solve polynomial arithmetic subtraction exercise
                    result = f - g
                    answer["answer"] = result.coefficients
                    pass
                elif task == "multiplication":
                    result = f * g
                    answer["answer"] = result.coefficients
                
                elif task == "long_division":
                    q, r = pa.polynomial_LD(f, g)
                    if q is None:
                        answer["answer-q"] = None
                        answer["answer-r"] = None
                    else:
                        answer["answer-q"] = q.coefficients
                        answer["answer-r"] = r.coefficients
                    
                elif task == "extended_euclidean_algorithm":
                    a, b, d = pa.poly_extended_euclidean_algorithm(f, g)
                    answer["answer-a"] = a.coefficients
                    answer["answer-b"] = b.coefficients
                    answer["answer-gcd"] = d.coefficients
                    
            elif task == "irreducibility_check":
                    f = pa.Polynomial(exercise["f"], p)
                    answer["answer"] = pa.poly_irreducibility_check(f)
                
            elif task == "irreducible_element_generation":
                n = exercise["degree"]
                poly = pa.poly_generate_irreducible(p, n)
                answer["answer"] = poly.coefficients

        else: # exercise["type"] == "finite_field_arithmetic"
            
            h = pa.Polynomial(exercise["polynomial_modulus"], p)
            
            if task in ["addition", "subtraction", "multiplication", "division"]:
                f = pa.Polynomial(exercise["f"], p)
                g = pa.Polynomial(exercise["g"], p)
                # Check what task within the finite field arithmetic tasks we need to perform
                if exercise["task"] == "addition":
                    g = pa.Polynomial(exercise["g"], p)
                    result = f + g
                    result = ffa.poly_mod_reduction(result, h)
                    answer["answer"] = result.coefficients
                    
                elif task == "subtraction":
                    g = pa.Polynomial(exercise["g"], p)
                    result = f - g
                    result = ffa.poly_mod_reduction(result, h)
                    answer["answer"] = result.coefficients
                    
                elif task == "multiplication":
                    g = pa.Polynomial(exercise["g"], p)
                    result = ffa.finite_field_multiply(f, g, h)
                    answer["answer"] = result.coefficients
                    
                elif task == "division":
                    g = pa.Polynomial(exercise["g"], p)
                    result = ffa.finite_field_division(f, g, h)
                    if result is None:
                        answer["answer"] = None
                    else:
                        answer["answer"] = result.coefficients
                        
            elif task == "inversion":
                f = pa.Polynomial(exercise["f"], p)
                f_inv = ffa.finite_field_inversion(f,h) #
                if f_inv is None:
                    answer["answer"] = None
                else:
                    answer["answer"] = f_inv.coefficients
            elif task == "primitivity_check":
                f = pa.Polynomial(exercise["f"], p)
                is_prim = ffa.is_primitive(f, h, p)
                answer["answer"] = is_prim
            elif task == "primitive_element_generation":
                n = h.degree()
                prim_elem = ffa.primitive_generation(h, p)
                answer["answer"] = prim_elem.coefficients
            # Solve finite field arithmetic addition exercise
        # et cetera
    except Exception as e:
    # Handle any errors gracefully
        if "answer-q" in answer or "answer-r" in answer:
            answer["answer-q"] = None
            answer["answer-r"] = None
        elif "answer-a" in answer:
            answer["answer-a"] = None
            answer["answer-b"] = None
            answer["answer-gcd"] = None
        else:
            answer["answer"] = None
    # Open file at answer_location for writing, creating the file if it does not exist yet
    # (and overwriting it if it does already exist).
    with open(answer_location, "w") as answer_file:
        # Serialize Python answer data (stored in answer) to JSON answer data and write it to answer_file
        json.dump(answer, answer_file, indent=4)

# You can call your function from here
# Please do not *run* code outside this block
# You can however define other functions or constants
if __name__ == '__main__':
    solve_exercise('Simple/Exercises/exercise0.json', 'Simple/Answers/answer0.json')