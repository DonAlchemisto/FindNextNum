

'''
from vid by singing banana
https://www.youtube.com/watch?v=scQ51q_1nhw

or mathologer's vid:
https://www.youtube.com/watch?v=4AuV93LOPcE

Works if pattern is a polynomial
1st we take the diff until we reach constant diff
We can then apply Newton's formula, or just add
diff -> Babbage machine
'''

import numpy as np
from collections import OrderedDict
from functools import partial

factorial = np.vectorize(np.math.factorial)
n_choose_k = lambda n, k: int(factorial(n) / (factorial(k) * factorial(n - k)))
# Not the most efficient


def find_next_num(seq):
    '''The Babbage difference engine way'''
    seq = np.array(seq)
    seq_len = seq.shape[0]
    
    if seq_len < 3:
        raise ValueError('Provided seq must be at least of length 3')

    diff_order = OrderedDict({0: seq})
    found_sol = False
    for diff in range(1, seq_len):
        d = np.diff(diff_order[diff - 1])
        if np.alen(d) < 2:
            break
        
        diff_order[diff] = d
        if np.alen(np.unique(d)) == 1:
            # all constants
            found_sol = True
            break

    if found_sol:
        order = diff
        result = 0
        for i in reversed(range(order + 1)):
            result += diff_order[i][-1]
        return result

    else:
        return



class PolySeq:
    '''
    TODO:
    - A seq of Powers of 2 -> is its own 1st difference, never reaches a constant
      Add this case

    - Use the multiplicative binomial coefficient formula
      instead of n_choose_k for higher efficiency
      
    '''
    def __init__(self, seq):
        '''
        attrs:
          - solvable: bool -> according to whether a constant difference was reached or not

        methods:
          - get_nth_term:
              args: n [int] -> starting at zero
              returns: int

          - get_next_term:
            calls get_nth_term with the next term index
              returns: int
        '''
        # Establish the required difference
        self.seq = np.array(seq)
        seq_len = np.alen(self.seq)
        
        if seq_len < 3:
            raise ValueError('Provided seq must be at least of length 3')

        self.diff_order = {0: seq}
        self.solvable = False
        for diff in range(1, seq_len):
            d = np.diff(self.diff_order[diff - 1])
            if np.alen(d) < 2:
                break
            
            self.diff_order[diff] = d
            if np.alen(np.unique(d)) == 1:
                # all constants
                self.solvable = True
                break

        if self.solvable:
            self.order = diff
            self._get_nth_term = self._construct_polynomial_formula()


    def _get_nth_term(self, n):
        '''Will only be called if self.solvable == False'''
        assert self.solvable == False
        raise ValueError('The provided sequence is unsolvable with this length at least')


    def get_nth_term(self, n):
        try:
            return self.seq[n]
        except IndexError:
            # Need to calculate it
            return self._get_nth_term(n)

    def get_next_term(self):
        return self.get_nth_term(len(self.seq))
        

    def _construct_polynomial_formula(self):
        def poly_formula(x, coeffs):
            xs = np.array([n_choose_k(x, i) for i in range(len(coeffs))])
            return np.sum(coeffs * xs).astype(np.int64)
            
        coeffs = np.fromiter((val[0] for val in self.diff_order.values()), dtype=np.int64)
        # Should simplify if they have a greatest common divisor

        # calculate exponents -> quick and easy way
        return partial(poly_formula, coeffs=coeffs)
        


if __name__ == '__main__':
    # Define some sequences to start with
    # pentagonal nums
    pent_nums = [0, 1, 5, 12, 22] # next should be 35

    # hexagonal nums
    hex_nums = [0, 1, 6, 15, 28] # next should be 45

    # First approach (Babbage engine)
    next_num = find_next_num(pent_nums)

    # Second approach (Newton's)
    solver = PolySeq(pent_nums)
    assert solver.get_nth_term(5) == 35

    solver = PolySeq(hex_nums)
    assert solver.get_nth_term(5) == 45

    
    seq = [1, 2, 4, 8, 16, 31, 57, 99, 163]
    solver = PolySeq(seq)
    assert solver.get_nth_term(len(seq)) == 256
