

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
    def __init__(self, seq):
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
        

    def _construct_polynomial_formula(self):
        '''
        Will only be called if self.solvable == True
        Evaluates the Newton's forward difference formula:
        sum[for k in range(0, order)] of: (kth_diff[0th elem] * something)
        something = ((x * (x-1)) * ((x * (x-1) * (x - 2)) * ... * (x - k + 1)) / factorial(k)
        The formula is useful in calculating the n'th term directly
            without having to calculate the terms in the middle
        Returns a callable that can be called with the required nth term
        '''
        # the diffs
        diffs = np.array([d[0] for d in self.diff_order.values()])
        factorials = factorial( np.arange(self.order + 1) )
##        def get_x_terms(k):
            #xs = [1, x]
        
        # Get the terms to be subtracted from xs
        x_diff_terms = [0, 0] # These are what gets subtracted from x
        terms = [(1 - i) for i in range(self.order - 1)]  
        x_diff_terms += terms
        x_diff_terms = np.array(x_diff_terms)

        other_terms = diffs / factorials
        
        def func(x, x_diff_terms, other_terms):
            'First construct the diffs then get the cumprod of them'
            xs = [1,] + [x - i for i in x_diff_terms[1:]]
            xs = np.array(xs)
            xs = np.cumprod(xs)
            return int(np.sum(xs * other_terms))

        return partial(func, x_diff_terms=x_diff_terms, other_terms=other_terms)

        

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
    assert solver.get_nth_term(len(seq)) == 256  # FAILS
