"""Wrapper class for the C++ based SMS solver
"""

import ctypes as ct
import sys
import pysms

smslib = ct.CDLL("libsms.so")

# load functions into aliases
sms_create_solver = smslib.create_solver
sms_add_literal = smslib.add_literal
sms_next_solution = smslib.next_solution
sms_destroy_solver = smslib.destroy_solver

# specify function signatures
sms_create_solver.argtypes = [ct.c_int]
sms_create_solver.restype = ct.c_void_p
sms_next_solution.argtypes = [ct.c_void_p]
sms_next_solution.restype = ct.POINTER(ct.c_int)
sms_add_literal.argtypes = [ct.c_void_p, ct.c_int]
sms_add_literal.restype = None
sms_destroy_solver.argtypes = [ct.c_void_p]
sms_destroy_solver.restype = None


class Solver:
    def __init__(self, vertices : int = 2, clauses=[]):
        self.sms_solver = sms_create_solver(max(vertices, 2)) # <2 vertices is an error
        for clause in clauses:
            self.addClause(clause)

    def __del__(self):
        sms_destroy_solver(self.sms_solver)

    @classmethod
    def fromBuilder(cls, builder : pysms.graph_builder.GraphEncodingBuilder):
        return cls(builder.n, builder)

    def addClause(self, clause):
        for lit in clause:
            sms_add_literal(self.sms_solver, lit)
        sms_add_literal(self.sms_solver, 0)

    def __iter__(self):
        return self

    def __next__(self):
        sol = sms_next_solution(self.sms_solver)
        if sol:
            pass
            m = sol[0]
            G = []
            for i in range(m):
                G.append((sol[2*i+1], sol[2*i+2]))
            return G
        else:
            raise StopIteration()

if __name__ == "__main__":
    # as an example, enumerate all graphs with n vertices when run standalone
    n = int(sys.argv[1])
    t = 0
    solver = Solver(n)
    for g in solver:
        print(g)
        t += 1
    del solver
    print(f"there are {t} non-isomorphic graphs on {n} vertices")
