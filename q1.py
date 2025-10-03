"""
sudoku_solver.py

Implement the function solve_sudoku(grid: List[List[int]]) -> List[List[int]] using a SAT solver from PySAT.
"""

from pysat.formula import CNF
from pysat.solvers import Solver
from typing import List

def solve_sudoku(grid: List[List[int]]) -> List[List[int]]:
    """Solves a Sudoku puzzle using a SAT solver. Input is a 2D grid with 0s for blanks."""
    
    #This P is a propositional variable which encodes our logic, it gives a number which uniquely identifies r, c, d and hence the SAT Solver is able to interpret the proposition for that particular as true if the number exists in the given row and column and vice versa. r is denoting rows which go from 0 to n-1. c is denoting columns which go from 0 to n-1. d is the value which if present in the cell in column c and row r makes the proposition for P(r,c,d) true. 
    #The uniqueness of number encoded for specific r, c, d can be seen by uniqueness similar encoding a number in base 9.
    def P(r, c, d):
        return r * 81 + c * 9 + d
    
    cnf = CNF()
    
    #This below nested "for" loop structure ensures that for a given cell, no two digits are simultaneously present in a single cell.
    for r in range(9):
        for c in range(9):
            for d1 in range(1, 10):
                for d2 in range(d1+1, 10):
                    cnf.append([-P(r, c, d1), -P(r, c, d2)])

    #This below nested "for" loop structure ensures that for a given cell, at least one digit is present in the cell.
    for r in range(9):
        for c in range(9):
            cnf.append([P(r, c, d) for d in range(1, 10)])

    #This below nested "for" loop structure ensures that in all the rows, no digit is present multiple times in a row. 
    for r in range(9):
        for d in range(1, 10):
            for c1 in range(9):
                for c2 in range(c1 + 1, 9):
                    cnf.append([-P(r, c1, d), -P(r, c2, d)])

    #This below nested "for" loop structure ensures that in all the rows, each unique digit is present in at least one column in the row.
    for r in range(9):
        for d in range(1, 10):
            cnf.append([P(r, c, d) for c in range(9)])

    #This below nested "for" loop structure ensures that in all the columns, no digit is present multiple times in a column.
    for c in range(9):
        for d in range(1, 10):
            for r1 in range(9):
                for r2 in range(r1 + 1, 9):
                    cnf.append([-P(r1, c, d), -P(r2, c, d)])
    
    #This below nested "for" loop structure ensures that in all the columns, each unique digit is present in at least one row in the column.
    for c in range(9):
        for d in range(1, 10):
            cnf.append([P(r, c, d) for r in range(9)])

    #This below nested "for" loop structure ensures that in all the subgrids, no digit is present multiple times in the subgrid.
    for br in range(0, 9, 3):
        for bc in range(0, 9, 3):
            cells = [(r, c) for r in range(br, br + 3) for c in range(bc, bc + 3)]
            for d in range(1, 10):
                for i in range(len(cells)):
                    for j in range(i + 1, len(cells)):
                        r1, c1 = cells[i]
                        r2, c2 = cells[j]
                        cnf.append([-P(r1, c1, d), -P(r2, c2, d)])
    
    #This below nested "for" loop structure ensures that in all the subgrids, each unique digit is present in at least one cell in the subgrid.
    for br in range(0, 9, 3):
        for bc in range(0, 9, 3):
            cells = [(r, c) for r in range(br, br + 3) for c in range(bc, bc + 3)]
            for d in range(1, 10):
                cnf.append([P(r, c, d) for (r, c) in cells])

    #This below nested "for" loop structure accomodates the fact that some digits were already present in the grid.
    for r in range(9):
        for c in range(9):
            if grid[r][c] != 0:
                cnf.append([P(r, c, grid[r][c])])

    #Using the SAT solver to solve the CNF clauses generated above.
    solution = [[0] * 9 for _ in range(9)]
    with Solver(name="glucose3") as solver:
        solver.append_formula(cnf.clauses)
        if not solver.solve():
            return grid
        model = solver.get_model()
        true_vars = set(x for x in model if x > 0)
        for r in range(9):
            for c in range(9):
                for d in range(1, 10):
                    if P(r, c, d) in true_vars:
                        solution[r][c] = d 
                        break
    return solution