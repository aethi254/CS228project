"""
Sokoban Solver using SAT (Boilerplate)
--------------------------------------
Instructions:
- Implement encoding of Sokoban into CNF.
- Use PySAT to solve the CNF and extract moves.
- Ensure constraints for player movement, box pushes, and goal conditions.

Grid Encoding:
- 'P' = Player
- 'B' = Box
- 'G' = Goal
- '#' = Wall
- '.' = Empty space
"""

from pysat.formula import CNF
from pysat.solvers import Solver

# Directions for movement
DIRS = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

class SokobanEncoder:
    def __init__(self, grid, T):
        """
        Initialize encoder with grid and time limit.

        Args:
            grid (list[list[str]]): Sokoban grid.
            T (int): Max number of steps allowed.
        """
        self.grid = grid
        self.T = T
        self.N = len(grid)
        self.M = len(grid[0])

        self.goals = []
        self.boxes = []
        self.player_start = None

        # TODO: Parse grid to fill self.goals, self.boxes, self.player_start
        self._parse_grid()

        self.num_boxes = len(self.boxes)
        self.cnf = CNF()

    def _parse_grid(self):
        """Parse grid to find player, boxes, and goals."""
        # TODO: Implement parsing logic
        for i in range(self.N):
            for j in range(self.M):
                k = self.grid[i][j]
                if k == 'P':
                    self.player_start = (i, j)
                elif k == 'B':
                    self.boxes.append((i, j))
                elif k == 'G':
                    self.goals.append((i, j))

    # ---------------- Variable Encoding ----------------
    def var_player(self, x, y, t):
        """
        Variable ID for player at (x, y) at time t.
        """
        # TODO: Implement encoding scheme
        return 1 + t * (self.N * self.M) + x * self.M + y

    def var_box(self, b, x, y, t):
        """
        Variable ID for box b at (x, y) at time t.
        """
        # TODO: Implement encoding scheme
        return (1 + b) * ((self.T + 1) * self.N * self.M) + t * (self.N * self.M) + x * self.M + y + 1

    # ---------------- Encoding Logic ----------------
    def encode(self):
        """
        Build CNF constraints for Sokoban:
        - Initial state
        - Valid moves (player + box pushes)
        - Non-overlapping boxes
        - Goal condition at final timestep
        """
        # TODO: Add constraints for:
        # 1. Initial conditions
        # 2. Player movement
        # 3. Box movement (push rules)
        # 4. Non-overlap constraints
        # 5. Goal conditions
        # 6. Other conditions
        
        # Divide cells into walls and non-walls
        cellswithoutwalls = [(i,j) for i in range(self.N) for j in range(self.M) if self.grid[i][j] != '#']
        cellswithwalls = [(i,j) for i in range(self.N) for j in range(self.M) if self.grid[i][j] == '#']
        
        # 1. Initial conditions

        # Player starts at initial position
        x, y = self.player_start
        self.cnf.append([self.var_player(x, y, 0)])

        # Each box starts at its initial position
        for box, (boxx, boxy) in enumerate(self.boxes):
            self.cnf.append([self.var_box(box, boxx, boxy, 0)])

        # Player must be at exactly one position at time 0 (among valid cells)
        playerpos0 = [self.var_player(i, j, 0) for i,j in cellswithoutwalls]
        # At least one position
        if playerpos0:
            self.cnf.append(playerpos0)
        # At most one position (no two positions simultaneously)
        for a in range(len(playerpos0)):
            for b in range(a + 1, len(playerpos0)):
                self.cnf.append([-playerpos0[a], -playerpos0[b]])
        
        # Each box must be at exactly one position at time 0 (among valid cells)
        for box in range(self.num_boxes):
            boxpos0 = [self.var_box(box, i, j, 0) for i,j in cellswithoutwalls]
            # At least one position
            if boxpos0:
                self.cnf.append(boxpos0)
            # At most one position per box
            for a in range(len(boxpos0)):
                for b in range(a + 1, len(boxpos0)):
                    self.cnf.append([-boxpos0[a], -boxpos0[b]])
        
        # No two boxes can occupy the same cell at time 0
        for (i,j) in cellswithoutwalls:
            cellboxes = [self.var_box(box, i, j, 0) for box in range(self.num_boxes)]
            for a in range(len(cellboxes)):
                for b in range(a + 1, len(cellboxes)):
                    self.cnf.append([-cellboxes[a], -cellboxes[b]])
        
        # Player and boxes cannot occupy the same cell at time 0
        for (i,j) in cellswithoutwalls:
            for b in range(self.num_boxes):
                self.cnf.append([-self.var_player(i, j, 0), -self.var_box(b, i, j, 0)])

        # Nothing can be in wall cells at time 0
        for (i,j) in cellswithwalls:
            # Player cannot be in walls
            self.cnf.append([-self.var_player(i, j, 0)])  
            for box in range(self.num_boxes):
                # Boxes cannot be in walls
                self.cnf.append([-self.var_box(box, i, j, 0)])  

        # 2. Player movement

        for t in range(self.T+1):
            # Player must be at exactly one position at each timestep
            playert = [self.var_player(i, j, t) for (i, j) in cellswithoutwalls]
            # At least one position
            if playert:
                self.cnf.append(playert)
            # At most one position
            for a in range(len(playert)):
                for b in range(a + 1, len(playert)):
                    self.cnf.append([-playert[a], -playert[b]])
        
            # Player cannot be in wall cells
            for (i, j) in cellswithwalls:
                self.cnf.append([-self.var_player(i, j, t)])

            # Player and boxes cannot occupy the same cell at any timestep
            for (i, j) in cellswithoutwalls:
                for b in range(self.num_boxes):
                    self.cnf.append([-self.var_player(i, j, t), -self.var_box(b, i, j, t)])

        # If player is at (i,j) at time t, they must move to an adjacent valid cell at time t+1.
        for t in range(self.T):
            for (i, j) in cellswithoutwalls:
                validnextpos = []
                for d, (dx, dy) in DIRS.items():
                    newi, newj = i + dx, j + dy
                    if (newi, newj) in cellswithoutwalls: 
                        validnextpos.append(self.var_player(newi, newj, t + 1))
                if validnextpos:
                    self.cnf.append([-self.var_player(i, j, t)] + validnextpos)

        for t in range(self.T):
            for (i, j) in cellswithoutwalls:  # Player's current position
                for d, (dx, dy) in DIRS.items():
                    newi, newj = i + dx, j + dy  # Player's next position
                    if (newi, newj) not in cellswithoutwalls:
                        continue
                    
                    for b in range(self.num_boxes):
                        new2i, new2j = newi + dx, newj + dy  # Box's destination if pushed
                        
                        if (new2i, new2j) in cellswithoutwalls:
                            # Push is possible, so we have to ensure that there are no other boxes at destination
                            nototherboxesatpushdest = []
                            for b2 in range(self.num_boxes):
                                if b2 != b:
                                    nototherboxesatpushdest.append(-self.var_box(b2, new2i, new2j, t))
                            
                            # If player moves from (i,j) to (newi,newj) and box b is at (newi,newj), then destination must be free and box moves to (new2i,new2j)
                            clause = ([-self.var_player(i, j, t), -self.var_player(newi, newj, t + 1), 
                                       -self.var_box(b, newi, newj, t)] 
                                      + nototherboxesatpushdest + [self.var_box(b, new2i, new2j, t + 1)])
                            self.cnf.append(clause)
                        
                        else:
                            # Push is impossible, so player cannot move into (newi,newj) if box b is there
                            self.cnf.append([-self.var_player(i, j, t), -self.var_player(newi, newj, t + 1), 
                                           -self.var_box(b, newi, newj, t)])

        # 3. Box movement (push rules)

        for t in range(self.T + 1):
            # Each box must be at exactly one position at each timestep
            for b in range(self.num_boxes):
                boxt = [self.var_box(b, i, j, t) for (i, j) in cellswithoutwalls]
                if boxt:
                    self.cnf.append(boxt)  # At least one position
                # At most one position per box
                for a in range(len(boxt)):
                    for c in range(a + 1, len(boxt)):
                        self.cnf.append([-boxt[a], -boxt[c]])

            # Boxes cannot be in wall cells
            for (i, j) in cellswithwalls:
                for b in range(self.num_boxes):
                    self.cnf.append([-self.var_box(b, i, j, t)])

            # No two boxes can occupy the same cell
            for (i, j) in cellswithoutwalls:
                for b1 in range(self.num_boxes):
                    for b2 in range(b1 + 1, self.num_boxes):
                        self.cnf.append([-self.var_box(b1, i, j, t), -self.var_box(b2, i, j, t)])

        # Box can stay in place or move to adjacent cell
        for t in range(self.T):
            for b in range(self.num_boxes):
                for (i, j) in cellswithoutwalls:
                    validnextposbox = []
                    validnextposbox.append(self.var_box(b, i, j, t + 1))  # Stay in place
                    # Or move to adjacent cells if pushed
                    for d, (dx, dy) in DIRS.items():
                        newi, newj = i + dx, j + dy
                        if (newi, newj) in cellswithoutwalls:
                            validnextposbox.append(self.var_box(b, newi, newj, t + 1))
                    if validnextposbox:
                        # If box is at (i,j) at time t, it must be at one of these positions at t+1
                        self.cnf.append([-self.var_box(b, i, j, t)] + validnextposbox)

        # Box can only move if pushed that is if box moves from (i,j) to (newi,newj), then player must have been at (previ,prevj) at time t and at (i,j) at time t+1
        for t in range(self.T):
            for b in range(self.num_boxes):
                for (i, j) in cellswithoutwalls:
                    for d, (dx, dy) in DIRS.items():
                        newi, newj = i + dx, j + dy
                        if (newi, newj) not in cellswithoutwalls:
                            continue
                        previ, prevj = i - dx, j - dy
                        if (previ, prevj) in cellswithoutwalls:
                            # If box moves from (i,j) to (newi,newj), it must have been pushed, so we have to ensure that the player was at (previ,prevj) at time t and (i,j) at time t+1
                            self.cnf.append([-self.var_box(b, i, j, t), -self.var_box(b, newi, newj, t + 1), 
                                            self.var_player(previ, prevj, t)])
                            self.cnf.append([-self.var_box(b, i, j, t), -self.var_box(b, newi, newj, t + 1), 
                                            self.var_player(i, j, t + 1)])
                        else:
                            # If no push possible, box cannot move in this direction
                            self.cnf.append([-self.var_box(b, i, j, t), -self.var_box(b, newi, newj, t + 1)])
        
        # 5. Goal conditions

        for b in range(self.num_boxes):
            boxatgoals = [self.var_box(b, gi, gj, self.T) for (gi, gj) in self.goals]
            if boxatgoals:
                self.cnf.append(boxatgoals)  # Box must be at some goal

        return self.cnf
        
        # Non-Overlap and other conditions are already handled above

def decode(model, encoder):
    """
    Decode SAT model into list of moves ('U', 'D', 'L', 'R').

    Args:
        model (list[int]): Satisfying assignment from SAT solver.
        encoder (SokobanEncoder): Encoder object with grid info.

    Returns:
        list[str]: Sequence of moves.
    """
    if not model:
        return -1
    
    N, M, T = encoder.N, encoder.M, encoder.T
    cellswithoutwalls = [(i,j) for i in range(N) for j in range(M) if encoder.grid[i][j] != '#']
    
    # Find player positions at each timestep by checking which variables are true
    playerpositions = []
    for t in range(T + 1):
        for (i, j) in cellswithoutwalls:
            varid = encoder.var_player(i, j, t)
            if varid in model:
                playerpositions.append((i, j))
                break
    
    # Convert position sequence to move directions
    moves = []
    for t in range(len(playerpositions) - 1):
        currposition = playerpositions[t]
        nextposition = playerpositions[t + 1]
        
        # Calculate direction of movement of the player
        di, dj = nextposition[0] - currposition[0], nextposition[1] - currposition[1]
        
        if (di, dj) == (-1, 0):
            moves.append('U')
        elif (di, dj) == (1, 0):
            moves.append('D')
        elif (di, dj) == (0, -1):
            moves.append('L')
        elif (di, dj) == (0, 1):
            moves.append('R')
    
    return moves


def solve_sokoban(grid, T):
    """
    DO NOT MODIFY THIS FUNCTION.

    Solve Sokoban using SAT encoding.

    Args:
        grid (list[list[str]]): Sokoban grid.
        T (int): Max number of steps allowed.

    Returns:
        list[str] or "unsat": Move sequence or unsatisfiable.
    """
    encoder = SokobanEncoder(grid, T)
    cnf = encoder.encode()

    with Solver(name='g3') as solver:
        solver.append_formula(cnf)
        if not solver.solve():
            return -1

        model = solver.get_model()
        if not model:
            return -1

        return decode(model, encoder)