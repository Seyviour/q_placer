from bisect import bisect_left
from collections import defaultdict
from functools import lru_cache
import itertools
from tracemalloc import start
import numpy as np
from collections import deque



#During recursion, for each gate, need to know:
# 1. Current grid gate is in
# 2. current location of gate
# 3. pads connected to gate
# 4. 


#During each recursion step, need to build: 
# 1. New adjacency matrix
# 2. New b matrixes
# 3. Need to be able to transfer this data
# 4. Need to know what gates are in grid cell

class QpSolver():
    def __init__(self) -> None:
        #self._netdict = defaultdict(list)
        self.num_gates = None
        self.num_nets = None
        self.adjacency_matrix = None
        self._A_matrix = {}
        self._pad_net_dict = {}
        self._gate_pad_dict = defaultdict(list)
        self._pad_locations = {}
        self._net_gate_dict = defaultdict(list)
        self._grid_length = 100
        self._grid_heigth = 100
        self._gate_locations = {}
        self._gate_grid_span = {}
        self._EPSILON = 0.00000000
        self._row_size = 0
        self._column_size = 0
    

    def read_netlist(self, netlist_file_name):
        """
        Read netlist from input file. 
        This method implicitly updates self._gate_pad_dict, self._net_gate_dict and self._pad_locations.
        This method also calls self._net_dict_to_adjacency_matrix which in turn updates self.adjacency_matrix
        Format for netlist file *should* be in readme

        parameters: 
            net_list_file_name: string|path

        returns:
            None
        """

        with open(netlist_file_name, 'r') as f:
            line1 = next(f)
            num_gates, num_nets = [int(val) for val in line1.strip().split()]
            self.num_gates = num_gates
            self.num_nets = num_nets
            
            for _ in range (num_gates):

                this_line = next(f)
                gate_id, num_nets_connected, *connected_nets = [int(val) for val in this_line.strip().split()]

                for net in connected_nets:
                    self._net_gate_dict[net].append(gate_id)
            
            next_line = next(f)
            num_pads = int(next_line.strip())
            
            for _ in range(num_pads):
                this_line = next(f)
                pad_id, net_connected, x_coord, y_coord = [int(val) for val in this_line.strip().split()]

                self._pad_locations[pad_id] = (x_coord, y_coord)
                
                for gate in self._net_gate_dict[net_connected]:
                    self._gate_pad_dict[gate].append(pad_id)

            self._net_dict_to_adjacency_matrix()

    def _net_dict_to_adjacency_matrix(self):
        """
        update adjacency matrix with read netlist (self._net_gate_dict)

        parameters:
            None
        
        returns:
            None
        """

        self.adjacency_matrix = np.zeros((self.num_gates, self.num_gates), dtype=np.int32)

        for net_id, gate_list in self._net_gate_dict.items():
            for gate1, gate2 in itertools.combinations(gate_list, 2):
                self.adjacency_matrix[gate1-1][gate2-1] = 1
                self.adjacency_matrix[gate2-1][gate1-1] = 1



    def _calculate_b_matrixes(self, span_gates, pseudo_pad_dict):

        """
        build bx and by vectors for a group of gates

        parameters:
            span_gates: (iterable) of gates to build vectors for
            pseudo_pad_dict: Dictionary containing propagated pads and pseudo pads for each gate in span_gates 
        """

        n = len(span_gates)

        bx = np.zeros((n, 1))
        by = np.zeros((n, 1))

        for idx, gate in enumerate(span_gates):
            pads = pseudo_pad_dict[gate]
            num_pads = len(pads)
            bx[idx] =  x if (x := sum(pad[0] for pad in pads)) is not None else 0 
            by[idx] = y if (y := sum(pad[1] for pad in pads)) is not None else 0

        return bx, by

    
    def _propagate_pad_location(self, boundary, pad_location):

        """
        propagate pads/pseudopads to boundary of gridspan

        parameters:
            boundary: Tuple span of grid cells under consideration - (start_x, end_x, start_y, end_y)
            pad_location: True location of pad to be propagated - Tuple(x_coord, y_coord)
        
        returns:
            (pad_x, pad_y) -> pad/pseudo_pad location after propagation
        """
        
        start_x, end_x, start_y, end_y = self._physical_boundary_from_span(boundary)
        pad_x, pad_y = pad_location

        # Handle the case where a gate has been *assigned* to a different grid span, but is still
        # *physically* in the grid_span under consideration 
        if start_x <= pad_x < end_x and start_y <= pad_y < end_y: 
            if end_y - start_y >= end_x - start_x:
                pad_x = end_x - self._EPSILON
            else:
                pad_y = end_y - self._EPSILON

        if pad_x <= start_x: 
            pad_x = start_x + self._EPSILON
        elif pad_x > end_x:
            pad_x = end_x - self._EPSILON
        
        
        if pad_y <= start_y:
            pad_y = start_y + self._EPSILON
        elif pad_y > end_y:
            pad_y = end_y - self._EPSILON

        return (pad_x, pad_y)


    def _build_pad_pseudo_pad_dict(self, grid_span, span_gates, outside_gates):

        """
        Builds a dict of pads and pseudopads assigned to each gate in span_gates

        parameters:
            grid_span: Tuple -> span of grid cells under consideration
            span_gates: List -> list of gates *assigned* to this grid_span
            outside_gates: List -> list of gates not assigned to this grid_span
        
        returns:
            pad_dict: Dictionary -> dict(gate: list of pad/pseudopad locations) for each gate

        """

        start_row, end_row, start_col, end_col = grid_span

        pad_dict = defaultdict(list)

        for gate1 in span_gates:
            for gate2 in outside_gates:
                if self.adjacency_matrix[gate1-1][gate2-1] == 1:
                    gate_location = self._gate_locations[gate2]
                    propagated_location = self._propagate_pad_location(grid_span, gate_location)
                    pad_dict[gate1].append(propagated_location)
        
        for gate in span_gates: 
            for pad_id in self._gate_pad_dict[gate]:
                pad_location = self._pad_locations[pad_id]
                propagated_location = self._propagate_pad_location(grid_span, pad_location)
                pad_dict[gate].append(propagated_location)
        
        return pad_dict


    def _build_adjacency_matrix(self, span_gates):
        """
        Builds the adjacency matrix for a grid_span/list of gates in a grid_span
        Not to be confused with self._net_dict_to_adjacency_matrix

        parameters:
            span_gates: List -> list of gates in gridspan under consideration

        returns:
            adjacency_matrix: numpy array of shape (n, n) where n is the number of gates in span_gates

        """

        n = len(span_gates)
        C = np.zeros((n,n))

        idxs = [n for n in range(n)]
        for gate1_idx, gate2_idx in itertools.combinations(idxs, 2):

            gate1 = span_gates[gate1_idx]
            gate2 = span_gates[gate2_idx]

            if self.adjacency_matrix[gate1-1][gate2-1] == 1: 
                C[gate1_idx][gate2_idx] = 1
                C[gate2_idx][gate1_idx] = 1 
        
        return C

    
    def _calculate_A_matrix(self, span_gates, pad_dict):
        """
        Build the A matrix for a pass of the Quadratic Placer on a list of gates in a gridspan
        under consideration

        parameters: 
            span_gates: List -> list of gates in grid_span under consideration
            pad_dict: Dictionary -> dictionary (gate: list of associated pads) for each gate in span_gates

        returns:
            A: numpy matrix of shape (n, n) where n is the number of items in span_gates
        """

        C = self._build_adjacency_matrix(span_gates)

        A = -C

        for idx, gate in enumerate(span_gates):
            num_pads = len(pad_dict[gate])
            A[idx][idx] = num_pads + sum(C[idx])
    
        return A


    def _build_grid_matrixes(self, gridspan, span_gates, outside_gates):
        """
        Build matrixes and vectors necessary for the Quadratic placer to
        estimate the locations of gates in gridspan

        parameters:
            gridspan: span of grid cells under consideration
            span_gates: List -> iterable type of gates in gridspan
            outside_gates: List -> iterable of gates outside gridspan
        
        returns:
            A, bx, by
        """

        pad_dict = self._build_pad_pseudo_pad_dict(gridspan, span_gates, outside_gates)
        A = self._calculate_A_matrix(span_gates, pad_dict)
        bx, by = self._calculate_b_matrixes(span_gates, pad_dict)

        return A, bx, by


    def repartition(self, gridspan):

        """
        Apply recursive partitioning to gates assigned to this gridspan to find
        better positions for them

        parameters:
            gridspan: Tuple -> span of grid cells under considertion
                gridspan = (start_row, end_row, start_col, end_col)
        
        returns:
            None
        
        side_effects:
            updates the locations of gates assigned to this gridspan i.e. it updates 
            self._gate_locations for gates assigned to this grid_span
        """

        row_start, row_end, col_start, col_end = gridspan
        span_gates, outside_gates = self._find_gates_in_span(gridspan)

        A, bx, by = self._build_grid_matrixes(gridspan, span_gates, outside_gates)

        X = np.linalg.solve(A, bx)
        Y = np.linalg.solve(A, by)

        # print(bx, by)

        a = list(zip(span_gates, X, Y))
        mid = len(span_gates)//2

        if col_end-col_start >= row_end-row_start:
            a.sort(key = lambda x: (x[1], x[2], x[0]))
            split_col = (col_start + col_end)//2

            for gate, x_coord, y_coord in a[:mid]:
                self._gate_grid_span[gate] = (row_start, row_end, col_start, split_col)

            for gate, x_coord, y_coord in a[mid:]:
                self._gate_grid_span[gate] = (row_start, row_end, split_col, col_end)
        
        else:
            a.sort(key = lambda x: (x[2], x[1], x[0]))
            split_row = (row_start + row_end) // 2

            for gate, x_coord, y_coord in a[:mid]:
                self._gate_grid_span[gate] = (row_start, split_row, col_start, col_end)
            
            for gate, x_coord, y_coord in a[mid:]:
                self._gate_grid_span[gate] = (split_row, row_end, col_start, col_end)

        
        for gate, x_coord, y_coord in a: 
            self._gate_locations[gate] = (x_coord[0], y_coord[0])



    def _physical_boundary_from_span(self, grid_span):

        """
        calculate the "physical" boundaries of a gridspan

        parameters:
            grid_span: Tuple -> span of grid cells to consider

        returns:
            grid_boundaries: Tuple -> "physical" coordinates of grid_span boundaries
        """

        row_start, row_end, col_start, col_end = grid_span

        rs = self._row_size
        cs = self._column_size

        return ( col_start*cs, col_end*rs, row_start * rs, row_end * rs)

    
    def _span_split(self, grid_span):

        """
        split a grid span in two:
            if cols_spanned >= rows_spanned: split with columns
            else: split rows
        
        if the span is already a unit span(defines exactly a grid cell), return the grid_span
        
        parameters:
            grid_span: Tuple -> span of grid rows/columns under consideration
        
        returns:
            splits: List -> two new grid_spans if grid_span is not unit, else single-item list of grid_span passed in 
        """
        
        row_start, row_end, col_start, col_end = grid_span
        
        row_span, col_span = self._calculate_span_size(grid_span)

        split = []

        if col_span == row_span == 1:
            split =  [grid_span]
        
        elif col_span >= row_span:
            col_split = (col_start + col_end)//2
            split.append((row_start, row_end, col_start, col_split))
            split.append((row_start, row_end, col_split, col_end))
        
        else:
            row_split = (row_start + row_end)//2
            split.append((row_start, row_split, col_start, col_end))
            split.append((row_split, row_end, col_start, col_end))
        
        return split    

    def _calculate_span_size(self, gridspan):
        """
        calculates the number of grid rows and columns spanned by a gridspan

        parameters: 
            gridspan: yktv. MFJPM
        
        returns:
            spansize: Tuple (no_rows_spanned, no_cols_spanned)
        """
        row_start, row_end, col_start, col_end = gridspan
        return (row_end - row_start, col_end-col_start)


    def _find_gates_in_span(self, gridspan):

        """
        Filters the gates assigned to a gridspan and returns them in a list

        parameters:
            grid_span: Tuple -> span of grid rows/columns under consideration

        returns:
            span_gates: List -> list of gates *assigned* to gridspan
        """

        row_start, row_end, col_start, col_end = gridspan

        span_gates = []
        outside_gates = []

        span = (row_start, row_end, col_start, col_end)

        for gateid in range(1, self.num_gates + 1):
            if self._gate_grid_span[gateid] == span:
                span_gates.append(gateid)
            else:
                outside_gates.append(gateid)

        return span_gates, outside_gates

    def recursive_repartition(self, gridsize, earlystop = None):

        """
        partition the gates into a grid of gridsize

        parameters:
            gridsize: Tuple -> (num_grid_rows, num_grid_cols)

            earlystop: Tuple -> (num_rows_spanned, num_cols_spanned)
                        If earlystop is defined, a gridspan of size early stop is not recursed into
            
            returns:
                None
            
            side_effects:
                updates gage locations in self._gate_locations
        """

        r, c = gridsize

        self._column_size = self._grid_length / c
        self._row_size = self._grid_heigth / r 

        #row_start, row_end, col_start, col_end = boundaries
        boundaries = (0, r, 0, c)

        for gateid in range(1, self.num_gates+1):
            self._gate_grid_span[gateid] = boundaries

        recurse_stack = deque()
        recurse_stack.append(boundaries)

        while recurse_stack:
            this_boundary = recurse_stack.pop()
            span_size = self._calculate_span_size(this_boundary)
            splits = self._span_split(this_boundary)
            if len(splits) == 2: 
                self.repartition(this_boundary)
                if span_size == earlystop:
                    continue
                else:
                    recurse_stack.extendleft(splits)


        #I'm going to assume span_gates is sorted for now. Really there isn't any reason it shouldn't
        #In any case, it should be a problem for self.find_gates_in_span

    def write_locations(self, filename):

        """
        write gate locations to file

        parameters:
            filename: string|path -> file to write gate locations to
        """
        with open(filename, 'w') as f:
            for gateid in range(1, self.num_gates+1):
                x, y = self._gate_locations[gateid]
                f.write(f"{gateid} {x} {y} \n")



