#!/usr/bin/python3

## originally from /mnt/chaelab/rachelle/pav_haplotype/src/scd.py

import sys
import math
import pandas as pd
import scipy
import scipy.stats
# from scipy.stats import binom
from functools import wraps

# args = sys.argv
# fname = args[1]
# fout = args[2]
# # fname = "/mnt/chaelab/rachelle/panNLRome/results/MAGIC/MAGIC.mosaic-200.noCan-0.vdw.ct.tsv"
# # fout = "/mnt/chaelab/rachelle/panNLRome/results/MAGIC/MAGIC.mosaic-200.noCan-0.vdw.ct.scd.tsv"
# print("fname:", fname)
# print("fout:", fout)


###############
##  METRICS  ##
###############

## r2 from Equations 1,2,3 of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1560400/
def r2(P_AB, P_A, P_B):
    numerator = (P_AB - (P_A * P_B)) ** 2
    denominator = (P_A * P_B) * (1 - P_A) * (1 - P_B)
    return numerator/denominator

def r2_max(P_A, P_B):
    numerator = P_A * (1 - P_B)
    denominator = P_B * (1 - P_A)
    return numerator/denominator

def r2_min(P_A, P_B):
    numerator = P_A * P_B
    denominator = (1 - P_A) * (1 - P_B)
    return numerator/denominator


#################
##  SCD CLASS  ##
#################

class InvalidCellID(Exception):
    """
    Exception raised for when cell ID is invalid
    """
    def __init__(self, cell, *args, **kwargs):
        self.cell = cell
        super().__init__(f"Invalid cell ID: {self.cell}")

class SCD():
    def __init__(self, Cd, Ce, focal_cell):
        ## Cd and Ce are instance methods
        self._Cd = Cd
        self._Ce = Ce
        self.focal_cell = focal_cell
        self.Cd = min(1, max(0, Cd(focal_cell)))
        self.Ce = min(1, max(0, Ce(focal_cell)))
        self.Ds = None
        self._scd()
    @property
    def Cd_function(self):
        return self._Cd.__name__
    @property
    def Ce_function(self):
        return self._Ce.__name__
    def _scd(self):
        self.Ds = self.Cd + ((1-self.Cd) * self.Ce)
        return

## to use this decorator, the first argument MUST be a cell ID
def check_valid_cell(method):
    @wraps(method)
    def _impl(self, cell, *method_args, **method_kwargs):
        self.check_cell(cell)
        return method(self, cell, *method_args, **method_kwargs)
    return _impl

class ContingencyTable():
    def __init__(self):
        self._N = "blank"
    def __iter__(self):
        """
        Generator of all values in contingency table
        """
        raise NotImplementedError("Subclass needs to define this")
    def _calculate_N(self):
        """
        Returns total number of data points in contingency table
        """
        raise NotImplementedError("Subclass needs to define this")
    def is_valid_cell(self, cell):
        """
        Returns bool of whether cell ID is valid
        """
        raise NotImplementedError("Subclass needs to define this")
    def check_cell(self, cell):
        """
        Raises InvalidCellID if cell ID is invalid
        """
        if not self.is_valid_cell(cell):
            raise InvalidCellID(cell)
    def cell_IDs(self):
        """
        Returns generator of all cell IDs
        """
        raise NotImplementedError("Subclass needs to define this")
    @property
    def N(self):
        if self._N is "blank":
            self._N = self._calculate_N()
        return self._N
    def min_cells(self):
        """
        Returns list of cell IDs that have the smallest value in the contingency table
        """
        raise NotImplementedError("Subclass needs to define this")
    @check_valid_cell
    def opposite_cell(self, cell):
        """
        Takes a cell ID and returns the cell ID of the cell opposite it
        (up to subclasses to define what "opposite" is).
        Raises InvalidCellID if cell ID is invalid.
        """
        raise NotImplementedError("Subclass needs to define this")
    @check_valid_cell
    def contributing_counts(self, cell):
        """
        Takes a cell ID and outputs list of counts (int) contributing to the cell.
        Raises InvalidCellID if cell ID is invalid.
        """
        raise NotImplementedError("Subclass needs to define this")
    @check_valid_cell
    def contributing_frequencies(self, cell):
        """
        Takes a cell ID and outputs list of frequencies (float) contributing to the cell.
        Raises InvalidCellID if cell ID is invalid.
        """
        return [count/self.N for count in self.contributing_counts(cell)]
    @check_valid_cell
    def cell_count(self, cell):
        """
        Takes a cell ID and outputs count frequency (int) of cell.
        Raises InvalidCellID if cell ID is invalid.
        """
        raise NotImplementedError("Subclass needs to define this")
    @check_valid_cell
    def cell_frequency(self, cell):
        """
        Takes a cell ID and outputs count frequency (int) of cell.
        Raises InvalidCellID if cell ID is invalid.
        """
        return self.cell_count(cell)/self.N
    ####################
    ##  MISC METRICS  ##
    ##    functions   ##
    ####################
    def r2(self):
        """
        Returns NaN if r2 returns ZeroDivisionError else r2
        """
        cell = next(self.cell_IDs())
        try:
            return r2(self.cell_frequency(cell), *self.contributing_frequencies(cell))
        except ZeroDivisionError:
            return float("NaN")
    def binom(self, cell):
        expected_freq = scipy.prod(self.contributing_frequencies(cell))
        true_binom = scipy.stats.binom.cdf(
            self.cell_count(cell), self.N, expected_freq)
        return true_binom
    #####################
    ##  SCD functions  ##
    #####################
    @check_valid_cell
    def scd(self, cell, Cd = "Cd_binom", Ce = "Ce_pielou"):
        """
        Returns SCD object (that calculates single-cell depletion) for a given cell.
        """
        return SCD(getattr(self, Cd), getattr(self, Ce), cell)
        # Cd_val = getattr(self, Cd)(cell)
        # Ce_val = getattr(self, Ce)(cell)
        # return Cd_val + ((1-Cd_val) * Ce_val)
    @check_valid_cell
    def Cd_binom(self, cell):
        """
        Calculates Cd using binom(observed focal freq)/binom(expected focal freq).
        Raises InvalidCellID if cell ID is invalid.
        """
        expected_freq = scipy.prod(self.contributing_frequencies(cell))
        # true_binom = scipy.stats.binom.cdf(
        #     self.cell_count(cell), self.N, expected_freq)
        true_binom = self.binom(cell)
        expected_binom = scipy.stats.binom.cdf(
            (self.N * expected_freq), self.N, expected_freq)
        return true_binom/expected_binom
    @check_valid_cell
    def Ce_pielou(self, cell, minimum_val = 10**-5):
        """
        Calculates Ce using 1 - Pielou's evenness index.
        Raises InvalidCellID if cell ID is invalid.
        (Adds 'minimum_val' to each value so that no cell has a value of 0 as 0 values are ignored by
        Shannon's diversity index [which is used by Pielou's evenness index])
        """
        ## Wikipedia: When all types in the dataset of interest are equally common, all pi values equal 1 / R, and the Shannon index hence takes the value ln(R).
        all_vals = [i for i in self]
        ## remove one cell with same value as focal cell
        all_vals.remove(self.cell_count(cell))
        all_vals = [val + minimum_val for val in all_vals]
        ## divide Shannon's diversity index (Hs, entropy) by theoretical max Hs (ln(R))
        evenness = 1 - (scipy.stats.entropy(all_vals) / scipy.log(len(all_vals)))
        return evenness
    @check_valid_cell
    def Ce_opposite_extreme(self, cell):
        """
        Checks how maxed out the opposite cell is (the more maxed out, the more significant the score)
        """
        opp_cell = self.opposite_cell(cell)
        opp_cell_expected_freq = scipy.prod(self.contributing_frequencies(opp_cell))
        opp_cell_binom = scipy.stats.binom.cdf(
            self.cell_count(opp_cell), self.N, opp_cell_expected_freq)
        ## we do this because we calculated the cdf and what we really want is the
        ## area under the right side of the observed frequency of the cdf curve
        ## (i.e. how much greater than expected it is)
        opp_cell_extremeness = 1 - opp_cell_binom
        return opp_cell_extremeness
    @check_valid_cell
    def Ce_opposite_expectedness(self, cell):
        """
        Checks how close to expected value the opposite cell is
        (the closer to expected, the more significant the score)
        (assume focal cell is 0 and work backwards to get expected counts from observed relative freq)
        (if opposite cell's expected frequency is less than half of both )
        """
        contributing_frequencies = self.contributing_frequencies(cell)
        if all(x >= 0.4 for x in contributing_frequencies):
            return 1
        opp_cell_expected_freq = 1 - sum(contributing_frequencies)
        opp_cell_expected_count = self.N * opp_cell_expected_freq
        opp_cell = self.opposite_cell(cell)
        opp_cell_count = self.cell_count(opp_cell)
        opp_cell_freq = self.cell_frequency(opp_cell)
        # ## square difference
        # squared_diff = (opp_cell_expected_freq - opp_cell_freq) ** 2
        # return squared_diff
        # ## binom difference
        # observed_binom = scipy.stats.binom.cdf(
        #     self.cell_count(opp_cell), self.N, opp_cell_expected_freq)
        # expected_binom = scipy.stats.binom.cdf(
        #     (self.N * opp_cell_expected_freq), self.N, opp_cell_expected_freq)
        # return (observed_binom - expected_binom) ** 2
        # return ((opp_cell_count - opp_cell_expected_count) / max(opp_cell_expected_count, opp_cell_count)) ** 2
        ## TO FIX: doesn't give bad score to complete correlation :(
        return ((opp_cell_count - opp_cell_expected_count) / (self.N - self.cell_count(cell))) ** 2
        # return opp_cell_freq + opp_cell_expected_freq*(1 - opp_cell_freq)
    @check_valid_cell
    def Ce_fisher(self, cell):
        """
        Replaces focal cell value with value of opposite cell then calculates Fisher's Exact
        (Effectively eliminates anything that would've been significant for Fisher's Exact)
        (Wait hang on this is just Pielou...)
        """
        pass
    @check_valid_cell
    def Ce_difference_of_squares(self, cell):
        """
        ???
        """
        pass


## 2 dimensional contingency table using list of lists as input
class ContingencyTable2D(ContingencyTable):
    def __init__(self, l):
        ## takes 2 dimensional list
        self.data = l
        super().__init__()
    def __iter__(self):
        for row in self.data:
            for col in row:
                yield col
        return
    def _calculate_N(self):
        return sum(i for i in self)
    @property
    def dim(self):
        """
        Returns dimensions of table in form of tuple(rows, cols)
        """
        rows = len(self.data)
        if rows == 0:
            return (rows, 0)
        else:
            return (rows, len(self.data[0]))
    def is_valid_cell(self, cell):
        if (hasattr(cell, "__iter__")
            and len(cell) == 2
            and all(isinstance(val, int) for val in cell)
            and len(self.data) > cell[0]
            and len(self.data[cell[0]]) > cell[1]):
            return True
        return False
    def cell_IDs(self):
        for i, row in enumerate(self.data):
            for j in range(len(row)):
                yield (i, j)
        return
    @check_valid_cell
    def row(self, cell):
        return cell[0]
    @check_valid_cell
    def col(self, cell):
        return cell[1]
    def min_cells(self):
        """
        Returns list of coordinates of cells with the smallest value in the contingency table
        """
        cells = [(0,0)]
        min_val = self.data[0][0]
        for row_i, row in enumerate(self.data):
            for col_i, col in enumerate(row):
                if col > min_val: continue
                elif col == min_val: cells.append((row_i, col_i))
                else:
                    cells = [(row_i, col_i)]
                    min_val = col
        return cells
    @check_valid_cell
    def opposite_cell(self, cell):
        """
        Takes cell coordinates and outputs coordinates of cell opposite it
        """
        nrow, ncol = self.dim
        return (nrow - cell[0] - 1,
                ncol - cell[1] - 1)
    @check_valid_cell
    def row_count(self, cell):
        return sum(self.data[self.row(cell)])
    @check_valid_cell
    def col_count(self, cell):
        return sum(row[self.col(cell)] for row in self.data)
    @check_valid_cell
    def row_frequency(self, cell):
        return self.row_count(cell)/self.N
    @check_valid_cell
    def col_frequency(self, cell):
        return self.col_count(cell)/self.N
    @check_valid_cell
    def contributing_counts(self, cell):
        """
        Takes cell coordinates and outputs frequencies of categories contributing to it
        """
        return [self.row_count(cell), self.col_count(cell)]
    @check_valid_cell
    def contributing_frequencies(self, cell):
        """
        Takes cell coordinates and outputs frequencies of categories contributing to it
        """
        return super().contributing_frequencies(cell)
    @check_valid_cell
    def cell_count(self, cell):
        """
        Takes a cell coordinates and outputs count frequency (int) of cell
        """
        row_i, col_i = cell
        return self.data[row_i][col_i]
    ####################
    ##  MISC METRICS  ##
    ##    functions   ##
    ####################
    def fisher_exact(self):
        """
        Returns NaN if dimension of contingency table is not 2x2,
        else returns p-value of Fisher's Exact test (scipy.stats.fisher_exact)
        """
        if self.dim != (2, 2): return float("NaN")
        return scipy.stats.fisher_exact(self.data, alternative = "two-sided")[1]
    def chi2_contingency(self):
        """
        Returns NaN if any cell is 0,
        else returns p-value of output from scipy.stats.chi2_contingency
        """
        for val in self:
            if val == 0: return float("NaN")
        return scipy.stats.chi2_contingency(self.data)[1]


## 2x2 contingency table using list of lists as input
class ContingencyTable2x2(ContingencyTable2D):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    ####################
    ##  MISC METRICS  ##
    ##    functions   ##
    ####################
    def count_if_focal_is_0(self, cell, focal_cell):
        """
        Returns NaN if focal cell cannot be 0. Else calculates expected count in desired cell.
        """
        focal_contributing_frequencies = self.contributing_frequencies(focal_cell)
        if sum(focal_contributing_frequencies) > 1:
            return float("NaN")
        ## if cell is opposite focal cell
        if all(self.row(cell) != self.row(focal_cell) and self.col(cell) != self.col(focal_cell)):
            return self.N - sum(self.contributing_counts(focal_cell))
        ## if same row
        elif self.row(cell) == self.row(focal_cell):
            return self.row_count(focal_cell)
        ## if same column
        else:
            return self.col_count(focal_cell)
    @check_valid_cell
    def frequency_if_focal_is_0(self, cell, focal_cell):
        return self.count_if_focal_is_0(cell, focal_cell)/self.N
    ##########################
    ##  MISC HYPOTHETICALS  ##
    ##       functions      ##
    ##########################
    @check_valid_cell
    def ct_if_focal_is_0(self, cell):
        row_i, col_i = self.row(cell), self.col(cell)
        new_data = [[col for col in row] for row in self.data]
        new_data[row_i][col_i] = 0
        new_data[row_i][(col_i+1)%2] = self.row_count(cell)
        new_data[(row_i+1)%2][col_i] = self.col_count(cell)
        new_data[(row_i+1)%2][(col_i+1)%2] = self.N - self.row_count(cell) - self.col_count(cell)
        return self.__class__(new_data)
    @check_valid_cell
    def ct_if_independent(self, cell):
        """
        Value of cell provided will be adjusted under presumption of independence
        """
        row_i, col_i = self.row(cell), self.col(cell)
        new_data = [[col for col in row] for row in self.data]
        ref_row = self.data[(row_i+1) % 2]
        try:
            ref_ratio = ref_row[col_i]/ref_row[(col_i+1) % 2]
            new_data[row_i][col_i] = ref_ratio * self.data[row_i][(col_i+1) % 2]
            return self.__class__(new_data)
        except ZeroDivisionError:
            return None
    #####################
    ##  SCD functions  ##
    #####################
    @check_valid_cell
    def Ce_r2(self, cell):
        """
        Applies r2 to cell opposite focal cell & normalises by r2 of expected count
        (we first assume focal cell is 0 and work backwards from observed total frequency of each
         categorical variable to get the expected values in each remaining cell)
        (then we work 'forward' to infer the categorical variable frequencies where we would obtain
         these expected values if both variables are independent (allowing for the focal cell to be
         whatever value is needed to achieve an outcome that meets this assumption))
        (THEN we calculate r2 for the opposite cell using hypothetical contributing frequencies and
         observed opposite cell count)
        (eliminates complete association)
        (really only works on 2x2 tables AND where contributing frequencies to focal cell are < 0.5)
        """
        ct_if_0 = self.ct_if_focal_is_0(cell)
        ct_if_0_independent = ct_if_0.ct_if_independent(cell)
        if ct_if_0_independent is None:
            # return float("NaN")
            return 1
        opp_cell = self.opposite_cell(cell)
        new_data_for_calc = ct_if_0_independent.data
        new_data_for_calc[self.row(opp_cell)][self.col(opp_cell)] = self.cell_count(opp_cell)
        ct_for_calc = self.__class__(new_data_for_calc)
        return ct_for_calc.r2()



# ## test
# x = ContingencyTable2D([[20, 0],[20, 20]])
# y = x.scd((0,1))
# x.scd(10)


# #################
# ##  FUNCTIONS  ##
# #################
# ## effectively copied from pav_fisher.R into python3

# ## data should be an iterable of ContingencyTable objects
# ## focal_cell should be a function that takes a ContingencyTable object and outputs a cell ID
# def single_cell_depletion_gen(data, progress_increment = 10000,
#                               focal_cell = lambda x:"naccs1y2y",
#                               Cd = "Cd_binom", Ce = "Ce_pielou"):
#     output = []
#     for i, ct in enumerate(data):
#         if i % progress_increment == 0:
#             print(i)
#         output.append(ct.scd(focal_cell(ct)), Cd = Cd, Ce = Ce)
#     return output


# # def single_cell_depletion_smallest(df, contingency = [["naccs1y2n", "naccs1y2y"],["naccs1n2n", "naccs1n2y"]],
# #                                    progress_increment = 10000, focal_cell = "naccs1y2y",
# #                                    Cd = "Cd_binom", Ce = "Ce_pielou"):
# #     def make_cell_val(row_i):
# #         def cell_val(cell):
# #             return df.iloc[row_i][cell]
# #         return cell_val
# #     def make_contingency_table(row_i):
# #         output = []
# #         for row in contingency:
# #             new_row = []
# #             for colname in row:
# #                 new_row.append(df.iloc[row_i][colname])
# #             output.append(new_row)
# #         return output
# #     ct_categories_d = {ct_category: (row_i, col_i) for row_i, row in enumerate(contingency)
# #                        for col_i, ct_category in enumerate(row)}
# #     ct_categories = list(ct_categories_d.keys())
# #     output = []
# #     for i in range(len(df)):
# #         if i % progress_increment == 0:
# #             print(i)
# #         cell_val = make_cell_val(i)
# #         ## get cells with smallest values
# #         smallest_val = cell_val(ct_categories[0])
# #         smallest_cell = {ct_categories[0]}
# #         for ct_category in ct_categories[1:]:
# #             val = cell_val(ct_category)
# #             if (val == smallest_val):
# #                 smallest_cell = smallest_cell.union({ct_category})
# #             elif (val < smallest_val):
# #                 smallest_val = val
# #                 smallest_cell = {ct_category}
# #         ## calculate scd for cell with smallest value
# #         if (len(smallest_cell) == 1):
# #             smallest_cell = smallest_cell.pop()
# #             ct_depleted = reorder_contingency_table_to_depleted_at_B(make_contingency_table_dict(i),
# #                                                                      smallest_cell)
# #             output.append([smallest_cell,
# #                            scd(ct_depleted, smallest_cell),
# #                            Cd(ct_depleted),
# #                            Ce(ct_depleted)])
# #         else:
# #             output.append([','.join(smallest_cell), 'NA', 'NA', 'NA'])
# #         ## TODO: calculate scd for focal cell
# #     return output


# def single_cell_depletion(df, contingency = {'A': "naccs1y2n", 'B': "naccs1y2y",
#                                              'C': "naccs1n2n", 'D': "naccs1n2y"}):
#     def make_cell_val(row_i):
#         def cell_val(cell):
#             return df.iloc[row_i][contingency[cell]]
#         return cell_val
#     def make_contingency_table_dict(row_i):
#         output_d = {}
#         for cell in ('A', 'B', 'C', 'D'):
#             output_d[cell] = df.iloc[row_i][contingency[cell]]
#         return output_d
#     output = []
#     for i in range(len(df)):
#         if i % 10000 == 0:
#             print(i)
#         cell_val = make_cell_val(i)
#         smallest_val = cell_val('A')
#         smallest_cell = {'A'}
#         for cell in ('B', 'C', 'D'):
#             val = cell_val(cell)
#             if (val == smallest_val):
#                 smallest_cell = smallest_cell.union({cell})
#             elif (val < smallest_val):
#                 smallest_val = val
#                 smallest_cell = {cell}
#         ## calculate scd for cell with smallest value
#         if (len(smallest_cell) == 1):
#             smallest_cell = smallest_cell.pop()
#             ct_depleted = reorder_contingency_table_to_depleted_at_B(make_contingency_table_dict(i),
#                                                                      smallest_cell)
#             output.append([smallest_cell,
#                            scd(ct_depleted, 'B'),
#                            Cd(ct_depleted),
#                            Ce(ct_depleted)])
#         else:
#             output.append([','.join(smallest_cell), 'NA', 'NA', 'NA'])
#     return output

# def reorder_contingency_table_to_depleted_at_B(cont_table, depleted_cell):
#     if (depleted_cell == 'B'):
#         return cont_table
#     elif (depleted_cell == 'A'):
#         ## new order should be C, A, D, B
#         cell_order = ('C', 'A', 'D', 'B')
#     # elif (depleted_cell == 'B'):
#     #     cell_order = ('A', 'B', 'C', 'D')
#     elif (depleted_cell == 'C'):
#         cell_order = ('A', 'C', 'B', 'D')
#     elif (depleted_cell == 'D'):
#         cell_order = ('C', 'D', 'A', 'B')
#     else:
#         print(f"Unknown cell specified ('{depleted_cell}'). Must be 'A', 'B', 'C', or 'D'.")
#         print(depleted_cell)
#         return
#     default_order = ('A', 'B', 'C', 'D')
#     return {default_order[i]: cont_table[cell_order[i]]
#             for i in range(len(default_order))}

# ## cont_table should be a dictionary of length 4, with cell values indexe as A (upper left), B (upper right), C (lower left), and D (lower right).
# ## depleted cell should be either "A", "B", "C", or "D" indicating the cell to assess for depletion (this function only really works if the depleted cell is the cell with the smallest value)
# def scd(cont_table, depleted_cell):
#     ## reorder cont_table if depleted_cell is anything but B so we can base all further calculations on the assumption that B is the depleted cell
#     cont_table = reorder_contingency_table_to_depleted_at_B(cont_table, depleted_cell)
#     depletion = Cd(cont_table)
#     evenness = Ce(cont_table)
#     return We(evenness, depletion)

# ## depletion coefficicent, where depleted cell is B
# ## takes a named vector (where names are A, B, C, or D corresponding to contingency table upper left, right, bottom left, right cells)
# def Cd(ct):
#     N = sum(ct.values())
#     ## expected proportion of B (i.e. P(A or B) * P(B or D) )
#     expected_B_proportion = ( (sum([ct['A'], ct['B']])/N) *
#                               (sum([ct['B'], ct['D']])/N) )
#     ## expected value of B (i.e. P(A or B) * P(B or D) * <sum of all cells> )
#     expected_depleted_cell_val = int(N * expected_B_proportion) ## round down to int cuz count data
#     ## depletion score if B is expected value
#     expected_depletion_score = scipy.stats.binom.cdf(expected_depleted_cell_val,
#                                                      N, expected_B_proportion)
#     ## depletion score given actual ct values
#     actual_depletion_score = scipy.stats.binom.cdf(ct['B'], N, expected_B_proportion)
#     ## ratio of actual/expected depletion scores
#     depletion = actual_depletion_score/expected_depletion_score
#     ## if >1: B value is greater than expected given frequency of A+B and B+D
#     return depletion

# ## expectedness coefficient, where depleted cell is B; calculates "species evenness" for cells A, C, D
# ## takes a dict (where indices are A, B, C, or D corresponding to contingency table upper left, right, bottom left, right cells)
# def Ce_pielou(ct):
#     ## Wikipedia: When all types in the dataset of interest are equally common, all pi values equal 1 / R, and the Shannon index hence takes the value ln(R).
#     evenness = 1 - (scipy.stats.entropy([ct[cell] for cell in ('A', 'C', 'D')]) /
#                     math.log(3))
#     return evenness

# ## expectedness coefficient, where depleted cell is B; calculates deviation from expected frequency
# def Ce_binom(ct):
#     contributing_cells = {'A': ('B', 'C'), 'C': ('A', 'D'), 'D': ('B', 'C')}
#     ## expected proportion
#     N = sum(ct.values())
#     N_without_B = sum([ct['A'], ct['C'], ct['D']])
#     expected_proportions = {cell: ( (sum([ct[cell], ct[adj_cells[0]]])/N) *
#                                     (sum([ct[cell], ct[adj_cells[1]]])/N) )
#                             for cell, adj_cells in contributing_cells.items()}
#     ## without B
#     sum_expected_proportions_without_B = sum(list(expected_proportions.values()))
#     expected_proportions_without_B = {cell: exp_prop / sum_expected_proportions_without_B
#                                       for cell, exp_prop in expected_proportions.items()}
#     ## calculate two-sided binom, then make it one-sided and inverse (i.e. area under curve between observed vs peak)
#     binom_tests = {cell: (1 - scipy.stats.binom_test(ct[cell], N_without_B, exp_prop, alternative = "two-sided"))/2
#                    for cell, exp_prop in expected_proportions_without_B.items()}
#     print(expected_proportions)
#     print(expected_proportions_without_B)
#     print(binom_tests)
#     ## combine binom results
#     expectedness = scipy.prod(list(binom_tests.values()))
#     return expectedness
    

# Ce = Ce_pielou

# ## function that weights Ce based on Cd
# def We(ce, cd):
#     return ( (cd*cd) + ((1-cd) * ce) )

# def D_diversity(dat):
#     N = sum(dat)
#     output = 1
#     for count in dat:
#         freq = count/N
#         output *= freq**(-freq)
#     return output

# ###############
# ##  PROCESS  ##
# ###############

# #############
# ##  MAGIC  ##
# #############

# # ## parse
# # df_ct = pd.read_table(fname, sep = '\t', header = 0,
# #                       converters = {"og1": lambda x:str(x),
# #                                     "og2": lambda x:str(x)})

# # ## process
# # df_scd = single_cell_depletion(df_ct, contingency = {'A': "PA", 'B': "PP", 'C': "AA", 'D': "AP"})

# ## write
# def to_int(x):
#     try: return int(x)
#     except: return x

# # print("Writing")
# # with open(fout, "w+") as f:
# #     _sink = f.write('\t'.join(list(df_ct.columns) +
# #                               ["depleted_cell", "scd.depleted", "cd.depleted", "ce.depleted"]) + '\n')
# #     for i in range(len(df_ct)):
# #         if i % 10000 == 0:
# #             print(i)
# #         tmp = df_ct.iloc[i].values.tolist()
# #         tmp = tmp[:2] + list(map(to_int, tmp[2:]))
# #         _sink = f.write('\t'.join(map(str, tmp + df_scd[i])) + '\n')
