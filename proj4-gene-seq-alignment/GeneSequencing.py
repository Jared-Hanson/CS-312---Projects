#!/usr/bin/python3

# from PyQt5.QtCore import QLineF, QPointF


import math
import numpy as np
import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

DIAGONAL = 1
LEFT = 2
UP = 3

# Basic object stored in each index in the matrix. Holds cost value, direction from previous node, and if it was a sub
class index_item(object):
    def __init__(self, value, direction, sub):
        self.value = value
        self.direction = direction
        self.sub = sub


class GeneSequencing:

    def __init__(self):
        pass

    # Function to determine if an index has the same value in both sequences or if it differs.
    def diff(self, i, j, sequence1, sequence2):
        if sequence1[i - 1] == sequence2[j - 1]:
            return MATCH
        else:
            return SUB
    # Function retreaves the letter associated with a certain index, Handles edge cases cleanly.
    def getLetter(self, i, sequence):
        if i == 0:
            return '-'
        else:
            return sequence[i - 1]

    # This function takes the two sequences and uses the values stored in the matrix to compute the cost (found in
    # the corner index) and the string to present in the GUI.
    def trace_back(self, E, sequence1, sequence2):
        i = len(sequence1)
        j = len(sequence2)
        curNode = E[i][j]
        cost = curNode.value
        sequenceTop = []
        sequenceLeft = []

        while curNode.direction != "None":
            if curNode.direction == LEFT:
                sequenceTop.append(self.getLetter(i, sequence1))
                sequenceLeft.append('-')
                i = i - 1
                curNode = E[i][j]
            elif curNode.direction == UP:
                sequenceLeft.append(self.getLetter(j, sequence2))
                sequenceTop.append('-')
                j = j - 1
                curNode = E[i][j]
            else:
                sequenceTop.append(self.getLetter(i, sequence1))
                sequenceLeft.append(self.getLetter(j, sequence2))
                i = i - 1
                j = j - 1

                curNode = E[i][j]
        str1 = ""
        str2 = ""
        sequenceTop.reverse()
        sequenceLeft.reverse()
        return cost, str1.join(sequenceTop), str2.join(sequenceLeft)

    # Basic universal align algorithm. Implementation basically pulled from the book. Uses previously implemented
    # matrix E. Calls traceback function (used for bod banded and un-banded) to retrieve the aligned strings.
    def find_alignment(self, sequence1, sequence2, E):
        for i in range(0, len(sequence1) + 1):
            E[i][0] = index_item(i * INDEL, LEFT, 0)
        for j in range(1, len(sequence2) + 1):
            E[0][j] = index_item(j * INDEL, UP, 0)
        E[0][0].direction = "None"

        # Double loop with three possible directional cases to consider each time.
        for j in range(1, len(sequence2) + 1):
            for i in range(1, len(sequence1) + 1):
                # Left case:
                if E[i - 1][j].value + INDEL < E[i - 1][j - 1].value + self.diff(i, j, sequence1, sequence2) and \
                        E[i - 1][j].value + INDEL < E[i][j - 1].value + INDEL:
                    E[i][j] = index_item(E[i - 1][j].value + INDEL, LEFT, 0)
                # Right case:
                elif E[i][j - 1].value + INDEL < E[i - 1][j - 1].value + self.diff(i, j, sequence1, sequence2) and E[i][
                    j - 1].value + INDEL < E[i - 1][j].value + INDEL:
                    E[i][j] = index_item(E[i][j - 1].value + INDEL, UP, 0)
                # Diagonal case (check if sub or match)
                else:
                    E[i][j] = index_item(E[i - 1][j - 1].value + self.diff(i, j, sequence1, sequence2), DIAGONAL,
                                         self.diff(i, j, sequence1, sequence2) == SUB)

        return self.trace_back(E, sequence1, sequence2)

    # Bounded algorithm implementation. Uses the same previously initialized matrix.
    # For loops only run for the bounded indexes, edge cases are handled differently. Returns INF if the two
    # sequences are comparatively to large or small for each other.
    def find_alignment_bounded(self, sequence1, sequence2, E):
        if abs(len(sequence1) - len(sequence2)) > MAXINDELS:
            return float('inf'), "No Alignment Possible", "No Alignment Possible"
        for i in range(0, MAXINDELS + 1):
            E[i][0] = index_item(i * INDEL, LEFT, 0)
        for j in range(1, MAXINDELS + 1):
            E[0][j] = index_item(j * INDEL, UP, 0)
        E[0][0].direction = "None"

        # Double loop with three possible directional cases to consider each time.
        for j in range(1, len(sequence2) + 1):
            for i in range(j - MAXINDELS, j + MAXINDELS + 1):
                if 0 < i < len(sequence1) + 1:
                    if E[i - 1][j].value + INDEL < E[i - 1][j - 1].value + self.diff(i, j, sequence1, sequence2) and \
                            E[i - 1][j].value + INDEL < E[i][j - 1].value + INDEL:
                        E[i][j] = index_item(E[i - 1][j].value + INDEL, LEFT, 0)

                    elif E[i][j - 1].value + INDEL < E[i - 1][j - 1].value + self.diff(i, j, sequence1, sequence2) and \
                            E[i][j - 1].value + INDEL < E[i - 1][j].value + INDEL:
                        E[i][j] = index_item(E[i][j - 1].value + INDEL, UP, 0)
                    else:
                        E[i][j] = index_item(E[i - 1][j - 1].value + self.diff(i, j, sequence1, sequence2), DIAGONAL,
                                             self.diff(i, j, sequence1, sequence2) == SUB)

        return self.trace_back(E, sequence1, sequence2)

    # Function to send sequences to correct function, either bounded or unbounded.
    # Sequences come in to function already of "align_length" size
    def align_single(self, sequence1, sequence2, bounded, E):
        if bounded:
            return self.find_alignment_bounded(sequence1, sequence2, E)
        else:
            return self.find_alignment(sequence1, sequence2, E)

    # This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
    # handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
    # you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment
    def align(self, sequences, table, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length
        results = []

        # This is where my matrix is initialized. The same matrix will be used throughout.
        basic_index_item = index_item(float('inf'), "uninitialized", 0)
        E = [[basic_index_item for x in range(align_length + 1)] for y in range(align_length + 1)]

        for i in range(len(sequences)):
            jresults = []
            for j in range(len(sequences)):
                if j < i:
                    s = {}
                else:

                    # My single function call returns all three needed pieces of information. It is called using the
                    # truncated strings of Align_Length. Final alignments are also truncated to length 100.
                    tempScore, align1, align2 = self.align_single(sequences[i][:align_length],
                                                                  sequences[j][:align_length], banded, E)
                    score = tempScore
                    alignment1 = align1[:100].format(i + 1,
                                                     len(sequences[i]), align_length, ',BANDED' if banded else '')
                    alignment2 = align2[:100].format(j + 1,
                                                     len(sequences[j]), align_length, ',BANDED' if banded else '')

                    s = {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
                    table.item(i, j).setText('{}'.format(int(score) if score != math.inf else score))
                    table.repaint()
                jresults.append(s)
            results.append(jresults)
        return results
