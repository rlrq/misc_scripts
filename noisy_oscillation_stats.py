#!/usr/bin/python3

## reads a series of values and extracts the following information:
## - 1/2 oscillation period (per 1/2 cycle, defined as the # of values that consecutively change in the same direction)
## - 1/2 oscillation amplitude (per 1/2 cycle, defined as difference between the last & first value within each 1/2 cycle period)
## - 1/2 oscillation min (per 1/2 cycle, defined as the smaller value bounding a 1/2 cycle)
## - 1/2 oscillation max (per 1/2 cycle, defined as the larger value bounding a 1/2 cycle)
## - 1/2 oscillation direction (per 1/2 cycle, direction of change of consecutive values; +/-)
## - smoothed 1/2 oscillation period (per 1/2 cycle, defined as 1/2 oscillation period, but allowing for changes in opposite direction less than a predefined magnitude)
## - smoothed 1/2 oscillation amplitude (per 1/2 cycle, defined as the difference between the largest and smallest value within as smoothed 1/2 oscillation period)
## - smoothed 1/2 oscillation min (per 1/2 cycle, smallest value in a smoothed 1/2 oscillation period)
## - smoothed 1/2 oscillation max (per 1/2 cycle, largest value in a smoothed 1/2 oscillation period)
## - smoothed 1/2 oscillation direction (per 1/2 cycle, overall direction of change smoothed 1/2 oscillation period; +/-)

