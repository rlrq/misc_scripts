#!/usr/bin/python3

## reads a series of values and extracts the following information:
## - 1/2 oscillation period (per 1/2 cycle, defined as the # of values that consecutively change in the same direction)
## - 1/2 oscillation amplitude (per 1/2 cycle, defined as difference between the last & first value within each 1/2 cycle period)
## - 1/2 oscillation min (per 1/2 cycle, defined as the smaller value bounding a 1/2 cycle)
## - 1/2 oscillation max (per 1/2 cycle, defined as the larger value bounding a 1/2 cycle)
## - 1/2 oscillation direction (per 1/2 cycle, direction of change of consecutive values; +/-)
## - 1/2 oscillation rate (per 1/2 cycle, 1/2 oscillation amplitude divided by 1/2 oscillation period)
## - smoothed 1/2 oscillation period (per 1/2 cycle, defined as 1/2 oscillation period, but allowing for changes in opposite direction less than a predefined magnitude)
## - smoothed 1/2 oscillation amplitude (per 1/2 cycle, defined as the difference between the largest and smallest value within as smoothed 1/2 oscillation period)
## - smoothed 1/2 oscillation min (per 1/2 cycle, smallest value in a smoothed 1/2 oscillation period)
## - smoothed 1/2 oscillation max (per 1/2 cycle, largest value in a smoothed 1/2 oscillation period)
## - smoothed 1/2 oscillation direction (per 1/2 cycle, overall direction of change smoothed 1/2 oscillation period; +/-)
## - smoothed 1/2 oscillation rate (per 1/2 cycle, smoothed 1/2 oscillation amplitude divided by smoothed 1/2 oscillation period)
## - offset 1/2 oscillation period (per 1/2 cycle, defined as the # of values that are consecutively on the same side of the average value within a predefined window size that includes the first value of this period)
## - offset 1 oscillation amplitude (per 1 cycle, defined as the difference between the largest and smallest value between 2 consecutive offset 1/2 oscillation periods)
## - offset 1/2 oscillation min (per 1/2 cycle, smallest value in a offset 1/2 oscillation period)
## - offset 1/2 oscillation max (per 1/2 cycle, largest value in a offset 1/2 oscillation period)
## - offset 1/2 oscillation direction (per 1/2 cycle, whether the offset 1/2 oscillation period is above or below average; +/-)
## - offset 1 oscillation rate (per 1 cycle, offset 1 oscillation amplitude divided by 2 of the sum of the contributing offset 1/2 oscillation periods)


#################
##  TEST DATA  ##
#################

## input file, with header, where each column is one measured variable, and each row is a different timepoint (in order)
f_input = "/media/HDD3/rachelle/scd/results/oscillation/test.variant-freq.subset.txt"

## parse input file
with open(f_input, 'r') as f:
    header_inpt = f.readline().rstrip().split('\t')
    raw_inpt = [line.rstrip().split('\t') for line in f]

## function to subset column
def subset_column(colname, header, data):
    col_i = header.index(colname)
    return [float(entry[col_i]) for entry in data]

## get P2 data
var_dat = subset_column("P2", header_inpt, raw_inpt)


#################
##  FUNCTIONS  ##
#################

def neg_zero_pos(val):
    return 0 if val == 0 else (1 if val > 0 else -1)


#######################
##  STATS FUNCTIONS  ##
#######################

########  1/2 OSCILLATION  ########

## 1/2 oscillation period (per 1/2 cycle, defined as the # of values that consecutively change in the same direction)
## 1/2 oscillation direction (per 1/2 cycle, direction of change of consecutive values; +/-)
def half_periods_direction(data):
    periods = [] ## tracks oscillation period length
    periods_dir = [] ## tracks oscillation period direction
    curr_period_len = 1
    last_direction = neg_zero_pos(data[1] - data[0])
    for i, v in enumerate(data[2:]):
        i = i+2
        curr_direction = neg_zero_pos(data[i] - data[i-1])
        if curr_direction != last_direction:
            periods.append(curr_period_len)
            periods_dir.append(last_direction)
            curr_period_len = 1
            last_direction = curr_direction
        else:
            curr_period_len += 1
    periods.append(curr_period_len)
    periods_dir.append(curr_direction)
    return periods, periods_dir

## 1/2 oscillation amplitude (per 1/2 cycle, defined as difference between the last & first value within each 1/2 cycle period)
def half_amplitude(data, periods):
    output = []
    i_start = 0
    for period in periods:
        i_end = i_start + period
        output.append(data[i_end] - data[i_start])
        i_start = i_end
    return output

## 1/2 oscillation min (per 1/2 cycle, defined as the smaller value bounding a 1/2 cycle)
def half_min(data, periods):
    output = []
    i_start = 0
    for period in periods:
        i_end = i_start + period
        output.append(min(data[i_start], data[i_end]))
        i_start = i_end
    return output

## 1/2 oscillation max (per 1/2 cycle, defined as the larger value bounding a 1/2 cycle)
def half_max(data, periods):
    output = []
    i_start = 0
    for period in periods:
        i_end = i_start + period
        output.append(max(data[i_start], data[i_end]))
        i_start = i_end
    return output

## 1/2 oscillation rate (per 1/2 cycle, 1/2 oscillation amplitude divided by 1/2 oscillation period)
def half_rate(periods, amplitudes):
    return [amplitudes[i]/periods[i] for i in range(len(periods))]
