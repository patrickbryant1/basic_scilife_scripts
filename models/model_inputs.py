#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import Counter
import random

def split_on_h_group(complete_df, percentage):
    '''Split data into train, valid and test
    '''

    counts = Counter(complete_df['H_group_x'])
    unique_groups = [*counts.keys()]
    #1. Shuffle keys of counted_groups, 2 ensures same random shuffling each time
    unique_groups = np.array(unique_groups)
    random.Random(2).shuffle(unique_groups)

    train_groups = []
    valid_groups = []
    test_groups = []

    train_c = 0
    valid_c = 0
    test_c = 0
    total_entries = len(complete_df['H_group_x'])
    for i in range(0, len(unique_groups)):
        if train_c < total_entries*percentage:
            train_c += counts[unique_groups[i]]
            train_groups.append(unique_groups[i])

        if train_c >= total_entries*percentage:
            if valid_c <= total_entries*((1-percentage)/2):
                valid_c += counts[unique_groups[i]]
                valid_groups.append(unique_groups[i])
            else:
                test_c += counts[unique_groups[i]]
                test_groups.append(unique_groups[i])


    return train_groups, valid_groups, test_groups
