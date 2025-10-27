import sys
from typing import List, Dict, Iterable, Tuple


# Please do not remove package declarations because these are used by the autograder.

# Insert your suffix_array function here, along with any subroutines you need
def suffix_array(text: str) -> List[int]:
    """
    Generate the suffix array for the given text.
    """
    suffix_array: List[Tuple[str, int]] = list()

    for i in range(len(text)):
        suffix_tuple: Tuple[str, int] = (text[i:], i)
        suffix_array.append(suffix_tuple)

    suffix_array.sort()
    
    suffix_indexes: List[int] = list()
    
    for i in range(len(suffix_array)):
        suffix_indexes.append(suffix_array[i][1])

    return suffix_indexes