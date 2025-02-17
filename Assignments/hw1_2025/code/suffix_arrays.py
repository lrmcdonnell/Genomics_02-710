from typing import Tuple, List


def find_longest_repeated_substring(seq: str) -> List[str]:
    """
    The function returns the longest common substring for the input string. It allows overalpping
    substrings. E.g. in 'ATATAT' the longest common substring is 'ATAT'
    :param seq: str: input sequence
    :return: list of str: list of unique longest common substring in alphabetical order, if there is
                          a single substring, the list should contain one element .
    """

    return


def find_most_common_substring_of_length_n(seq: str, n: int) -> Tuple[int, List[str]]:
    """
    The function finds the most common substring of a given length n. If there are multiple
    substrings of length n with highest frequency, the function should return all of them. For the
    sequence: `ATATAT` and n=3 the correct resuts: 2, [ATA, TAT].
    :param seq: str: input sequence
    :param n: int: the length of targe substrings
    :return: number of occurances, and list of all substrings in alphabetical order
    """

    return


if __name__ == "__main__":
    with open("hiv.txt") as f:
        seq = f.readline()
    lrs = find_longest_repeated_substring(seq)
    n = 8
    freq, substrings = find_most_common_substring_of_length_n(seq, n)

    print(f"The first longest common substring has {gc_content(lrs[0]):.3f} GC-content.")
    print("The longest common substring(s) is/are:")
    for s in lrs:
        print(s)

    print(f"The most common substring(s) of length {n} occured {freq} times.")

    print(f"The most common substring(s) of length {n} is/are:")
    for s in substrings:
        print(s)
