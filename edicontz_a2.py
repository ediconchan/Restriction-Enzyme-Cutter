"""
Author : Edicon Chan
Date   : 2023-11-14
Assignment2 : Cut sequences with restriction enzymes
Linting: Pylint (9.90/10), flake8
"""

from sys import argv


def file_handler(fasta_file: str, enzyme_file: str) -> tuple:
    """"Takes a fasta file and restriction enzyme file as input
    and returns extracted strings from the files"""

    # Open nucleotide fasta files and extract the entire nucleotide string
    with open(fasta_file, 'r', encoding='utf8') as nucleotide_file:
        sequence = ""
        seq_name = ""
        for line in nucleotide_file:
            # If it starts with a caret, save the header name
            if line.startswith(">"):
                seq_name = line[1:]
            # Append each line into one long nucldeotide string
            else:
                sequence = sequence + line.rstrip()

    # Open enzyme file and extract the names and patterns
    with open(enzyme_file, 'r', encoding='utf8') as re_file:
        enzyme_list = re_file.readlines()
        # Gets rid of "\n" in the enzyme file
        for line in range(len(enzyme_list)):
            enzyme_list[line] = enzyme_list[line].rstrip()

    return sequence, enzyme_list, seq_name


def pattern_searcher(file: str, pattern: str) -> list:
    '''Searches a string with nucleotides for a given pattern
    and returns the index of every single occurance of the pattern'''

    index_list = []
    pattern_length = len(pattern)
    for nuc_index in range(len(file) - pattern_length + 1):
        # Search the file with a sliding window the size of the given pattern
        if file[nuc_index: nuc_index + pattern_length] == pattern:
            # If sliding window matches the pattern, record index
            index_list.append(nuc_index)
    return index_list


def re_site_index(start_idx: int, re_string: str) -> int:
    """Given a restriction site, and its index in the fasta file,
    returns the index of the nucleotide after the cutting site"""

    return start_idx + re_string.find("^")


def fragments(file: str, index_list: list, enzyme_site: str) -> list:
    """Given a list of indexes where the restriction sites are located,
    and a restriction site cutting pattern, return a list of
    fragments and index value obtained from cutting the nucleotide string"""

    starting_index = 0
    fragment_list = []
    
    # For every restriction site index in the list
    for rs in index_list:
        # Cut sequence according to index sites provided and
        # update starting index for next fragment in the sequence
        splice_index = re_site_index(rs, enzyme_site)
        fragment = file[starting_index: splice_index]
        fragment_list.append(fragment)
        starting_index = splice_index

    # Append last fragment to the list
    fragment = file[starting_index:]
    fragment_list.append(fragment)

    return fragment_list


def fragment_list_organizer(fragment_list: list, group_size: int) -> list:
    """Given a list of fragments and the requested group size output
    organize a fragment list containing a list of subdivided fragment groups"""

    final_frag_list = []
    for frag in fragment_list:
        final_frag_list.append(fragment_group_organizer(frag, group_size))
    return final_frag_list


def fragment_group_organizer(frag: str, group_size: int) -> list:
    """Create a list consisting of the fragment divided into group_size"""

    group_list = []
    for i in range(0, len(frag), group_size):
        group_list.append(frag[i: i + group_size])
    return group_list


def enzyme_cutter(enzyme_list: list, sequence: str,
                  dash_separator: str, group_size: int,
                  row_length: int) -> None:
    """Given a list of restriction enzymes, determine how many
    cutting sites are in the given sequence and output
    details to console"""

    for enzyme in enzyme_list:

        # Extract the name and the restriction site from the list
        enzyme_name, enzyme_site = enzyme.split(";")
        enzyme_pattern = enzyme_site.replace("^", "")

        # Determine how many cutting sites are found with
        # the given restriction enzyme pattern
        rs_list = pattern_searcher(sequence, enzyme_pattern)
        cutting_site_num = len(rs_list)

        if cutting_site_num == 0:
            # Print that no cutting sites were found
            no_site_output(enzyme_name, dash_separator)
        else:
            if cutting_site_num == 1:
                # Print that 1 cutting site found
                one_site_output(enzyme_name, enzyme_site)
            else:
                # Print that X cutting sites were found
                multi_site_output(enzyme_name, enzyme_site, cutting_site_num)

            # Create a list of fragments cut at each restriction site
            frags = fragments(sequence, rs_list, enzyme_site)
            fragment_list = fragment_list_organizer(frags, group_size)

            # Print out the fragments obtained
            fragment_output(fragment_list, group_size,
                            row_length, dash_separator)


def fragment_output(fragment_list: list, group_size: int,
                    row_length: int, dash_separator: str) -> None:
    """Outputs the entire nucleotide code divided into fragments.
    Output is formatted so that length is displayed for each fragment
    and the index number of each row is updated. Maximum length
    of row output and group size of nucleotides can be adjusted"""

    # Keeps track of the index at the start of each row
    nuc_index = 0
    # Counter keeps track of nucleotides in total
    nuc_counter = 0

    for frag in fragment_list:
        # For each fragment, determine and print the length
        fragment_length = (len(frag) - 1) * group_size + len(frag[-1])
        print(f"Length- {fragment_length}")

        for i in range(0, len(frag), row_length):
            # Divide fragment into chunks and print
            # according to row length and group size
            chunks = frag[i: i + row_length]
            print(f"{nuc_index+1}\t" + " ".join(chunks))
            nuc_index += len(chunks) * group_size

        nuc_index = fragment_length + nuc_counter
        nuc_counter += fragment_length
    print("\n")
    print(dash_separator)


def head_output(fasta_file: str, enzyme_file: str,
                seq_name: str, seq_num: int,
                dash_separator: str) -> None:
    """Print file names used and sequence name and base count to console"""

    print("\n")
    print(f"Restriction enzyme analysis of sequence from file {fasta_file}." +
          f"\nCutting with enzymes found in file {enzyme_file}")
    print(dash_separator)
    print(f"Sequence name: {seq_name}" +
          f"Sequence is {seq_num} bases long.")
    print(dash_separator)


def no_site_output(enzyme_name, dash_separator) -> None:
    """Print that no cutting sites were detected in sequence"""

    print(f"There are no sites for {enzyme_name}.\n")
    print(dash_separator)


def one_site_output(enzyme_name, enzyme_site) -> None:
    """Same as multi_site_output(), corrected for
    singular/plural grammar output to console"""

    print("There is 1 cutting site for " +
          f"{enzyme_name}, cutting at {enzyme_site}")
    print("There are 2 fragments:\n")


def multi_site_output(enzyme_name, enzyme_site, cutting_site_num) -> None:
    """Print number of cutting sites, fragments,
    as well as the enzyme details to console"""

    print(f"There are {cutting_site_num} cutting sites for " +
          f"{enzyme_name}, cutting at {enzyme_site}")
    print(f"There are {cutting_site_num + 1} fragments:\n")


def main():
    """Main function that simulates restriction enzyme digest of sequences"""

    # Input file variables
    fasta_file = argv[1]
    enzyme_file = argv[2]

    # Output format variables
    group_size = 10
    row_length = 60//group_size
    dash_separator = "-"*70

    # Extract sequence, enzymes and names from input files
    sequence, enzyme_list, seq_name = file_handler(fasta_file, enzyme_file)
    seq_num = len(sequence)

    # Output information obtained from the fasta and enzyme files
    head_output(fasta_file, enzyme_file,
                seq_name, seq_num, dash_separator)

    # Simulate restriction digest of sequence with
    # different enzymes and output results
    enzyme_cutter(enzyme_list, sequence,
                  dash_separator, group_size, row_length)


if __name__ == "__main__":
    main()
