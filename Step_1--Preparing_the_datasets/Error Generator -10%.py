import sys
import random
import math

from Bio import SeqIO


# Check whether the number of system command line arguments is 3 (the program name itself is always included)
if len(sys.argv) != 3:
    # print usage message to "standard error"
	print("Usage:", file=sys.stderr)
	print(f"{sys.argv[0]} <FASTA-input> <FASTA-output>", file=sys.stderr)
    # exit with 1 to indicate "we didn't do the work"
	sys.exit(1)

# Give the command line arguments names
FASTA_input = sys.argv[1]
FASTA_output = sys.argv[2]

#https://stackoverflow.com/questions/25023233/how-to-save-python-screen-output-to-a-text-file
OutputFile = sys.stdout
sys.stdout = open(FASTA_output, "w")

# the "parse" input strategy, returns a collection (a list) of the resulting sequence objects
FASTA_records = list(SeqIO.parse(FASTA_input, "fasta"))


#https://stackoverflow.com/questions/34338788/python-replace-3-random-characters-in-a-string-with-no-duplicates
#https://stackoverflow.com/questions/59763933/what-is-the-difference-between-the-random-choices-and-random-sample-function#:~:text=sample()%20does%20not%20produce,in%20the%20case%20of%20random.
#https://www.geeksforgeeks.org/python-random-sample-function/
#https://www.w3schools.com/python/ref_random_choices.asp
#https://www.w3schools.com/python/ref_random_choice.asp

replacement_chars='ACGT'
for coi5p in FASTA_records[0:]:
    coi5p_dic = dict(enumerate(coi5p.seq))
    valid_positions = [key for key in coi5p_dic if coi5p_dic[key].isalpha()]
    coi5p.seq = list(coi5p.seq)
    ten_percent = math.ceil(0.1 * len(coi5p.seq))
    random_positions = random.choices(valid_positions, k = ten_percent)

    random_chars = random.choices(replacement_chars, k = len(random_positions))
    char_counter = 0 #position of the character that will replace existing character

    for position in random_positions:
        #check if the replacement matches the existing one and generate another one if needed
        while coi5p.seq[position] == random_chars[char_counter]:
            random_chars[char_counter] = random.choice(replacement_chars)
        coi5p.seq[position] = random_chars[char_counter]
        char_counter = char_counter + 1


    print(">"+coi5p.description)
    print("".join(coi5p.seq))

sys.stdout.close()
sys.stdout = OutputFile

