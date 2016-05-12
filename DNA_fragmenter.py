__author__ = 'Kevin'

# NGS technology reads produce shorter reads anywhere from 25-500bp


from tkinter import *
from tkinter import filedialog
from random import randint

genome = ''
outputfile = "FRAGMENTED.txt"

root = Tk()
root.resizable(0, 0)

flength = int(input('Initial fragmenting length: '))

ftimes = int(input('How many times do you wish to fragment: '))

flower = int(input('Lower bound: '))

fupper = int(input('Upper bound: '))


def open_genome():
    filename = filedialog.askopenfilename()
    file = open(filename)
    x = file.read()
    return x

def replace_non_bases(sequence):
    sequence = sequence.replace('\n', '').replace(' ', '')
    print(sequence)
    x = []
    for i in sequence:
        if not i.isdigit():
            x.append(i)
    sequence = ''.join(x)
    return sequence


def fragment_DNA(string, length):
    list = []
    for j in range(0, len(string), length):
        list.append(string[0+j:length+j])
    return list

def run_fragment():
    global flength, ftimes, flower, fupper

    genome = open_genome()

    genome = replace_non_bases(genome)

    print(genome)
    genome_length = (len(genome))
    print(flength)
    fhalf = flength/2


    genome1 = fragment_DNA(genome[int(fhalf):], flength)
    genome2 = list(fragment_DNA(genome, flength))
    genome2.extend(genome1)

    # starting off at half.
    # randomize fragmenting
    for y in range(0, ftimes):
        genome3 = list(fragment_DNA(genome, randint(flower, fupper)))
        print(genome3)
        genome2.extend(genome3)

    with open(outputfile, 'w') as file:
        for item in genome2:
            file.write("%s\n" % item)

    x = 1

    with open(outputfile, 'w') as file:
        for item in genome2:
            file.write("%s %s\n" % (x, item))
            x += 1

root.assemble_button = Button(root, text="Fragment", command=lambda: run_fragment())
root.assemble_button.pack(side=BOTTOM)

root.mainloop()

