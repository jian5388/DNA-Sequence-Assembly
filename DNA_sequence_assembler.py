from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
import os

MAX_NUCLEOTIDES = 5000
NOFILE = 2

filename = ''
reads = dict() 
assembled_reads = ""
issequencesaved = NOFILE

root = Tk()
root.resizable(0, 0)


def read_data():
    """
    Reads DNA fragments into a dictionary, in where the sequence number is the key, and the DNA fragment is the value.
    
    Args:
        NONE

    Returns:
        filename - The name of the input file.

    Raises:
        KeyError: Raises an exception.
    """
    global filename
    global issequencesaved

    filename = filedialog.askopenfilename()
    try:
        if filename.endswith('.txt'):
            if os.stat(filename).st_size == 0:
                filename = ''
                issequencesaved = NOFILE
                messagebox.showerror("Error Message", "File is empty")
            else:
                with open(filename, 'r') as f:
                    for line in f:
                        (key, val) = line.split()
                        reads[int(key)] = val

        elif not filename:
            pass
        elif filename.endswith('.txt'):
            filename = ''
            issequencesaved = NOFILE
            messagebox.showerror("Error Message", "Incorrect file extension")

    except Exception:
        filename = ''
        issequencesaved = NOFILE
        messagebox.showerror("Error Message", "Incorrect file format or file does not exist")

    return filename


def save_sequence():
    """
      Saves the assembled sequence to a designated location.
      
      Args:
          NONE

      Returns:
          NONE

      Raises:
          KeyError: Raises an exception.
     """
    global filename
    global issequencesaved

    if filename == '':
        messagebox.showerror("Error Message", "No import file")
    elif issequencesaved == False:
        f = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
        if f is None:
            return
        f.write(assembled_reads)
        f.close()
        issequencesaved = True
    else:
        messagebox.showerror("Error Message", "No assembled sequence")



def get_overlap(left, right):
    """
      Gets the amount of nucleotides by which the right read overlaps the left read.
      
      Args:
          left - the left read.
          right - the right read.
    
      Returns:
         len(left[i:]) -  The amount by which the right strand overlaps the left strand.
         sys.maxsize - Python's max int value. Indicates that the right strand is a useless or insignificant read.
         0 - Indicates no overlap.

      Raises:
          KeyError: Raises an exception.
      """
    right_strand = right
    for i in range(len(left)):
        if left[i:] == right[:len(left)-i]:
            return len(left[i:])
        elif left[i:i+len(right_strand)] == right_strand:
            return sys.maxsize
    return 0


def get_all_overlaps(strands, useless_reads):
    """
      Creates an overlap matrix (essentially a 2D dictionary) containing all the overlap combinations.
      
      Args:
          strands - The dictionary of strands created in read_data()
          useless_reads - Insignificant reads marked by sys.maxsize in get_overlap()
    
      Returns:
          overlap_matrix - a 2D dictionary containing all overlap combinations.

      Raises:
          KeyError: Raises an exception.
      """
    overlap_matrix = dict()
    for key1, sequence1 in strands.items():
        for key2, sequence2 in strands.items():
            if key1 == key2:
                continue
            if key1 not in overlap_matrix:
                overlap_matrix[key1] = dict()
            overlap_matrix[key1][key2] = get_overlap(sequence1, sequence2)
            if overlap_matrix[key1][key2] == sys.maxsize:
                useless_reads.append(key2)
    return overlap_matrix


def delete_useless_reads(overlap_matrix, useless_reads):
    """
      Deletes all reads marked insignificant by sys.maxsize in get_overlap()
      
      Args:
          overlap_matrix - a 2D dictionary containing all overlap combinations.
          useless_reads - reads marked insignificant by sys.maxsize in get_overlap()
    
      Returns:
          NONE

      Raises:
          KeyError: Raises an exception.
      """
    for junk in useless_reads:
        del overlap_matrix[junk]
        for key in overlap_matrix:
            del overlap_matrix[key][junk]


def find_first_read(overlap_matrix):
    """
      Finds the first read (i.e. the leftmost read) in the sequence. The read that results in the least amount of
      overlaps when placed on the right of other reads is considered the first read.
      
      Args:
          overlap_matrix - a 2D dictionary containing all overlap combinations.
    
      Returns:
          templist.index(min(templist)) + 1 - The index to the dictionary key containing the first read.
    
      Raises:
          KeyError: Raises an exception.
      """
    templist = dict()
    sum_overlaps = 0
    for i in overlap_matrix:
        for j in overlap_matrix[i]:
            sum_overlaps += overlap_matrix[j][i]
        templist[i] = sum_overlaps
    print("Templist = ", templist)
    return min(templist, key=templist.get)


def find_next_read(overlap_matrix):
    """
      Finds the next read by iterating through the overlap_matrix. The read resulting in the largest overlap of the
      current read is considered the next read.
      
      Args:
          overlap_matrix - a 2D dictionary containing all overlap combinations.
    
      Returns:
          next_read - The index to the dictionary key containing the next read.
    
      Raises:
          KeyError: Raises an exception.
      """
    largest = 0
    next_read = 0
    for i in overlap_matrix:
        if overlap_matrix[i] > largest:
            largest = overlap_matrix[i]
            next_read = i
    return next_read


def find_order(first_read, overlap_matrix):
    """
      Finds the order of reads for assembly.
      
      Args:
          first_read - the first (i.e. the leftmost read) in the sequence.
          overlap_matrix - a 2D dictionary containing all overlap combinations.

      Returns:
          order = A list containing the order of assembly.
    
      Raises:
          KeyError: Raises an exception.
      """
    order = []
    for i in range(0, len(overlap_matrix)):
        if i == 0:
            order.append(first_read)
        elif i == 1:
            next_read = find_next_read(overlap_matrix[first_read])
            order.append(next_read)
        else:
            next_read = find_next_read(overlap_matrix[next_read])
            order.append(next_read)
    return order


def assemble_reads(read_order):
    """
      Assembles the reads using the order of reads provided by find_order()
      
      Args:
          read_order - A list containing the order of assembly.
    
      Returns:
          assembled_reads - the assembled genome/sequence.
    
      Raises:
          KeyError: Raises an exception.
      """
    global assembled_reads
    for i in range(0, (len(read_order) - 1)):
        if i == 0:
            assembled_reads = reads[read_order[i]]
        lap = get_overlap(reads[read_order[i]], reads[read_order[i + 1]])
        temp = reads[read_order[i + 1]][lap:]
        assembled_reads += temp
    return assembled_reads


def run_assembly():
    """
      Runs the assembly by utilizing required functions.
      
      Args:
          NONE
    
      Returns:
          assembled_reads - the assembled genome/sequence.
    
      Raises:
          KeyError: Raises an exception.
      """

    if filename != '':
        useless_reads = []

        matrix_of_overlaps = (get_all_overlaps(reads, useless_reads))

        useless_reads = list(set(useless_reads))
        useless_reads.sort()

        delete_useless_reads(matrix_of_overlaps, useless_reads)

        print("X = ", useless_reads)

        print("MODIFIED OVERLAPS = ", matrix_of_overlaps)

        first = (find_first_read(matrix_of_overlaps))
        print("FIRST = ", first)

        order = find_order(first, matrix_of_overlaps)
        print("ORDER = ", order)

        global assembled_reads
        assembled_reads = "".join(assemble_reads(order))
        print("ASSEMBLED = ", assembled_reads)

        global issequencesaved
        issequencesaved = False

        # REFERENCE = atgagccaagttccgaacaaggattcgcggggaggatagatcagcgcccgagaggggtgagtcggtaaagagcattggaacgtcggagatacaactcccaagaaggaaaaaagagaaagcaagaagcggatgaatttccccataacgccagtgaaactctaggaaggggaaagagggaaggtggaagagaaggaggcgggcctcccgatccgaggggcccggcggccaagtttggaggacactccggcccgaagggttgagagtaccccagagggaggaagccacacggagtagaacagagaaatcacctccagaggaccccttcagcgaacagagagcgcatcgcgagagggagtagaccatagcgataggaggggatgctaggagttgggggagaccgaagcgaggaggaaagcaaagagagcagcggggctagcaggtgggtgttccgccccccgagaggggacgagtgaggcttatcccggggaactcgacttatcgtccccacatagcagactcccggaccccctttcaaagtgaccgaggggggtgactttgaacattggggaccagtggagccatgggatgctcctcccgattccgcccaagctccttccccccaagggtcgcccaggaatggcgggaccccactctgcagggtccgcgttccatcctttcttacctgatggccggcatggtcccagcctcctcgctggcgccggctgggcaacattccgaggggaccgtcccctcggtaatggcgaatgggacccacaaatctctctagcttcccagagagaagcgagagaaaagtggctctcccttagccatccgagtggacgtgcgtcctccttcggatgcccaggtcggaccgcgaggaggtggagatgccatgccgacccgaagaggaaagaaggacgcgagacgcaaacctgcgagtggaaacccgctttattcactggggtcgacaactctggggagaggagggagggtcggctgggaagagtatatcctatgggaatccctggcttccccttatgtccagtccctccccggtccgagtaaagggggactccgggactccttgcatgctggggacgaagccgcccccgggcgctcccctcgttccaccttcgagggggttcacacccccaacctgcgggccggctattcttctttcccttctctcgtcttcctcggtcaacctcctaagttcctcttcctcctccttgctgaggttctttccccccgccgatagctgctttctcttgttctcgagggccttccttcgtcggtgatcctgcctctccttgtcggtgaatcctcccctggaaggcctcttcctaggtccggagtctacttccatctggtccgttcgggccctcttcgccgggggagccccctctccatccttatctttctttccgagaattcctttgatgtttcccagccagggatgttcatcctcaagtttcttgattttcttcttaaccttccggaggtctctctcgagttcctctaacttctttcttccgctcacccactgctcgagaacctcttctctccccccgcggtttttccttccttcgggccggctcatcttcgactagaggcgacggtcctcagtactcttactcttttctgtaaagaggagactgctggccctgtcgcccaagttcgag
        # ASSEMBLED = atgagccaagttccgaacaaggattcgcggggaggatagatcagcgcccgagaggggtgagtcggtaaagagcattggaacgtcggagatacaactcccaagaaggaaaaaagagaaagcaagaagcggatgaatttccccataacgccagtgaaactctaggaaggggaaagagggaaggtggaagagaaggaggcgggcctcccgatccgaggggcccggcggccaagtttggaggacactccggcccgaagggttgagagtaccccagagggaggaagccacacggagtagaacagagaaatcacctccagaggaccccttcagcgaacagagagcgcatcgcgagagggagtagaccatagcgataggaggggatgctaggagttgggggagaccgaagcgaggaggaaagcaaagagagcagcggggctagcaggtgggtgttccgccccccgagaggggacgagtgaggcttatcccggggaactcgacttatcgtccccacatagcagactcccggaccccctttcaaagtgaccgaggggggtgactttgaacattggggaccagtggagccatgggatgctcctcccgattccgcccaagctccttccccccaagggtcgcccaggaatggcgggaccccactctgcagggtccgcgttccatcctttcttacctgatggccggcatggtcccagcctcctcgctggcgccggctgggcaacattccgaggggaccgtcccctcggtaatggcgaatgggacccacaaatctctctagcttcccagagagaagcgagagaaaagtggctctcccttagccatccgagtggacgtgcgtcctccttcggatgcccaggtcggaccgcgaggaggtggagatgccatgccgacccgaagaggaaagaaggacgcgagacgcaaacctgcgagtggaaacccgctttattcactggggtcgacaactctggggagaggagggagggtcggctgggaagagtatatcctatgggaatccctggcttccccttatgtccagtccctccccggtccgagtaaagggggactccgggactccttgcatgctggggacgaagccgcccccgggcgctcccctcgttccaccttcgagggggttcacacccccaacctgcgggccggctattcttctttcccttctctcgtcttcctcggtcaacctcctaagttcctcttcctcctccttgctgaggttctttccccccgccgatagctgctttctcttgttctcgagggccttccttcgtcggtgatcctgcctctccttgtcggtgaatcctcccctggaaggcctcttcctaggtccggagtctacttccatctggtccgttcgggccctcttcgccgggggagccccctctccatccttatctttctttccgagaattcctttgatgtttcccagccagggatgttcatcctcaagtttcttgattttcttcttaaccttccggaggtctctctcgagttcctctaacttctttcttccgctcacccactgctcgagaacctcttctctccccccgcggtttttccttccttcgggccggctcatcttcgactagaggcgacggtcctcagtactcttactcttttctgtaaagaggagactgctggccctgtcgcccaagttcgag

    else:
        messagebox.showerror("Error Message", "No input file or incorrect file format")


def instructions_page():
    """
      Outputs instructions page upon 'Instructions' button being pressed.
      
      Args:
          NONE
    
      Returns:
          NONE
    
      Raises:
          KeyError: Raises an exception.
      """
    instructions = Toplevel()

    instructions.title("Instructions")
    # http://www.python-course.eu/tkinter_text_widget.php <-- to help with text
    instructions.title = Message(instructions, text="Welcome to the DNA sequence assembler by Kevin Jian. "
                                              "This program takes shogunned DNA strands in the form of "
                                              "the nucleotide bases Cytosine(C), Guanine(G), Adenine(A), "
                                              "and Thymine(T), and assembles them by looking for regions "
                                              "of overlap. When you are ready, Import your file by "
                                              "clicking file:Import on the top left of the main screen. "
                                              "\nThe file must be a .txt and in the format of: sequence #, "
                                              "(space), sequence, (newline).\n\n"
                                              "Example of correct format:\n"
                                              "1 ATCG\n2 ATGG\n3 AGGG\n"
                                              "\nWhen everything is set, press 'Assemble'. If you wish to display "
                                              "The sequence, press 'display' and you will be prompted to "
                                              "input the amount of strands you wish to display, the max "
                                              "number of strands is 5,000. After the sequence is assembled, "
                                              "export the assembled sequence by clicking file:export.")

    # Welcome to the DNA sequence assembler by Kevin Jian. This software takes shotgunned DNA strands in the form of the nucleotide bases Cytosine(C), Guanine(G), Adenine(A), and Thymine(T), and assembles them by looking for regions of overlap. When you are ready, Import your file by clicking file:Import on the top left of the main screen. Note: the format of the file is very specific. Each individual DNA sequence must be on its own line. When everything is set, press 'Assemble'.
    instructions.title.pack()


def display_pg():
    """
      Outputs assembly page upon the 'Assemble' button being pressed. The assembly page contains an input field and a
      display button. The user, if they desire, has the option to input an x amount of nucleotides to be displayed,
      with the maximum being 5000. Once the user presses the display button, x nucleotides is displayed on screen.
      
      Args:
          NONE
    
      Returns:
          NONE
    
      Raises:
          KeyError: Raises an exception.
      """

    if filename != '' and issequencesaved != NOFILE:
        # run_assembly()
        #
        # global issequencesaved
        # issequencesaved = FILENOTSAVED

        assemblypage = Toplevel()
        assemblypage.resizable(0, 0)

        assemblypage.L1 = Label(assemblypage, text=("Nucleotides to be displayed (max: {})"
                                                    .format(len(assembled_reads))))
        assemblypage.L1.pack(side=TOP)

        assemblypage.E1 = Entry(assemblypage)
        assemblypage.E1.pack(side=LEFT)

    else:
        messagebox.showerror("Error Message", "No sequence to display")

    def display():
        strands = int(assemblypage.E1.get())

        if strands > MAX_NUCLEOTIDES:
            strands = MAX_NUCLEOTIDES

        if 0 < strands <= len(assembled_reads):
            displaypage = Toplevel()
            displaypage.resizable(0, 0)

            scrollbar = Scrollbar(displaypage)
            sequence = Text(displaypage, font="TkDefaultFont", height=10, width=50)

            scrollbar.pack(side=RIGHT, fill=Y)
            sequence.pack(side=LEFT, fill=Y)

            scrollbar.config(command=sequence.yview)
            sequence.config(yscrollcommand=scrollbar.set)

            sequence.insert(END, assembled_reads[0:strands])
            sequence.config(state=DISABLED)
            sequence.bind("<1>", lambda event: sequence.focus_set())

        elif strands > len(assembled_reads):
            messagebox.showerror("Error Message", "Input is greater than nucleotide count")

        else:
            messagebox.showerror("Error Message", "Input is negative")

    assemblypage.B1 = Button(assemblypage, text="Display", command=display)
    assemblypage.B1.pack(side=BOTTOM)


class NGSHomeScreen:
    def __init__(self, master):
        """
            Handles all functionality for the GUI.

            Args:
                NONE
          
            Returns:
                NONE

            Raises:
                KeyError: Raises an exception.
            """

        master.geometry('405x200+200+200')
        menu = Menu(master)
        master.config(menu=menu)
        filemenu = Menu(master)
        importmenu = Menu(master)

        menu.add_cascade(label="File", menu=filemenu)
        filemenu.add_cascade(label="Import", menu=importmenu)
        filemenu.add_cascade(label="Save Sequence As...", command=save_sequence)
        importmenu.add_command(label="DNA Strands",  command=lambda: read_and_display_filename())

        def read_and_display_filename():
            """
              Reads the user inputted file and displays the filename on GUI.
              
              Args:
                  NONE
            
              Returns:
                  This is a description of what is returned.
            
              Raises:
                  KeyError: Raises an exception.
              """
            read_data()
            global filename
            self.infile_name_display.configure(text="Input File: %s" % filename)

        self.title = Label(master, text="DNA Sequence Assembler\n")
        self.title.pack()

        self.infile_name_display = Label(master, text="Input File: ")
        self.infile_name_display.pack()

        # def check_saved_file():
        #     if issequencesaved == False:
        #         messagebox.showerror("Warning", "Assembled sequence was not saved")
        #     elif messagebox.askokcancel("Warning", "Do you want to exit?"):
        #         quit()

        # self.exit_button = Button(master, text="Exit", command=lambda: check_saved_file())
        # self.exit_button.pack(side=BOTTOM)

        self.about_button = Button(master, text="Instructions", command=instructions_page)
        self.about_button.pack(side=BOTTOM)

        self.display_button = Button(master, text="Display", command=lambda: display_pg())
        self.display_button.pack(side=BOTTOM)

        self.assemble_button = Button(master, text="Assemble", command=lambda: run_assembly())
        self.assemble_button.pack(side=BOTTOM)

NGSHomeScreen(root)


def on_closing():
    if issequencesaved == False:
        messagebox.showerror("Warning", "Assembled sequence was not saved")
    elif messagebox.askokcancel("Warning", "Do you want to exit?"):
        quit()

root.protocol("WM_DELETE_WINDOW", on_closing)
root.mainloop()
