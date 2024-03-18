#!/usr/bin/env python
# Experimental user interface for peptide library generation
#

import tkinter as tk
import tkinter.filedialog
import tkinter.messagebox
import os
import tempfile
import locale
import operator
import pdb
from subprocess import call

locale.setlocale(locale.LC_ALL, "")

# logging
import sys

# sys.stderr = open('my_stderr.txt', 'w')
# sys.stdout = open('my_stdout.txt', 'w')

import Image as PIL
import ImageTk as piltk

# import pybel
# import openbabel
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import PepLibGen

# Check which PepLibGen is used - useful for building with py2app/py2exe
# print "PepLibGen", PepLibGen.__file__
from PepLibGen.StructGen import zinc_peptides
from PepLibGen.StructGen import synthrules
from PepLibGen.StructGen import StructGen as sg
from PepLibGen.StructGen import aminoacids as aa
from PepLibGen.Analysis import chem_analysis as ca

# Hack to make RDKIT play nice when packaged...
if "RDBASE" not in list(os.environ.keys()):
    os.environ["RDBASE"] = os.getcwd()

zinc_def_file = "zinc_output"
if not os.path.exists(zinc_def_file):
    zinc_def_file = os.path.join(os.path.dirname(sys.argv[0]), zinc_def_file)

# Warn if trying to generate more than this number of peptides
MAX_PEPTIDES = 1000000


class PepApp(tk.Frame):
    """
    Main class of pepgui.py, where the main application window is defined
    """

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.grid(sticky="NSEW")
        self.single_frame = SinglePepFrame()
        self.single_frame.grid(row=1, column=0)
        self.lib_frame = PepLibraryFrame()
        self.lib_frame.grid(
            row=1,
            column=0,
            columnspan=2,
            sticky=tk.N + tk.S + tk.E + tk.W,
            padx=10,
            pady=10,
        )
        self.lib_frame.grid_remove()
        self.chooseFrame()
        self.master.title("CycloPs")
        self.temp_holder, self.temp_filename = tempfile.mkstemp()
        # Make resizeable
        top = self.winfo_toplevel()
        top.resizable(0, 0)
        top.rowconfigure(0, weight=1)
        top.columnconfigure(0, weight=1)
        # self.rowconfigure(0, weight =1)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.columnconfigure(1, weight=1)
        self.topLevelMenu()
        # Variables
        self.on_off = dict(
            [(amino, tk.IntVar()) for amino in list(sg.all_aminos.keys())]
        )

    def topLevelMenu(self):
        """
        Create top level menu

        """
        top = self.winfo_toplevel()
        self.topMenu = tk.Menu(top)
        top["menu"] = self.topMenu
        self.optionsMenu = tk.Menu(self.topMenu)
        self.topMenu.add_cascade(
            label="Modify Amino-acid library", menu=self.optionsMenu
        )
        self.optionsMenu.add_command(label="Exit", command=self.quit)

        # Add option to edit available amino acids.
        self.optionsMenu.add_command(
            label="Choose L-Amino Acids", command=self.naturalAminoChooser
        )
        self.optionsMenu.add_command(
            label="Choose D-Amino Acids", command=self.dAminoChooser
        )
        self.optionsMenu.add_command(
            label="Choose special Amino Acids", command=self.specialAminoChooser
        )
        self.optionsMenu.add_command(label="ZINC amino acids", command=self.zinc_aminos)

    def destroy(self):
        tk.Frame.destroy(self)
        try:
            os.close(self.temp_handle)
            os.remove(self.temp_filename)
        except:
            pass

    def chooseFrame(self):
        """
        Buttons to choose between single peptide frame or library frame.

        Will not show both frames at once - hide one, and show one.
        """
        self.show_pep_frame = tk.Button(
            self, text="Peptide Generator", command=self.displayPepFrame, height=1
        )
        self.show_pep_frame.grid(row=0, column=0, sticky=tk.E + tk.N + tk.S)

        self.show_lib_frame = tk.Button(
            self,
            text="Peptide Library Generator",
            command=self.displayLibFrame,
            height=1,
        )

        self.show_lib_frame.grid(row=0, column=1, sticky=tk.W + tk.N + tk.S)

        self.show_pep_frame.config(relief=tk.SUNKEN)

    def displayPepFrame(self):
        """
        Hides the library generation frame, shows the peptide generation frame
        """
        self.lib_frame.grid_remove()
        self.single_frame.grid()
        self.show_pep_frame.config(relief=tk.SUNKEN)
        self.show_lib_frame.config(relief=tk.RAISED)

    def displayLibFrame(self):
        """
        Hides the peptide generation frame, shows the library generation frame
        """
        self.single_frame.grid_remove()
        self.lib_frame.grid()
        self.show_lib_frame.config(relief=tk.SUNKEN)
        self.show_pep_frame.config(relief=tk.RAISED)

    def aminoChooser(self, amino_type):
        """
        Presents window for selection of amino acids used to generate peptides
        """
        self.chooser = tk.Toplevel(padx=50)
        self.chooser.title("Amino Chooser")
        self.chooser.resizable(1, 0)

        self.chooser_frame = tk.LabelFrame(self.chooser)
        if amino_type == aa.aminos:
            self.chooser_frame.config(text="Natural Amino Acids")
        elif amino_type == aa.d_aminos:
            self.chooser_frame.config(text="D-Amino Acids")
        elif amino_type == aa.special_aminos:
            self.chooser_frame.config(text="Special Amino Acids")

        self.chooser_frame.grid(sticky=tk.E + tk.W)
        self.checkbuttons = dict([(amino, None) for amino in sg.all_aminos])

        for amino in sorted(sg.all_aminos.keys()):
            if amino in list(sg.aminodata.keys()):
                # print '%s is in aminodata' % (amino)
                self.on_off[amino].set(1)
            self.checkbuttons[amino] = tk.Checkbutton(
                self.chooser_frame,
                text=amino,
                variable=self.on_off[amino],
                command=self.add_remove_amino,
            )
            if amino in amino_type:
                self.checkbuttons[amino].grid(sticky=tk.W)

        all_toggle = lambda: self.toggleAminoSet(amino_type)
        self.choose_all = tk.Button(
            self.chooser_frame,
            command=all_toggle,
            text="Choose/Unchoose all.",
            height=1,
        )
        self.choose_all.grid()

    def toggleAminoSet(self, amino_set):
        # pdb.set_trace()
        state = [self.on_off[a].get() for a in list(amino_set.keys())]
        if not len(state) == len([s for s in state if s]):
            # If not all are on, turn them on
            [self.on_off[name].set(1) for name in amino_set]
        else:
            [self.on_off[name].set(0) for name in amino_set]
        self.add_remove_amino()

    def naturalAminoChooser(self):
        """
        Presents window to add/remove natural amino acids from working library
        """
        self.aminoChooser(aa.aminos)
        # print len(aa.aminos)

    def dAminoChooser(self):
        """
        Presents window to add/remove d amino acids from working library
        """
        self.aminoChooser(aa.d_aminos)
        # print len(aa.d_aminos)

    def specialAminoChooser(self):
        """
        Presents window to add/remove d amino acids from working library
        """
        self.aminoChooser(aa.special_aminos)
        # print len(aa.special_aminos)

    def add_remove_amino(self):
        """
        Adds/removes amino acid from library on checkbox tick.
        Syncs state of checkboxes with amino-acid library.
        """
        for name in list(sg.all_aminos.keys()):
            if self.on_off[name].get() and name not in list(sg.aminodata.keys()):
                # print name
                sg.add_amino(name)
            elif not self.on_off[name].get() and name in list(sg.aminodata.keys()):
                sg.remove_amino(name)

        # print len(sg.aminodata)

    def zinc_aminos(self):
        test = ZincChooser(self)


class SinglePepFrame(tk.Frame):
    """
    Holds the interface to work with a single peptide
    """

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.grid(sticky="NSEW")
        self.createSinglePepWidgets()
        self.temp_holder, self.temp_filename = tempfile.mkstemp()

    def destroy(self):
        tk.Frame.destroy(self)
        try:
            os.close(self.temp_handle)
            os.remove(self.temp_filename)
        except:
            pass

    def createSinglePepWidgets(self):
        """
        Create widgets inside a frame that deal with generating a single peptide
        """

        # Variables
        self.seq_holder = tk.StringVar()  # Holds peptide sequence
        self.old_sequence = ""  # Used to set constraint to linear
        self.num_resis = tk.IntVar()  # Value used for list selection
        self.num_resis.set(5)
        # when new sequence is entered
        self.cur_constraint = tk.StringVar()  # Selected constraint
        self.num_menu_items = 0
        self.labels = []  # Keep references to menu labels for easy deletion
        self.smiles_holder = tk.StringVar()
        self.imagedata = ""  # Reference to image to hold it around...
        # Can be synthesised...
        self.can_synthesise = tk.StringVar()
        # Defaults to menus
        self.text_entry = tk.IntVar()
        self.text_entry.set(0)
        self.logp = tk.StringVar()
        self.druglike = tk.StringVar()

        # Sequence entry form
        self.seq_frame = tk.LabelFrame(self, text="Enter Sequence Here")
        self.seq_frame.grid(row=0, columnspan=3, sticky=tk.E + tk.W)
        self.seq_frame.columnconfigure(0, weight=1)
        self.SequenceEntryWidgetInitialise()
        self.SequenceEntryUpdate()
        # Peptide builder
        self.pep_build_button = tk.Button(
            self.seq_frame, text="Peptide Builder", command=self.CallPeptideBuilder
        )
        self.pep_build_button.grid()

        # Constraint menu
        self.cur_constraint.set("linear")
        self.constButton = tk.Menubutton(
            self.seq_frame, text="Select desired constraint", relief="raised"
        )
        self.constButton.grid(row=0, column=self.num_resis.get())
        self.constButton.menu = tk.Menu(self.constButton, tearoff=0)
        self.constButton["menu"] = self.constButton.menu
        self.constButton.menu.add_radiobutton(
            label="linear",
            value="linear",
            variable=self.cur_constraint,
            command=self.SmilesHandlerWrapper,
        )
        self.constButton.bind("<1>", self.PopulateConstraintMenu)

        # Label to display SMILES
        self.smiles_frame = tk.LabelFrame(self, text="Sequence SMILES")
        self.smiles_frame.grid(
            row=1, column=0, columnspan=3, sticky=tk.N + tk.E + tk.S + tk.W
        )
        self.smiles_frame.columnconfigure(0, weight=1)
        self.smiles_label = tk.Entry(
            self.smiles_frame,
            state="readonly",
            readonlybackground="white",
            textvariable=self.smiles_holder,
        )
        self.smiles_label.grid(sticky=tk.W + tk.E)

        # Label to display around synthesis rules
        self.synth_frame = tk.LabelFrame(self, text="Synthesisability")
        self.synth_frame.grid(
            row=2, column=0, columnspan=3, sticky=tk.N + tk.E + tk.S + tk.W
        )
        self.synth_frame.columnconfigure(0, weight=1)
        self.synth_label = tk.Label(
            self.synth_frame, textvariable=self.can_synthesise, justify="left"
        )
        self.synth_label.grid(sticky=tk.W + tk.E)

        # Label to display around partition coefficient
        self.logp_frame = tk.LabelFrame(
            self, text="Estimated Octanol-Water Partion Coefficient (logP)"
        )
        self.logp_frame.grid(row=3, column=0, columnspan=3, sticky=tk.W + tk.E)
        self.logp_frame.columnconfigure(0, weight=1)
        self.logp_label = tk.Entry(
            self.logp_frame,
            textvariable=self.logp,
            justify=tk.LEFT,
            state="readonly",
            readonlybackground="white",
            relief="flat",
        )
        self.logp_label.grid(sticky=tk.E + tk.W)

        # Label to display around Druglike criteria
        self.druglike_frame = tk.LabelFrame(self, text="Druglike Properties")
        self.druglike_frame.grid(row=4, column=0, columnspan=3, sticky=tk.E + tk.W)
        self.druglike_frame.columnconfigure(0, weight=1)

        # Label to display around Molecule Image
        self.image_frame = tk.LabelFrame(self, text="Structure")
        self.image_frame.grid(
            row=0, column=4, columnspan=3, rowspan=5, sticky=tk.N + tk.E + tk.S + tk.W
        )
        self.image_frame.columnconfigure(0, weight=1)

    def CallPeptideBuilder(self):
        """
        Starts the peptide builder
        """

        def close_handler():
            try:
                self.SmilesHandlerWrapper()
            finally:
                self.builder.destroy()

        self.builder = PeptideBuilder(self, seq_holder=self.seq_holder)
        self.builder.protocol("WM_DELETE_WINDOW", close_handler)

    def SequenceEntryWidgetInitialise(self):
        """
        Draws the sequence entry widgets

        Both entry box, and menu lists

        """
        # Text entry box
        self.sequenceEntry = tk.Entry(self.seq_frame, textvariable=self.seq_holder)

    def SequenceEntryUpdate(self):
        """
        Chooses which entry widget is displayed with self.text_entry
        """

        self.sequenceEntry.grid(row=0, sticky=tk.N + tk.S + tk.E + tk.W)
        self.sequenceEntry.bind("<KeyPress-Return>", self.SmilesHandler)

    def UpdateSequenceEntryBoxes(self, event):
        """
        Updates the list of aminos in the sequence entry boxes.
        """
        for i in self.sequence_menu:
            self.sequence_menu[i].menu.delete(0, tk.END)
            for amino in sg.aminodata:
                self.sequence_menu[i].menu.add_radiobutton(
                    label=amino,
                    value=amino,
                    variable=self.sequence_menu_variables[i],
                    command=self.UpdateSeq,
                )

    def UpdateSeq(self):
        """
        Sets the sequence StringVar from the entry menu lists
        """
        new_seq = ""
        for i in self.sequence_menu_variables:
            if self.sequence_menu_variables[i].get() != "Select Amino Acid":
                new_seq += "," + self.sequence_menu_variables[i].get()
        self.seq_holder.set(new_seq.strip(","))
        # print self.seq_holder.get()
        self.SmilesHandler("")

    def PopulateConstraintMenu(self, event):
        """
        Sets the options in the drop down list to be the available constraint
        options.
        """
        seq = self.seq_holder.get()
        seq = seq.strip("\n")
        if "," in seq or "-" in seq or "ZINC" in seq:
            seq = seq.split(",")

        # print seq
        all_constraints = sg.what_constraints(seq)
        # Make return options more readable in menu
        for lab in self.labels:
            self.constButton.menu.delete(lab)

        self.labels = []
        for con in all_constraints:
            if "SS" in con:
                lab = "Disulphide Bond " + con
                self.constButton.menu.add_radiobutton(
                    label=lab,
                    value=con,
                    variable=self.cur_constraint,
                    command=self.SmilesHandlerWrapper,
                )
                self.labels.append(lab)
            if "HT" in con:
                lab = "Head-Tail Bond"
                self.constButton.menu.add_radiobutton(
                    label=lab,
                    value=con,
                    variable=self.cur_constraint,
                    command=self.SmilesHandlerWrapper,
                )
                self.labels.append(lab)
            if "SC" in con:
                lab = "Side-Chain Bond " + con
                self.constButton.menu.add_radiobutton(
                    label=lab,
                    value=con,
                    variable=self.cur_constraint,
                    command=self.SmilesHandlerWrapper,
                )
                self.labels.append(lab)

    def SmilesHandlerWrapper(self):
        """
        Calls self SmilesHandler when constraint is chosen
        """
        self.SmilesHandler("")

    def SmilesHandler(self, event=""):
        """
        Sets the appropriate label text to the SMILES of the input peptide
        sequence.

        Also displays the molecule image and synthesis rules
        """
        # If this is the first time the SMILES has been written - set the
        # constraint to a default of 'linear' as it is the only one that
        # must really exist
        if self.seq_holder.get() != self.old_sequence:
            self.cur_constraint.set("linear")
        if (
            "," in self.seq_holder.get()
            or "-" in self.seq_holder.get()
            or "ZINC" in self.seq_holder.get()
        ):
            seq = self.seq_holder.get().replace("\n", "").split(",")
        else:
            seq = self.seq_holder.get()

        # print seq
        # print 'Constraint', self.cur_constraint.get()
        try:
            if self.cur_constraint.get() == "linear":
                # print 'Smiles',
                self.smiles_holder.set(sg.linear_peptide_smiles(seq))
            else:
                constraint = self.cur_constraint.get()
                # print 'Smiles, const', seq, constraint
                self.smiles_holder.set(
                    sg.constrained_peptide_smiles(seq, constraint)[2]
                )
            self.old_sequence = self.seq_holder.get()
            self.DrawMoleculeImage()
            self.update_synthesis_rules()
            self.update_logp()
            self.update_druglike()
            # print self.smiles_holder.get()
            mol = Chem.MolFromSmiles(self.smiles_holder.get())
            self.smiles_holder.set(Chem.MolToSmiles(mol, True))
        except sg.UndefinedAminoError:

            # print 'SMILES GETTER ERROR'
            self.smiles_holder.set("None")

    def DrawMoleculeImage(self):
        """
        Adds label displaying a 2D model of the (constrained) peptide
        """
        try:
            self.image_label.destroy()
        except AttributeError:  # If it doesn't already exist
            pass
        try:
            # self.mol_obj = pybel.readstring('smi',self.smiles_holder.get())
            # self.mol_obj.draw(show=False, filename=self.temp_filename)
            # image = PIL.open(self.temp_filename)
            self.mol_obj = Chem.MolFromSmiles(self.smiles_holder.get())
            image = Draw.MolToImage(self.mol_obj, size=(350, 350))
            self.largeimage = Draw.MolToImage(self.mol_obj, size=(1000, 1000))
            self.image = image  # Save a reference, so it can be saved
            self.imagedata = piltk.PhotoImage(image)
            self.image_label = tk.Label(self.image_frame, image=self.imagedata)
            self.image_label.grid(row=0, columnspan=2)
            # Add savebutton
            self.image_save = tk.Button(
                self.image_frame, text="Save Molecule", command=self.SaveMoleculeDrawing
            )
            self.image_save.grid(row=1, column=0, sticky=tk.E + tk.W)
            # Add 3D structure savebutton
            self.struct_save = tk.Button(
                self.image_frame, text="Write 3D structure", command=self.Write3D
            )
            self.struct_save.grid(row=1, column=1, columnspan=2)
        except IOError as e:  # Catch empty smiles strings
            pass

    def SaveMoleculeDrawing(self):
        """
        Saves molecule drawing to a file
        """
        try:
            self.largeimage.save(
                tkinter.filedialog.asksaveasfilename(
                    defaultextension=".png",
                    initialfile=self.seq_holder.get().strip()
                    + "-"
                    + self.cur_constraint.get(),
                    initialdir=os.path.expanduser("~"),
                )
            )
        except KeyError:
            self.largeimage.save(
                tkinter.filedialog.asksaveasfilename(
                    initialfile=self.seq_holder.get().strip()
                    + "-"
                    + self.cur_constraint.get()
                    + ".png",
                    initialdir=os.path.expanduser("~"),
                )
            )

    def Write3D(self):
        """
        Writes 3D structure to a file
        """
        # self.mol_obj.make3D()
        AllChem.EmbedMolecule(self.mol_obj)
        AllChem.UFFOptimizeMolecule(self.mol_obj)
        outfile = tkinter.filedialog.asksaveasfilename(
            defaultextension=".sdf",
            initialfile=self.seq_holder.get().strip() + "-" + self.cur_constraint.get(),
        )
        out = Chem.MolToMolBlock(self.mol_obj)
        try:
            handle = open(outfile, "w")
            handle.write(out)
        finally:
            handle.close()

    def update_synthesis_rules(self):
        """
        Display whether a peptide is synthesisable or not, and why
        """
        pep = self.seq_holder.get()
        self.synth_label.config(bg="white")
        self.synth_label.config(justify="left")
        smiles = self.smiles_holder.get()
        # pdb.set_trace()

        if self.cur_constraint.get() != "linear":
            bond_def = self.cur_constraint.get()[2:]
        else:
            bond_def = ""

        if "," in pep or "-" in pep or "ZINC" in pep:
            pep = pep.strip("\n").split(",")
            # print pep
            new_pep = "".join([sg.aminodata[resi]["Letter"] for resi in pep])
            if len(new_pep) == len(pep):
                pep = new_pep
                try:
                    test_res = synthrules.run_all_tests(smiles, pep, bond_def)
                    results = "\n".join(test_res)
                    self.synth_label.config(fg="red")
                    self.can_synthesise.set(results)
                except TypeError as e:
                    self.synth_label.config(fg="blue")
                    self.can_synthesise.set("PASSED")
            else:
                self.synth_label.config(fg="green")
                self.can_synthesise.set("Unknown")
        else:
            try:
                results = "\n".join(synthrules.run_all_tests(smiles, pep, bond_def))
                self.synth_label.config(fg="red")
                self.can_synthesise.set(results)
            except TypeError:
                self.synth_label.config(fg="blue")
                self.can_synthesise.set("PASSED")

    def update_logp(self):
        """
        Display molecule LogP value
        """
        self.logp.set(ca.log_partition_coefficient(self.smiles_holder.get()))

    def update_druglike(self):
        """
        Show which druglike rules the molecule passes and fails
        """
        try:
            for label in self.druglike_labels:
                label.destroy()
        except AttributeError:  # First run, doesn't exist
            pass
        passed, failed = ca.lipinski_trial(self.smiles_holder.get())

        self.passed = []
        self.failed = []
        self.druglike_labels = []
        for result in passed:
            self.passed.append(tk.StringVar())
            self.passed[-1].set("PASSED: " + result)
            lab = tk.Entry(
                self.druglike_frame,
                textvariable=self.passed[-1],
                fg="blue",
                state="readonly",
                readonlybackground="white",
                relief="flat",
            )
            lab.grid(sticky=tk.E + tk.W)
            self.druglike_labels.append(lab)
        for result in failed:
            self.failed.append(tk.StringVar())
            self.failed[-1].set("FAILED: " + result)
            lab = tk.Entry(
                self.druglike_frame,
                textvariable=self.failed[-1],
                fg="red",
                state="readonly",
                readonlybackground="white",
                relief="flat",
            )
            lab.grid(sticky=tk.E + tk.W)
            self.druglike_labels.append(lab)


class PepLibraryFrame(tk.Frame):
    """
    Holds the interface to work with a peptide library
    """

    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.createPepLibraryWidgets()

    def createPepLibraryWidgets(self):
        """
        Create widgets inside a frame that generate a library of peptides
        """
        # Main frame
        for row in range(1, 8):
            self.rowconfigure(row, weight=1)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

        # Variables
        self.seqs_from_file = tk.StringVar()
        self.lib_seq_holder = tk.StringVar()
        self.lib_len = tk.IntVar()
        self.gen_linears = tk.IntVar()
        self.num_lin = tk.IntVar()
        # self.num_lin.set('Linear peptides...                               ')
        self.gen_SS = tk.IntVar()
        self.num_SS = tk.IntVar()
        # self.num_SS.set('Disulphide bonded peptides...                     ')
        self.gen_HT = tk.IntVar()
        self.num_HT = tk.IntVar()
        # self.num_HT.set('Head-Tail bonded peptides...                       ')
        self.gen_SCSC = tk.IntVar()
        self.num_SCSC = tk.IntVar()
        # self.num_SCSC.set('Side-chain to side-chain bonded peptides...       ')
        self.gen_SCNT = tk.IntVar()
        self.num_SCNT = tk.IntVar()
        # self.num_SCNT.set('Side-chain to N-terminal bonded peptides...       ')
        self.gen_SCCT = tk.IntVar()
        self.num_SCCT = tk.IntVar()
        # self.num_SCCT.set('Side-chain to C-terminal bonded peptides...       ')
        self.gen_linears.set(1)
        self.smiles_gen = []  # Holder for all library smiles
        # Include/remove non synthesisable peptides
        self.synth_filter = tk.IntVar()
        # Include/remove druglike peptides
        self.druglike_only = tk.IntVar()
        # N-methylate peptides
        self.nmethylate = tk.IntVar()
        # Unreadable peptides found in input file
        self.unreadable_peptides = tk.IntVar()

        # Sequence entry form
        self.lib_entry_frame = tk.LabelFrame(self, text="Enter library pattern")
        self.lib_entry_frame.grid(row=0, columnspan=2, sticky=tk.N + tk.S + tk.E + tk.W)
        self.libEntry = tk.Entry(
            self.lib_entry_frame, textvariable=self.lib_seq_holder, width=60
        )
        self.libEntry.grid(sticky=tk.N + tk.S + tk.E + tk.W, columnspan=2)
        self.libEntry.bind("<Return>", self.countPosPeptides)
        self.libOptionText = tk.Label(
            self.lib_entry_frame,
            text="OR - choose a file containing peptide definitions "
            "(This will override library pattern and constraint choices.)",
        )
        self.libOptionText.grid(sticky=tk.N + tk.S + tk.E + tk.W, columnspan=2)
        self.fileEntry = tk.Entry(
            self.lib_entry_frame, textvariable=self.seqs_from_file, width=60
        )
        self.fileEntry.grid(sticky=tk.N + tk.S + tk.E + tk.W, columnspan=2)
        self.fileEntryButton = tk.Button(
            self.lib_entry_frame, text="Choose file", command=self.choosePeptideFile
        )
        self.fileEntryButton.grid(column=5, row=2)
        # When return is pressed, count the peptides
        self.fileEntry.bind("<Return>", self.tryPeptideFile)
        # When Delete is pressed, if the widget is empty, reactive the
        # library pattern and constraint widgets
        self.fileEntry.bind("<Delete>", self.reactivateLibraryPanel)
        self.fileEntry.bind("<BackSpace>", self.reactivateLibraryPanel)

        # Choose constraints...
        self.lib_constraint_frame = tk.LabelFrame(
            self, text="Select which types of peptides will be included in the library"
        )
        self.lib_constraint_frame.grid(
            row=1, rowspan=6, sticky=tk.N + tk.S + tk.E + tk.W
        )

        self.lin_button = tk.Checkbutton(
            self.lib_constraint_frame,
            text="Linear",
            variable=self.gen_linears,
            command=self.constraintSelectHandler,
        )
        self.lin_button.grid(sticky=tk.W)

        self.SS_button = tk.Checkbutton(
            self.lib_constraint_frame,
            text="Disulphide bonded",
            variable=self.gen_SS,
            command=self.constraintSelectHandler,
        )
        self.SS_button.grid(sticky=tk.W)

        self.HT_button = tk.Checkbutton(
            self.lib_constraint_frame,
            text="Head-tail bonded",
            variable=self.gen_HT,
            command=self.constraintSelectHandler,
        )
        self.HT_button.grid(sticky=tk.W)

        self.SCSC_button = tk.Checkbutton(
            self.lib_constraint_frame,
            text="Side-chain to side-chain bonded",
            variable=self.gen_SCSC,
            command=self.constraintSelectHandler,
        )
        self.SCSC_button.grid(sticky=tk.W)

        self.SCNT_button = tk.Checkbutton(
            self.lib_constraint_frame,
            text="Side-chain to N-terminus bond",
            variable=self.gen_SCNT,
            command=self.constraintSelectHandler,
        )
        self.SCNT_button.grid(sticky=tk.W)

        self.SCCT_button = tk.Checkbutton(
            self.lib_constraint_frame,
            text="Side-chain to C-terminus bonded",
            variable=self.gen_SCCT,
            command=self.constraintSelectHandler,
        )
        self.SCCT_button.grid(sticky=tk.W)

        # Display number of possible sequences
        self.all_info_frame = tk.LabelFrame(
            self, text="Maximum number of linear combinations"
        )
        self.all_info_frame.grid(row=7, columnspan=2, sticky=tk.N + tk.S + tk.E + tk.W)
        self.lib_info_panel = tk.Entry(
            self.all_info_frame,
            state="readonly",
            readonlybackground="white",
            textvariable=self.lib_len,
        )
        self.lib_info_panel.grid(column=0, row=0, sticky=tk.N + tk.S + tk.E + tk.W)

        self.warning_label = tk.Label(
            self.all_info_frame,
            text="Warning! Very large library - may crash while generating",
            fg="red",
        )
        # Generate SMILES button (next to possible sequences)
        self.generate_lib = tk.Button(
            self.all_info_frame,
            text="Generate Library SMILES anyway...",
            command=self.genLibrarySMILES,
        )
        self.generate_lib.grid(row=0, column=1, sticky=tk.N + tk.S + tk.E + tk.W)
        self.generate_lib.grid_remove()

        # Display number of generated peptides
        self.gen_info_frame = tk.LabelFrame(
            self, text="Number of peptides generated for library"
        )
        self.gen_info_frame.grid(
            column=1, row=1, rowspan=6, sticky=tk.N + tk.S + tk.E + tk.W
        )

        self.lin_text_panel = tk.Label(
            self.gen_info_frame, text="Linear peptides... ", anchor=tk.W
        )
        self.lin_text_panel.grid(sticky=tk.W, row=0, column=0)
        self.lin_info_panel = tk.Label(self.gen_info_frame, textvariable=self.num_lin)
        self.lin_info_panel.grid(sticky=tk.W, row=0, column=1)

        self.SS_text_panel = tk.Label(
            self.gen_info_frame, text="Disulphide bonded peptides... ", anchor=tk.W
        )
        self.SS_text_panel.grid(sticky=tk.W, row=1, column=0)
        self.SS_info_panel = tk.Label(self.gen_info_frame, textvariable=self.num_SS)
        self.SS_info_panel.grid(sticky=tk.W, row=1, column=1)

        self.HT_text_panel = tk.Label(
            self.gen_info_frame, text="Head-tail bonded peptides... ", anchor=tk.W
        )
        self.HT_text_panel.grid(sticky=tk.W, row=2, column=0)
        self.HT_info_panel = tk.Label(self.gen_info_frame, textvariable=self.num_HT)
        self.HT_info_panel.grid(sticky=tk.W, row=2, column=1)

        self.SCSC_text_panel = tk.Label(
            self.gen_info_frame,
            text="Side-chain to Side-chain bonded peptides... ",
            anchor=tk.W,
        )
        self.SCSC_text_panel.grid(sticky=tk.W, row=3, column=0)
        self.SCSC_info_panel = tk.Label(self.gen_info_frame, textvariable=self.num_SCSC)
        self.SCSC_info_panel.grid(sticky=tk.W, row=3, column=1)

        self.SCNT_text_panel = tk.Label(
            self.gen_info_frame,
            text="Side-chain to N-Terminus bonded peptides... ",
            anchor=tk.W,
        )
        self.SCNT_text_panel.grid(sticky=tk.W, row=4, column=0)
        self.SCNT_info_panel = tk.Label(self.gen_info_frame, textvariable=self.num_SCNT)
        self.SCNT_info_panel.grid(sticky=tk.W, row=4, column=1)

        self.SCCT_text_panel = tk.Label(
            self.gen_info_frame,
            text="Side-chain to C-Terminus bonded peptides... ",
            anchor=tk.W,
        )
        self.SCCT_text_panel.grid(sticky=tk.W, row=5, column=0)
        self.SCCT_info_panel = tk.Label(self.gen_info_frame, textvariable=self.num_SCCT)
        self.SCCT_info_panel.grid(sticky=tk.W, row=5, column=1)

        # Checkbutton to filter by synthrules
        self.synth_check = tk.Checkbutton(
            self,
            text="Filter out peptides that break synthesis rules",
            variable=self.synth_filter,
            command=self.constraintSelectHandler,
        )
        self.synth_check.grid(row=8, columnspan=2, sticky=tk.W)

        # Checkbutton to filter for druglike compounds
        self.druglike_check = tk.Checkbutton(
            self,
            text="Filter out non-druglike peptides "
            "(Most peptides over 3 residues in length are not druglike)",
            variable=self.druglike_only,
            command=self.constraintSelectHandler,
        )
        self.druglike_check.grid(row=9, columnspan=2, sticky=tk.W)

        # Checkbutton for N-methylating peptides
        self.nmethylation = tk.Checkbutton(
            self,
            text="N-methylate peptide backbone and n-terminus",
            variable=self.nmethylate,
            command=lambda: 0,
        )
        self.nmethylation.grid(row=10, columnspan=2, sticky=tk.W)

        # Write Indicator
        self.writelabel = tk.Label(self)
        self.writelabel.grid(row=14, sticky=tk.E + tk.W, columnspan=2)

        # Write SMILES button
        self.write_lib = tk.Button(
            self, text="Write SMILES to file", command=self.writeLibrarySMILES
        )
        self.write_lib.grid(columnspan=2, row=11, sticky=tk.N + tk.S + tk.E + tk.W)

        # Write 2D button
        self.write_png = tk.Button(
            self,
            text="Write 2D drawings of library to folder",
            command=self.drawLibrary,
        )
        self.write_png.grid(columnspan=2, row=12, sticky=tk.N + tk.S + tk.E + tk.W)

        # Write 3D buttons
        self.write_pdbs = tk.Button(
            self,
            text="Write 3D structures of library to folder",
            command=self.writeLibrary3D,
        )
        self.write_pdbs.grid(row=13, column=0, sticky=tk.N + tk.S + tk.E + tk.W)

        self.write_pdb = tk.Button(
            self,
            text="Write 3D structures of library to single .sdf file",
            command=self.writeLibrary3Dtofile,
        )
        self.write_pdb.grid(row=13, column=1, sticky=tk.N + tk.S + tk.E + tk.W)

    def tryPeptideFile(self, event):
        """Attempt to process the input peptide file, deactivate pattern
        and constraint widgets."""
        self.countPosPeptides("")
        self.libEntry.configure(state=tk.DISABLED)
        for checkbutton in [
            self.lin_button,
            self.SS_button,
            self.HT_button,
            self.SCSC_button,
            self.SCNT_button,
            self.SCCT_button,
        ]:
            checkbutton.configure(state=tk.DISABLED)

    def choosePeptideFile(self):
        """Set the file to generate peptides from and try to process it."""
        fhandle = tkinter.filedialog.askopenfile(
            title="Choose a peptide definition file"
        )
        fhandle.close()
        self.seqs_from_file.set(fhandle.name)
        self.tryPeptideFile("")

    def reactivateLibraryPanel(self, event):
        """If fileEntry label is blank, reactivate library pattern box."""
        if not self.seqs_from_file.get():
            self.libEntry.configure(state=tk.NORMAL)
            for checkbutton in [
                self.lin_button,
                self.SS_button,
                self.HT_button,
                self.SCSC_button,
                self.SCNT_button,
                self.SCCT_button,
            ]:
                checkbutton.configure(state=tk.NORMAL)

    def constraintSelectHandler(self):
        """
        Called when a constraint button is clicked to update peptide counts.

        Passes '' as the event parameter
        """
        self.countPosPeptides("")

    def countPosPeptides(self, event):
        """
        Displays the number of peptides to be generated in the GUI
        """
        self.generate_lib.grid_remove()
        if self.seqs_from_file.get():
            self.lib_len.set(self.countFilePeptides())
        else:
            self.lib_len.set(self.countPatternPeptides())
        self.warning_label.grid_forget()

        if self.lib_len.get() > MAX_PEPTIDES:
            self.warning_label.grid(sticky=tk.N + tk.S + tk.E + tk.W, column=0, row=1)
            self.generate_lib.grid()
            for var in [
                self.num_lin,
                self.num_SS,
                self.num_HT,
                self.num_SCSC,
                self.num_SCCT,
                self.num_SCNT,
            ]:
                var.set(0)
        else:
            self.genLibrarySMILES()

    def countPatternPeptides(self):
        """Counts peptides to be generated from the library pattern."""
        aminos = sg.return_available_residues()
        seq = self.lib_seq_holder.get().replace("x", "X")
        self.lib_seq_holder.set(seq)
        num_xs = self.lib_seq_holder.get().count("X")
        if (
            len(
                [
                    resi
                    for resi in self.lib_seq_holder.get()
                    if resi not in aminos and resi != "X"
                ]
            )
            == 0
        ):
            if num_xs > 0:
                return len(aminos) ** num_xs
            elif len(self.lib_seq_holder.get()) > 0:
                return 1
            else:
                return 0

    def countFilePeptides(self):
        line_num = None
        try:
            line_num = 0
            with open(self.seqs_from_file.get()) as pep_file:
                for line in pep_file:
                    if line.strip() and not line.startswith("#"):
                        line_num += 1
            return line_num
        except IOError:
            # Error message clickbox
            errmsg = "Could not read file %s" % self.seqs_from_file.get()
            tkinter.messagebox.showerror("File not found", errmsg)
            return 0

    def genLibrarySMILES(self):
        """
        Displays counts of types for SMILES library
        """
        if self.lib_len.get() == 0:
            pass
        else:
            # Set labels to counts as ... while counting
            blank = tk.StringVar()
            blank.set("...")
            self.lin_info_panel.config(textvariable=blank)
            self.SS_info_panel.config(textvariable=blank)
            self.HT_info_panel.config(textvariable=blank)
            self.SCSC_info_panel.config(textvariable=blank)
            self.SCNT_info_panel.config(textvariable=blank)
            self.SCCT_info_panel.config(textvariable=blank)

            # Decide whether to generate peptides from pattern or file
            if self.seqs_from_file.get():
                out_gen = sg.gen_library_from_file(
                    self.seqs_from_file.get(), ignore_errors=True
                )
            else:
                out_gen = sg.gen_structs_from_seqs(
                    sg.gen_all_matching_peptides(self.lib_seq_holder.get()),
                    ssbond=self.gen_SS.get(),
                    htbond=self.gen_HT.get(),
                    scctbond=self.gen_SCCT.get(),
                    scntbond=self.gen_SCNT.get(),
                    scscbond=self.gen_SCSC.get(),
                    linear=self.gen_linears.get(),
                )

            if self.synth_filter.get():
                out_gen = sg.filtered_output(out_gen, synthrules.synth_pass)
            if self.druglike_only.get():
                out_gen = sg.filtered_output(
                    out_gen, ca.lipinski_pass, key=operator.itemgetter(2)
                )

            # Count each type of constrained peptide, and update variables
            self.startWriting("COUNTING PEPTIDES")
            counts = sg.count_constraint_types(out_gen, ignore_errors=True)
            self.num_lin.set(counts["linear"])
            self.num_SS.set(counts["SS"])
            self.num_HT.set(counts["HT"])
            self.num_SCSC.set(counts["SCSC"])
            self.num_SCNT.set(counts["SCNT"])
            self.num_SCCT.set(counts["SCCT"])

            # print 'counted'
            self.doneWriting()
            # Set labels back to numbers
            self.lin_info_panel.config(textvariable=self.num_lin)
            self.SS_info_panel.config(textvariable=self.num_SS)
            self.HT_info_panel.config(textvariable=self.num_HT)
            self.SCSC_info_panel.config(textvariable=self.num_SCSC)
            self.SCNT_info_panel.config(textvariable=self.num_SCNT)
            self.SCCT_info_panel.config(textvariable=self.num_SCCT)

    def writeLibrary(self, format, one_file=False):
        """
        Writes the requested SMILES strings and names to a file
        """
        if self.lib_len.get() == 0:
            pass
        else:
            if self.seqs_from_file.get():
                out_gen = sg.gen_library_from_file(
                    self.seqs_from_file.get(), ignore_errors=True
                )
            else:
                peps = sg.gen_all_matching_peptides(self.lib_seq_holder.get())
                out_gen = sg.gen_structs_from_seqs(
                    peps,
                    ssbond=self.gen_SS.get(),
                    htbond=self.gen_HT.get(),
                    scctbond=self.gen_SCCT.get(),
                    scntbond=self.gen_SCNT.get(),
                    scscbond=self.gen_SCSC.get(),
                    linear=self.gen_linears.get(),
                )

            if self.synth_filter.get():
                out_gen = sg.filtered_output(
                    out_gen, synthrules.synth_pass, key=operator.itemgetter(0, 1, 2)
                )
            if self.druglike_only.get():
                out_gen = sg.filtered_output(
                    out_gen, ca.lipinski_pass, key=operator.itemgetter(2)
                )
            if self.nmethylate.get():
                out_gen = sg.nmethylate_peptides(out_gen)

            if format == "text":
                save = tkinter.filedialog.asksaveasfilename(
                    defaultextension=".txt",
                    initialdir=os.path.expanduser("~"),
                    initialfile="peptide_library",
                )
            elif one_file:
                save = tkinter.filedialog.asksaveasfilename(
                    defaultextension=".sdf",
                    initialdir=os.path.expanduser("~"),
                    initialfile="peptide_library",
                )
            else:
                save = tkinter.filedialog.askdirectory(
                    initialdir=os.path.expanduser("~")
                )
            self.startWriting("WRITING LIBRARY")
            try:
                if not save:
                    return
                if one_file:
                    out = sg.write_library(
                        out_gen, save, write=format, write_to_file=True
                    )
                else:
                    out = sg.write_library(out_gen, save, write=format)
            finally:
                self.doneWriting()
            if (out < self.lib_len.get()) and self.seqs_from_file.get():
                self.errorWritingFile()

    def writeLibrarySMILES(self):
        """'
        Calls self.writeLibrary to write the library to SMILES
        """
        self.writeLibrary("text")

    def writeLibrary3D(self):
        """'
        Calls self.writeLibrary to write the library to PDBs
        """
        self.writeLibrary("structure")

    def writeLibrary3Dtofile(self):
        """
        Calls self.writeLibrary to write the library to a single sdf
        """
        self.writeLibrary("structure", True)

    def drawLibrary(self):
        """
        Calls self.writeLibrary to write the library to 2D images
        """
        self.writeLibrary("draw")

    def startWriting(self, msg):
        """
        Indicator when the prgram is writing peptides
        """
        self.writelabel.config(text=msg, justify="center", relief="ridge")
        self.writelabel.update_idletasks()

    def doneWriting(self):
        """
        Remove writing indicator
        """
        # print 'DONE WRITE'
        self.writelabel.config(text="", relief="flat")

    def errorWritingFile(self):
        """If not all peptide lines in a file have been written"""
        msg = (
            "WARNING: Skipped writing some peptides in %s!" % self.seqs_from_file.get()
        )
        self.writelabel.config(text=msg, justify="center", relief="ridge")


class ZincChooser(tk.Toplevel):
    """
    GUI for adding ZINC non-natural amino-acids to the working amino library
    """

    def __init__(self, master):
        tk.Toplevel.__init__(self, master)
        self.resizable(0, 0)
        self.grid()
        # Title
        self.title("Choose ZINC amino-acids...")

        # labelframes for available zinc aminos
        self.zinc_available_lab = tk.LabelFrame(self, text="Available ZINC amino-acids")
        self.zinc_available_lab.grid(row=0, column=0)

        # labelframe for included zinc aminos
        self.zinc_included_lab = tk.LabelFrame(self, text="Included ZINC amino-acids")
        self.zinc_included_lab.grid(row=0, column=1)

        # Get ZINC peptides
        self.zinc_aminos = zinc_peptides.extract_peptides(zinc_def_file)

        # Available ZINC peptides variable
        self.available_zinc = tk.StringVar()
        self.available_zinc.set(" ".join(sorted(self.zinc_aminos.keys())))

        # Included ZINC peptides variable
        self.included_zinc = tk.StringVar()
        self.included_zinc.set(
            " ".join([amino for amino in sg.aminodata if amino not in sg.all_aminos])
        )

        # Set current listbox( included or available)
        self.cur_listbox = ""

        # Initialise
        self.zinc_listbox()
        self.zinc_included_listbox()
        self.add_remove_amino_button()
        self.structure_button()
        self.add_all_button()

    def zinc_listbox(self):
        """
        Listbox for choosing ZINC amino-acids
        """
        self.zinc_box = tk.Listbox(
            self.zinc_available_lab,
            activestyle="dotbox",
            selectmode=tk.SINGLE,
            listvariable=self.available_zinc,
            height=20,
        )
        self.zinc_box.grid()
        self.scroll_zinc_box = tk.Scrollbar(
            self.zinc_available_lab, orient=tk.VERTICAL, command=self.zinc_box.yview
        )
        self.zinc_box.config(yscrollcommand=self.scroll_zinc_box.set)
        self.scroll_zinc_box.grid(row=0, column=1, sticky=tk.N + tk.S)
        self.zinc_box.bind("<1>", self.set_cur)
        self.zinc_box.bind("<Double-Button-1>", self.add_remove_amino)

    def zinc_included_listbox(self):
        """
        Listbox for displaying ZINC amino-acids in library
        """
        self.zinc_included_box = tk.Listbox(
            self.zinc_included_lab,
            activestyle="dotbox",
            selectmode=tk.SINGLE,
            listvariable=self.included_zinc,
            height=20,
        )
        self.zinc_included_box.grid()
        self.scroll_included = tk.Scrollbar(
            self.zinc_included_lab,
            orient=tk.VERTICAL,
            command=self.zinc_included_box.yview,
        )
        self.zinc_included_box.config(yscrollcommand=self.scroll_included.set)
        self.scroll_included.grid(row=0, column=1, sticky=tk.N + tk.S)
        self.zinc_included_box.bind("<1>", self.set_cur)

    def set_cur(self, event):
        """
        Sets the current widget for adding/removing amino acids
        """
        self.cur_listbox = event.widget

    def structure_button(self):
        """
        Show a window with a drawing of the ZINC amino acid structure
        """
        self.structure_button = tk.Button(
            self, text="View selected amino", command=self.view_structure, padx="3.2m"
        )
        self.structure_button.grid(column=0, row=2, columnspan=2)

    def view_structure(self):
        """
        show_structure_button handler - wrapper around StructureWindow class
        """
        # print 'Called view structure'
        index = self.cur_listbox.index(tk.ACTIVE)
        name = self.cur_listbox.get(index)
        self.view_struct = StructureWindow(
            self, amino=self.zinc_aminos[name], name=name
        )

    def add_remove_amino_button(self):
        """
        Button to call the add-remove amino method
        """
        self.add_remove_button = tk.Button(
            self, text="Add/remove amino acid", command=self.add_remove_amino
        )
        self.add_remove_button.grid(column=0, row=1, columnspan=2)

    def add_all_button(self):
        """
        Button to add all ZINC aminos to working library
        """
        self.add_all_button = tk.Button(
            self, text="Add all ZINC amino-acids", command=self.add_all_aminos
        )
        self.add_all_button.grid(column=0, row=3, columnspan=2)

    def add_all_aminos(self):
        """
        Include all ZINC aminos in library
        """
        index = 0
        while True:
            name = self.zinc_box.get(index)
            if name not in sg.aminodata:
                sg.aminodata[name] = dict(self.zinc_aminos[name])
            if self.zinc_box.get(index) == self.zinc_box.get(tk.END):
                break
            index += 1
        cur_zincs = [amino for amino in sg.aminodata if amino not in sg.all_aminos]
        print(("Number of zincs in library", len(cur_zincs)))
        self.included_zinc.set(" ".join(sorted(cur_zincs)))

    def add_remove_amino(self, event=""):
        """
        Add or remove the selected ZINC amino-acid to/from the working library
        """
        cur_index = self.cur_listbox.index(tk.ACTIVE)
        name = self.cur_listbox.get(cur_index)
        # print cur_index, name
        if name in sg.aminodata:
            del sg.aminodata[name]
            cur_zincs = [amino for amino in sg.aminodata if amino not in sg.all_aminos]
            # print 'cur_zincs',cur_zincs
            self.included_zinc.set(" ".join(cur_zincs))

            # print [name for name in sg.aminodata]
            # print self.included_zinc.get()
            # print self.available_zinc.get()
        else:
            sg.aminodata[name] = dict(self.zinc_aminos[name])
            cur_zincs = [amino for amino in sg.aminodata if amino not in sg.all_aminos]
            # print 'cur_zincs', cur_zincs
            self.included_zinc.set(" ".join(cur_zincs))

            # print [name for name in sg.aminodata]
            # print self.included_zinc.get()
            # print self.available_zinc.get()


class StructureWindow(tk.Toplevel):
    """
    Window to display 2D structure of an amino-acid / peptide

    aminodef should contain information
    """

    def __init__(self, master, **kwargs):
        tk.Toplevel.__init__(self, master)
        self.temp_holder, self.temp_filename = tempfile.mkstemp()
        # print kwargs
        if "amino" in kwargs:
            try:
                self.name = kwargs["name"]
                self.smiles = kwargs["amino"]["SMILES"]
                self.title(self.name)
                self.display_structure()
            except KeyError:
                self.destroy()
        else:
            self.destroy()

    def display_structure(self):
        """
        Draw the amino acid and display in the window
        """
        self.mol_obj = Chem.MolFromSmiles(self.smiles)
        # self.mol_obj.draw(show=False, filename=self.temp_filename)
        image = Draw.MolToImage(self.mol_obj)
        # image = PIL.open(self.temp_filename)
        self.image = image  # Save a reference, so it can be saved
        self.imagedata = piltk.PhotoImage(image)
        self.image_label = tk.Label(self, image=self.imagedata)
        self.image_label.grid()

    def destroy(self):
        """
        Clean up tempfile
        """
        tk.Toplevel.destroy(self)
        try:
            os.close(self.temp_handle)
            os.remove(self.temp_filename)
        except:
            pass


class PeptideBuilder(tk.Toplevel):
    """
    Window that allows building a peptide by selecting residues from a listbox
    """

    def __init__(self, master, **options):
        self.seq = options["seq_holder"]
        del options["seq_holder"]
        # print 'Master', master
        tk.Toplevel.__init__(self, master, **options)
        # Variables
        self.aminos = tk.StringVar()
        self.aminos.set(" ".join(list(sg.aminodata.keys())))
        self.AminoListbox()
        self.PeptideBox()
        self.AddButton()
        self.ViewStructButton()
        # self.DoneButton()

    def AminoListbox(self):
        """
        Listbox of all aminos in sg.aminodata
        """
        self.amino_lab = tk.LabelFrame(self, text="Available amino acids")
        self.amino_lab.grid(row=0, column=0)
        self.amino_box = tk.Listbox(
            self.amino_lab,
            activestyle="dotbox",
            selectmode=tk.SINGLE,
            listvariable=self.aminos,
            height=20,
        )
        self.amino_box.grid()
        self.scroll_amino_box = tk.Scrollbar(
            self.amino_lab, orient=tk.VERTICAL, command=self.amino_box.yview
        )
        self.amino_box.config(yscrollcommand=self.scroll_amino_box.set)
        self.scroll_amino_box.grid(row=0, column=1, sticky=tk.N + tk.S)
        self.amino_box.bind("<Double-Button-1>", self.AddPeptideBox)

    def PeptideBox(self):
        """
        Displays sequence on added residues as text
        """
        self.pep_lab = tk.LabelFrame(self, text="Current Peptide")
        self.pep_lab.grid(row=0, column=1)
        self.display_box = tk.Text(self.pep_lab, bg="white", wrap=tk.WORD, width=30)
        self.display_box.grid()

    def AddPeptideBox(self, event=""):
        """
        Adds currently selected amino-acid to peptide box.

        Will update the sequence holder
        """
        index = self.amino_box.index(tk.ACTIVE)
        name = self.amino_box.get(index)

        if len(self.display_box.get(1.0, "end")) == 1:
            self.display_box.insert("end", name)
        else:
            self.display_box.insert("end", "," + name)
        self.seq.set(self.display_box.get(1.0, "end"))
        try:
            self.master.SmilesHandler()
        except:
            "Failed calling master SmilesHandler"

    def AddButton(self):
        """
        Button to call Add to peptide box
        """
        self.add_button = tk.Button(
            self, text="Add amino-acid", command=self.AddPeptideBox
        )
        self.add_button.grid(row=1)

    def ViewAminoStruct(self):
        """
        Draws the selceted amino-acis in a window
        """

        # print 'Called view '
        index = self.amino_box.index(tk.ACTIVE)
        name = self.amino_box.get(index)
        self.view_struct = StructureWindow(self, amino=sg.aminodata[name], name=name)

    def ViewStructButton(self):
        """
        Button to call ViewAminoStruct
        """
        self.struct_button = tk.Button(
            self, text="View amino-acid", command=self.ViewAminoStruct
        )
        self.struct_button.grid(row=2)


if __name__ == "__main__":
    app = PepApp()
    app.mainloop()
