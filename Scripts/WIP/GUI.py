########################################################################################
#  Script to generate GUI to use
#  Input: Count file that contains Ensemble IDs
#  Output: Fasta file with headers, containing Ensemble information and NCBI sequence
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
########################################################################################
import tkinter

from Fetch_transcript_seq import *
from tkinter import *
from tkinter.ttk import Combobox
from tkinter.ttk import Progressbar
from tkinter.scrolledtext import *
from tkinter import messagebox
from tkinter import filedialog
from os import path
from tkinter import Menu
from tkinter import ttk

window = Tk()
# window.wm_attributes("-transparentcolor", 'grey')
window.title("Neo-antigen pipeline")
window.geometry('640x360')

tab_control = ttk.Notebook(window)  # Adds tab control, NOTE: widgets should be stored on the tabs and not the window.
tab1 = ttk.Frame(tab_control)
tab_control.add(tab1, text='Import')
tab2 = ttk.Frame(tab_control)
tab_control.add(tab2, text='Process')
tab3 = ttk.Frame(tab_control)
tab_control.add(tab3, text='Save')
tab4 = ttk.Frame(tab_control)
tab_control.add(tab4, text='Info')
tab_control.pack(expand=1, fill='both')

label1_1 = Label(tab4, text='-\tThank you for using the Neo-antigen pipeline App.', anchor='w') # TODO should be padded to the left
label1_1.grid(column=0, row=0, padx=10, pady=[10, 0])

# TODO This should be labels.
textArea1_2 = Text(tab4, width=50, height=3, bg='#f0f0f0')
textArea1_2.insert(INSERT, 'This application has been commissioned by the \nTheoretical Biology & Bioinformatics group at \nUtrecht University.')
textArea1_2.grid(column=0, row=1, padx=10)
textArea1_2.config(state=tkinter.DISABLED)
sep = ttk.Separator(tab4).grid(sticky='ew', pady=[5,10])




# label1_2 = Label(tab1, text="""This application has been commissioned by
#                                the Theoretical Biology & Bioinformatics group
#                                at Utrecht University""")
# label1_2.grid(column=1, row=1, padx=10, pady=10)



window.mainloop()