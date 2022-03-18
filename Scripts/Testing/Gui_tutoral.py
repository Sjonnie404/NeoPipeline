from tkinter import *
from tkinter.ttk import Combobox
from tkinter.ttk import Progressbar
from tkinter.scrolledtext import *
from tkinter import messagebox
from tkinter import filedialog
from os import path
from tkinter import Menu
from tkinter import ttk

# TODO: Write code in modular fassion

def button1_clicked():
    input = textbox1.get()
    label1.configure(text='your input: ' + input)
    print('clicked button1')
    combobox1_output = combobox1.get()
    print('Selected value: ' + combobox1_output)
    print(checkbox1_state.get())
    print(checkbox2_state.get())
    if sel.get() == 1:                  # Note get will be integers, so for strings we need elif tree.
        print('selected radiobutton: A')
    elif sel.get() == 2:
        print('selected radiobutton: B')
        outcome = messagebox.askyesnocancel('Question for you', 'Do you like testing?')  # Note that the code gets frozen until you answer
        print(outcome)
    elif sel.get() == 3:
        print('selected radiobutton: C')
        messagebox.showinfo('Testing this box', 'Is it actually showing?!')


def rad1_clicked():
    print("User selected 'A'")


window = Tk()

window.title("Shane's app")
window.geometry('640x360')

tab_control = ttk.Notebook(window)  # Adds tab control, NOTE: widgets should be stored on the tabs and not the window.
tab1 = ttk.Frame(tab_control)
tab_control.add(tab1, text='First')
tab2 = ttk.Frame(tab_control)
tab_control.add(tab2, text='Second')
tab_control.pack(expand=1, fill='both')

textbox1 = Entry(tab1, width=25)
textbox1.grid(column=1, row=0)
textbox1.focus()

button1 = Button(tab1, text='Click Me',
                 bg='grey',
                 command=button1_clicked)  # Note that () shouldn't be added, else the function triggers on startup
button1.grid(column=2, row=0)

label1 = Label(tab1, text='Hello',
               font=("Commic Sans", 30))
label1.grid(column=3, row=1)

textbox2 = Entry(tab1, width=25, state='disabled')
textbox2.grid(column=1, row=5)
textbox2.focus()  # focus is put on the first call.

combobox1 = Combobox(tab1, state='readonly')
combobox1['values'] = (1, 2, 3, 4, 5)
#combobox1.current(0)
combobox1.grid(column=2, row=5)

checkbox1_state = BooleanVar()
checkbox1_state.set(True)
checkbox1 = Checkbutton(tab1, text='Please Choose', var=checkbox1_state)  # Quite ambiguous if you ask me.
checkbox1.grid(column=3, row=5)

checkbox2_state = BooleanVar()
checkbox2_state.set(True)
checkbox2 = Checkbutton(tab1, text='Please Choose again',var=checkbox2_state)
checkbox2.grid(column=3, row=6)

sel = IntVar()
rad1 = Radiobutton(tab1, text="A", value=1, variable=sel, command=rad1_clicked)  # Why do logic when user is selecting?
rad2 = Radiobutton(tab1, text="B", value=2, variable=sel)
rad3 = Radiobutton(tab1, text="C", value=3, variable=sel)  # Note values need to be integers.
rad1.grid(column=1, row=7)
rad2.grid(column=2, row=7)
rad3.grid(column=3, row=7)

textArea1 = ScrolledText(tab1, width=30, height=5)
textArea1.insert(INSERT, 'If you have any comments,\n please place them here')  # Can also use escape characters.
# textArea1.delete(1.0, END)  # Use this to clear text area's
textArea1.grid(column=1, row=8)

spin_default = IntVar()
spin_default.set(58)  # This doesn't seem to work.
spinBox1 = Spinbox(tab1, values=(57, 58, 1928), width=10, textvariable=spin_default)  # Or use a range: from_=57, to=1928
spinBox1.grid(column=1, row=9)

bar = Progressbar(tab1, length=200)
bar.grid(column=1, row=10)
bar['value'] = 100  # Setting the progressbar based on percentage.

#file = filedialog.askopenfilename()  # This automatically runs, and should be put under button logic.
# files = filedialog.askopenfilenames()  # This works with multiple files.

# Specify specific filetypes, good practice is to always have the 'all types' option.
# file = filedialog.askopenfilename(filetypes = (("Text files","*.txt"),("all files","*.*")))
#directory = filedialog.askdirectory()  # Can also be done with directories.
#file1 = filedialog.askopenfilename(initialdir= path.dirname(__file__))  # can also specify initial filepath

menu = Menu(window)
new_item = Menu(menu, tearoff=0)
new_item.add_command(label='New', command=button1_clicked())
new_item.add_separator()
new_item.add_command(label='Edit')
menu.add_cascade(label='File', menu=new_item)
window.config(menu=menu)


window.mainloop()