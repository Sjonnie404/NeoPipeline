from tkinter import *


def btn_clicked():
    print("Button Clicked")


window = Tk()

window.geometry("1440x1024")
window.configure(bg = "#ffffff")
canvas = Canvas(
    window,
    bg = "#ffffff",
    height = 1024,
    width = 1440,
    bd = 0,
    highlightthickness = 0,
    relief = "ridge")
canvas.place(x = 0, y = 0)

img0 = PhotoImage(file =f"../../Images/img0.png")
b0 = Button(
    image = img0,
    borderwidth = 0,
    highlightthickness = 0,
    command = btn_clicked,
    relief = "flat")

b0.place(
    x = 261, y = 236,
    width = 209,
    height = 104)

background_img = PhotoImage(file =f"../../Images/background.png")
background = canvas.create_image(
    -209.5, 98.0,
    image=background_img)

canvas.create_text(
    325.5, -308.5,
    text = "This is a test",
    fill = "#000000",
    font = ("None", int(48.0)))

canvas.create_text(
    79.5, -155.0,
    text = "Your text should appear here:",
    fill = "#000000",
    font = ("None", int(48.0)))

entry0_img = PhotoImage(file =f"../../Images/img_textBox0.png")
entry0_bg = canvas.create_image(
    575.0, -142.0,
    image = entry0_img)

entry0 = Entry(
    bd = 0,
    bg = "#c4c4c4",
    highlightthickness = 0)

entry0.place(
    x = 437, y = -176,
    width = 276,
    height = 66)

window.resizable(False, False)
window.mainloop()
