import tkinter as tk
from tkinter import ttk

answer = 0

root = tk.Tk()

def button_click(num1, num2, num3):
    num3.set(num1.get() + num2.get())

window_height = 600
window_width = 800

screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

centre_x = int(screen_width/2 - window_width/2)
centre_y = int(screen_height/2 - window_height/2)

root.geometry(f'{window_width}x{window_height}+{centre_x}+{centre_y}')
root.resizable(False, False)

number1 = tk.IntVar()
number2 = tk.IntVar()
number3 = tk.IntVar()

number1_textbox = ttk.Entry(root, textvariable=number1)
number1_textbox.pack()

number2_textbox = ttk.Entry(root, textvariable=number2)
number2_textbox.pack()

number3_textbox = ttk.Label(root, textvariable=number3)
number3_textbox.pack()

button = ttk.Button(root, text="add numbers", command=lambda: button_click(number1, number2, number3))
button.pack()

root.mainloop()
