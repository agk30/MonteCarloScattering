import tkinter as tk
import os
import threading

window = tk.Tk()

class myThread (threading.Thread):
    def __init__(self, threadID, name, counter):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.counter = counter
    def run(self):
        run_fortran()
    def exit(self):
        stop_fortran()


def run_fortran():
    os.system('D:/Dev/MonteCarloScattering/build/MCScattering.exe')
def stop_fortran():
    os.system('taskkill /F /IM MCScattering.exe')

thread1 = myThread(2, "Thread-1", 2)

button_start = tk.Button(master=window, text='Run', width=25,
    height=5,
    command=thread1.start)
button_start.pack()

button_stop = tk.Button(master=window, text='Stop', width=25,
    height=5,
    command=thread1.exit)
button_stop.pack()


window.mainloop()
