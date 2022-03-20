from tkinter import *

root = Tk()

clicked = StringVar()
clicked.set("Earth_iso")

model = OptionMenu(root, clicked, "Earth_iso", "Earth_MSIS")
model.grid(row=3,column=1)


Tminlabel = Label(root, text="Tmin")
Tminlabel.grid(row=0,column=0)
Tminent = Entry(root)
Tminent.insert(0,"0")
Tminent.grid(row=0,column=1)


def myClick():
    Tminval= Tminent.get()
    Tmin= int(Tminval)
    Tmin

def runfile():
    exec(GW_2D_propagation PYTHON.py)

save1= Button(root, text="Save", command=myClick)
save1.grid(row=0,column=3)
GW= Button(root, text="GW_2D_propagation", command=runfile)
GW.grid(row=25,column=1)

root.mainloop()