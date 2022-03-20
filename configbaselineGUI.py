from tkinter import *

root = Tk()

#model select------------------#
clicked = StringVar()
clicked.set("Earth_Isothermal")

model = OptionMenu(root, clicked, "Earth_Isothermal", "Earth_MSIS")
model.grid(row=24,column=1)
#------------------------------#
#Use input values for T
T1label = Label(root, text="Tmin")
T1label.grid(row=0,column=0)
T1ent = Entry(root)
T1ent.insert(0,"0")
T1ent.grid(row=0,column=1)

T2label = Label(root, text="Tmax")
T2label.grid(row=1,column=0)
T2ent = Entry(root)
T2ent.insert(0,"5000")
T2ent.grid(row=1,column=1)

TSlabel = Label(root, text="skipT")
TSlabel.grid(row=2,column=0)
TSent = Entry(root)
TSent.insert(0,"30")
TSent.grid(row=2,column=1)
#--------------------------------------#
#User inputs for X
X1label = Label(root, text="Xmin")
X1label.grid(row=3,column=0)
X1ent = Entry(root)
X1ent.insert(0,"-20000")
X1ent.grid(row=3,column=1)

X2label = Label(root, text="Xmax")
X2label.grid(row=4,column=0)
X2ent = Entry(root)
X2ent.insert(0,"20000")
X2ent.grid(row=4,column=1)
#--------------------------------------#
#--------------------------------------#
#User inputs for Z
Z1label = Label(root, text="Zmin")
Z1label.grid(row=5,column=0)
Z1ent = Entry(root)
Z1ent.insert(0,"0")
Z1ent.grid(row=5,column=1)

Z2label = Label(root, text="Zmax")
Z2label.grid(row=6,column=0)
Z2ent = Entry(root)
Z2ent.insert(0,"160000")
Z2ent.grid(row=6,column=1)
#--------------------------------------#
#Button Functions
def saveT():
    Tminval= T1ent.get()
    Tmaxval= T2ent.get()
    skipTval= TSent.get()
    global Tmin, Tmax, skipT
    Tmin= int(Tminval)  #Initial time
    Tmax= int(Tmaxval)  #Final time in seconds
    skipT= int(skipTval)    #Number of seconds to skip string results
    print(Tmin,Tmax,skipT)

def saveX():
    Xminval= X1ent.get()
    Xmaxval= X2ent.get()
    global Xmin,Xmax
    Xmin= int(Xminval)
    Xmax= int(Xmaxval)
    print(Xmin,Xmax)
   
def saveZ():
    Zminval= Z1ent.get()
    Zmaxval= Z2ent.get()
    global Zmin, Zmax
    Zmin= int(Zminval)
    Zmax= int(Zmaxval)
    print(Zmin,Zmax)

def runfile():
    exec(open('GW_2D_propagation.py').read())
#--------------------------------------#
#Buttons
saveTval= Button(root, text="Save", command=saveT)
saveTval.grid(row=2,column=3)

saveXval= Button(root, text="Save", command=saveX)
saveXval.grid(row=4,column=3)

saveZval= Button(root, text="Save", command=saveZ)
saveZval.grid(row=6,column=3)



GW= Button(root, text="GW_2D_propagation", command=runfile)
GW.grid(row=25,column=1)

root.mainloop()