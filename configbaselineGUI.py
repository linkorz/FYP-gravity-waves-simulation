#import relevant modules, functions
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
#User inputs for dx,dz, SpongeHeight
dxlabel = Label(root, text="dx")
dxlabel.grid(row=7,column=0)
dxent = Entry(root)
dxent.insert(0,"500")
dxent.grid(row=7,column=1)

dzlabel = Label(root, text="dz")
dzlabel.grid(row=8,column=0)
dzent = Entry(root)
dzent.insert(0,"500")
dzent.grid(row=8,column=1)

SHlabel = Label(root, text="SpongeHeight")
SHlabel.grid(row=9,column=0)
SHent = Entry(root)
SHent.insert(0,"20000")
SHent.grid(row=9,column=1)
#--------------------------------------#
#Gaussian wind shear inputs
u_maxlabel = Label(root, text="u_max")
u_maxlabel.grid(row=10,column=0)
u_maxent = Entry(root)
u_maxent.insert(0,"0")
u_maxent.grid(row=10,column=1)

u_zloclabel = Label(root, text="u_zloc")
u_zloclabel.grid(row=11,column=0)
u_zlocent = Entry(root)
u_zlocent.insert(0,"100000")
u_zlocent.grid(row=11,column=1)

u_siglabel = Label(root, text="u_sig")
u_siglabel.grid(row=12,column=0)
u_sigent = Entry(root)
u_sigent.insert(0,"10000")
u_sigent.grid(row=12,column=1)
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
    return Tmin, Tmax, skipT

def saveX():
    Xminval= X1ent.get()
    Xmaxval= X2ent.get()
    global Xmin,Xmax
    Xmin= int(Xminval)
    Xmax= int(Xmaxval)
    print(Xmin,Xmax)
    return Xmin, Xmax
   
def saveZ():
    Zminval= Z1ent.get()
    Zmaxval= Z2ent.get()
    global Zmin, Zmax
    Zmin= int(Zminval)
    Zmax= int(Zmaxval)
    print(Zmin,Zmax)
    return Zmin, Zmax

def savedxdzSH():
    dxval= dxent.get()
    dzval= dzent.get()
    SHval= SHent.get()
    global dx, dz, SpongeHeight
    dx= int(dxval)
    dz= int(dzval)
    SpongeHeight= int(SHval)
    print(dx,dz,SpongeHeight)
    return dx, dz, SpongeHeight

def savewind():
    u_maxval= u_maxent.get()
    u_zlocval= u_zlocent.get()
    u_sigval= u_sigent.get()
    global u_max, u_zloc, u_sig
    u_max= int(u_maxval)
    u_zloc= int(u_zlocval)
    u_sig= int(u_sigval)
    print(u_max,u_zloc,u_sig)
    return u_max, u_zloc, u_sig

def runfile():
    exec(open('GW_2D_propagation.py').read())
#--------------------------------------#
#Buttons
saveTval= Button(root, text="Save", command=saveT) #Save button for T inputs
saveTval.grid(row=2,column=3)

saveXval= Button(root, text="Save", command=saveX) #Save button for X inputs
saveXval.grid(row=4,column=3)

saveZval= Button(root, text="Save", command=saveZ) #Save button for Z inputs
saveZval.grid(row=6,column=3)

savedxdzSHval= Button(root, text="Save", command=savedxdzSH) #Save button for dx, dz, and SpongeHeight
savedxdzSHval.grid(row=9,column=3)

savewindval= Button(root, text="Save", command=savewind) #Save button for Z inputs
savewindval.grid(row=12,column=3)

GW= Button(root, text="GW_2D_propagation", command=runfile) #Button to run the main script
GW.grid(row=25,column=1)

root.mainloop()


