#import relevant modules, functions
from tkinter import *
from PIL import ImageTk,Image


root = Tk()
root.title("MAGNUS-P Master")

planetimage = ImageTk.PhotoImage(Image.open("Planets.JFIF"))
planetimagelabel = Label(image=planetimage)
planetimagelabel.grid(columnspan=20)

#model select------------------#
ModelLabel = Label(root, text="Atmospheric")
ModelLabel.grid(row=22,column=5)
ModelLabel.configure(font=("Arial",20))

ModelLabel1 = Label(root, text="Model")
ModelLabel1.grid(row=23,column=5)
ModelLabel1.configure(font=("Arial",20))

model = StringVar()
model.set("Earth_Isothermal")

modelinput = OptionMenu(root, model, "Earth_Isothermal", "Earth_MSIS")
modelinput.grid(row=24,column=5)
#------------------------------#
#Use input values for T
TimeSettingsLabel = Label(root, text="Time Settings")
TimeSettingsLabel.grid(row=1,column=1)
TimeSettingsLabel.configure(font=("Arial",20))

T1label = Label(root, text="Tmin")
T1label.grid(row=2,column=0)
T1ent = Entry(root)
T1ent.insert(0,"0")
T1ent.grid(row=2,column=1)

T2label = Label(root, text="Tmax")
T2label.grid(row=3,column=0)
T2ent = Entry(root)
T2ent.insert(0,"5000")
T2ent.grid(row=3,column=1)

TSlabel = Label(root, text="skipT")
TSlabel.grid(row=4,column=0)
TSent = Entry(root)
TSent.insert(0,"30")
TSent.grid(row=4,column=1)
#--------------------------------------#
#User inputs for X
DomainSettingsLabel = Label(root, text="Domain Settings")
DomainSettingsLabel.grid(row=6,column=1)
DomainSettingsLabel.configure(font=("Arial",20))

X1label = Label(root, text="Xmin")
X1label.grid(row=7,column=0)
X1ent = Entry(root)
X1ent.insert(0,"-20000")
X1ent.grid(row=7,column=1)

X2label = Label(root, text="Xmax")
X2label.grid(row=8,column=0)
X2ent = Entry(root)
X2ent.insert(0,"20000")
X2ent.grid(row=8,column=1)
#--------------------------------------#
#--------------------------------------#
#User inputs for Z
Z1label = Label(root, text="Zmin")
Z1label.grid(row=9,column=0)
Z1ent = Entry(root)
Z1ent.insert(0,"0")
Z1ent.grid(row=9,column=1)

Z2label = Label(root, text="Zmax")
Z2label.grid(row=10,column=0)
Z2ent = Entry(root)
Z2ent.insert(0,"160000")
Z2ent.grid(row=10,column=1)
#--------------------------------------#
#User inputs for dx,dz, SpongeHeight
ResolutionsLabel = Label(root, text="Resolutions")
ResolutionsLabel.grid(row=1,column=5)
ResolutionsLabel.configure(font=("Arial",20))

dxlabel = Label(root, text="dx")
dxlabel.grid(row=2,column=4)
dxent = Entry(root)
dxent.insert(0,"500")
dxent.grid(row=2,column=5)

dzlabel = Label(root, text="dz")
dzlabel.grid(row=3,column=4)
dzent = Entry(root)
dzent.insert(0,"500")
dzent.grid(row=3,column=5)

SHlabel = Label(root, text="SpongeHeight")
SHlabel.grid(row=4,column=4)
SHent = Entry(root)
SHent.insert(0,"20000")
SHent.grid(row=4,column=5)
#--------------------------------------#
#Gaussian wind shear inputs
WindLabel = Label(root, text="Wind")
WindLabel.grid(row=6,column=9)
WindLabel.configure(font=("Arial",20))

u_maxlabel = Label(root, text="u_max")
u_maxlabel.grid(row=7,column=8)
u_maxent = Entry(root)
u_maxent.insert(0,"0")
u_maxent.grid(row=7,column=9)

u_zloclabel = Label(root, text="u_zloc")
u_zloclabel.grid(row=8,column=8)
u_zlocent = Entry(root)
u_zlocent.insert(0,"100000")
u_zlocent.grid(row=8,column=9)

u_siglabel = Label(root, text="u_sig")
u_siglabel.grid(row=9,column=8)
u_sigent = Entry(root)
u_sigent.insert(0,"10000")
u_sigent.grid(row=9,column=9)
#--------------------------------------#
#Flags inputs
FlagsLabel = Label(root, text="Flags")
FlagsLabel.grid(row=6,column=5)
FlagsLabel.configure(font=("Arial",20))

IsTopSpongeLayer = IntVar()
IsTopSpongeLayer.set("0")
Radiobutton(root, text="Off", variable=IsTopSpongeLayer, value=0).grid(row=7,column=5)
Radiobutton(root, text="On", variable=IsTopSpongeLayer, value=1).grid(row=7,column=6)
IsTopSpongeLayerlabel = Label(root, text="IsTopSpongeLayer")
IsTopSpongeLayerlabel.grid(row=7,column=4)

IsViscosity = IntVar()
IsViscosity.set("1")
Radiobutton(root, text="Off", variable=IsViscosity, value=0).grid(row=8,column=5)
Radiobutton(root, text="On", variable=IsViscosity, value=1).grid(row=8,column=6)
IsViscositylabel = Label(root, text="IsViscosity")
IsViscositylabel.grid(row=8,column=4)

IsConduction = IntVar()
IsConduction.set("1")
Radiobutton(root, text="Off", variable=IsConduction, value=0).grid(row=9,column=5)
Radiobutton(root, text="On", variable=IsConduction, value=1).grid(row=9,column=6)
IsConductionlabel = Label(root, text="IsConduction")
IsConductionlabel.grid(row=9,column=4)

IsDiffusionImplicit = IntVar()
IsDiffusionImplicit.set("0")
Radiobutton(root, text="Off", variable=IsDiffusionImplicit, value=0).grid(row=10,column=5)
Radiobutton(root, text="On", variable=IsDiffusionImplicit, value=1).grid(row=10,column=6)
IsDiffusionImplicitlabel = Label(root, text="IsDiffusionImplicit")
IsDiffusionImplicitlabel.grid(row=10,column=4)
#---------------------------------------#
#Forcing inputs
ForcingLabel = Label(root, text="Forcing")
ForcingLabel.grid(row=1,column=9)
ForcingLabel.configure(font=("Arial",20))

ForcingThermal = BooleanVar()
ForcingThermal.set(False)
Radiobutton(root, text="True", variable=ForcingThermal, value=True).grid(row=2,column=9)
Radiobutton(root, text="False", variable=ForcingThermal, value=False).grid(row=2,column=10)
ForcingThermallabel = Label(root, text="ForcingThermal")
ForcingThermallabel.grid(row=2,column=8)

ForcingVerticalVelocity = StringVar()
ForcingVerticalVelocity.set(True)
Radiobutton(root, text="True", variable=ForcingVerticalVelocity, value=True).grid(row=3,column=9)
Radiobutton(root, text="False", variable=ForcingVerticalVelocity, value=False).grid(row=3,column=10)
ForcingVerticalVelocitylabel = Label(root, text="ForcingVerticalVelocity")
ForcingVerticalVelocitylabel.grid(row=3,column=8)
#---------------------------------#
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

def saveDomain():
    Xminval= X1ent.get()
    Xmaxval= X2ent.get()
    Zminval= Z1ent.get()
    Zmaxval= Z2ent.get()
    global Xmin,Xmax,Zmin,Zmax
    Xmin= int(Xminval)
    Xmax= int(Xmaxval)
    Zmin= int(Zminval)
    Zmax= int(Zmaxval)
    print(Xmin,Xmax,Zmin,Zmax)
    return Xmin, Xmax, Zmin, Zmax
   

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
saveTval.grid(row=5,column=1)

saveDomainval= Button(root, text="Save", command=saveDomain) #Save button for Domain settings inputs
saveDomainval.grid(row=11,column=1)

savedxdzSHval= Button(root, text="Save", command=savedxdzSH) #Save button for dx, dz, and SpongeHeight
savedxdzSHval.grid(row=5,column=5)

savewindval= Button(root, text="Save", command=savewind) #Save button for Z inputs
savewindval.grid(row=10,column=9)



GW= Button(root, text="Run MAGNUS-P", command=runfile) #Button to run the main script
GW.grid(row=25,column=5)

root.mainloop()


