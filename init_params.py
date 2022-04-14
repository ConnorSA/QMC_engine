## @package init_params
## @brief Python functions and scripts to obtain user input
# import sympy as sym
import tkinter as tk
# from PIL import ImageTk, Image
# import os
import numpy as np


## @brief Function to get a numerical value/parameter through command line.
## This is used when a gui cant be used or is unsuitable. \n
## Example code to get the bond length: \n
## x=get_input(parameter='bond length')
## @param parameter The name of the parameter/variable needed from the user.
def get_input(parameter='parameter'):
    while True:
        x=input('Please enter a value for the {0}: '.format(parameter))
        try:
            x=float(x)
            break
        except:
            print('You did not enter a proper numerical value for the {0}'.format(parameter))
    return x


## @param title GUI window title
title='User Input'

## @param bg_color background color of the gui window
bg_color='white'
## @param fg_color foreground color of the gui window
fg_color='grey'
## note: Mac OS may not support changing button color
## @param bg_button background color of the close button
bg_button='grey'
## @param fg_color foreground color of the close button
fg_button='black'

## @param font Font of the GUI text
font='Times New Roman'
## @param fontsize Size of GUI font
fontsize=18

## @param length Length of the GUI window
length=1000
## @param height Height of the GUI window
height=1000


#This creates the main window of the GUI
window=tk.Tk()
window.title(title)
window.geometry('{0}x{1}'.format(length,height))
window.configure(background=bg_color)

# Use this code if latex is needed in the gui
# This uses sympy to create an image of the latex for the GUI

# sym.init_printing(backcolor=color,forecolor=color)

# Name of sympy output image/location
# image1='./pics/image1.png'

# f1=r'Please enter the distance between two atoms. \\'
# f2=r'If no distance is entered, 1 {\AA} will be used. \\'
# f3=r'Press the "Close" button once done. '
# sym.preview(f1+f2+f3,viewer='file',filename=image1,euler=False)

# #Creates a Tkinter-compatible photo image, then used in the GUI window
# pic=image1
# img=ImageTk.PhotoImage(Image.open(pic))
# panel=tk.Label(window,image=img)
# panel.grid(column=0) #,row=0,padx=0,pady=5)

## @param txt0 The text shown in the top row of the GUI
txt0="The default simulation parameters are shown in the boxes in \
the right-hand column below. Please change any to your desired simulation \
parameters and then click 'Close' at the bottom."

## Initializing/Configuring the top row text
text0=tk.Label(window,text=txt0,font=(font,fontsize),wraplength=length-5) #,justify='center')
text0.configure(background=bg_color)
text0.grid(row=0,column=0,columnspan=2,sticky='w') #column=0,padx=5,pady=5)

## @brief Function for creating a row of user entry within the GUI
## @param window The tkinter window object that has been initialized
## @param tk_text Variable name (as string) for holding tkinter label
## @param tk_entry Variable name (as string) for holding tkinter entry object
## @param tk_input  Variable name (as string) for displaying the default parameter value
## @param row The row of the GUI to display the prompt
## @param dval Default value of the parameter/variable
## @param text Name of the parameter/variable
## @param font The name of the font to be used as given by \ref basis_functions
## @param fontsize The size of the font
def new_param(window,
              tk_text,tk_entry,tk_input,
              row,
              dval='1',
              text='param',
              font=font,fs=fontsize):

    # Text to ask user for a specific input
    vars()[tk_text]=tk.Label(window,text=text,font=(font,fs))
    vars()[tk_text].configure(background=bg_color)
    vars()[tk_text].grid(row=row,pady=5,column=0,sticky='e')  #padx=5,

    # Create a space for user text entry
    # Only need first line if you want to show a default value
    vars()[tk_input]=tk.StringVar(value=dval)
    vars()[tk_entry]=tk.Entry(window,width=9,textvariable=vars()[tk_input])
    vars()[tk_entry].grid(column=1,row=row,pady=5,sticky='w')
    return vars()[tk_text],vars()[tk_entry],vars()[tk_input]


## @param params The names of the parameters needed, given in the order of entry.
params=['Number of electrons','Number of atoms','Bond length if there are two atoms', \
'Number of trials for latin hypercube search','Total number of steps for MCMC run', \
'Biased optimizer','Number of steps', \
'Type of orbitals', \
'Proton number', \
'Number of OMP Threads', \
'Random seed for the latin hypecube', \
'Number of linear terms in the single-electron wavefunction per atom', \
'Number of dofs in the Jastrow (interaction) term', \
'Lengthscale for the finite difference method', \
'Inverse lengthscale of nuclear-electron interaction', \
'Inverse lengthscale of electron-electron interaction', \
'Amount of distance away from atoms to plot', \
'Number of points along axis in each direction to plot']

## @param def_values Default parameter values carried over from latin_driver.f90
def_values=[1,1,1.5,40,1000000,0,10,100,1,10,12345,1,7,0.01,1.0,1.0,5.0,20]

## @param floats Indices of parameters that are floats
floats=[2, 13, 14, 15, 16]
## @param integers Indices of parameters that are iNtegers
integers=[0, 1, 3, 4, 5, 6, 7,8, 9, 10, 11, 12, 17]

## @param N_params Number of parameters/variables
N_params=len(params)

txt1="The user must provide the following inputs: "
text1=tk.Label(window,text=txt1,font=(font,fontsize),wraplength=length-10)
text1.configure(background=bg_color)
text1.grid(row=1,columnspan=3,sticky='w') # pady=5,padx=5

for i,param in enumerate(params[:5]):
    vars()['text'+str(i)],vars()['entry'+str(i)],vars()['input'+str(i)]= \
    new_param(window,'text'+str(i),'entry'+str(i),'input'+str(i),row=i+2, \
    dval=str(def_values[i]),text=param+': ')


txt2="Biased Optimizer Options (not enabled by default):"
text2=tk.Label(window,text=txt2,font=(font,fontsize),wraplength=length-10)
text2.configure(background=bg_color)
text2.grid(row=7,column=0,sticky='w',pady=5) # pady=5,padx=5

## Create a radio button for user to choose between two options
txt2_1="Enable the biased optimizer?"
text2_1=tk.Label(window,text=txt2_1,font=(font,fontsize),wraplength=length-10)
text2_1.configure(background=bg_color)
text2_1.grid(row=8,column=0,sticky='e') # pady=5,padx=5

bopt=tk.IntVar()
bopt.set(0)
tk.Radiobutton(window,text='No',variable=bopt,value=0,indicatoron=1,\
background=bg_button).grid(row=8,column=1,sticky='w',padx=50)
tk.Radiobutton(window,text='Yes',variable=bopt,value=1,indicatoron=1,\
background=bg_button).grid(row=8,column=1,sticky='w')
# tk.Button(window,text='Reset',command=lambda:bopt.set(0)).grid(row=8,column=1,sticky='e')
entry5=bopt

i=6
vars()['text'+str(i)],vars()['entry'+str(i)],vars()['input'+str(i)]= \
new_param(window,'text'+str(i),'entry'+str(i),'input'+str(i),row=i+4, \
dval=str(def_values[i]),text=params[6]+': ')


txt3="Optional user inputs (default options are shown):"
text3=tk.Label(window,text=txt3,font=(font,fontsize),wraplength=length-10)
text3.configure(background=bg_color)
text3.grid(row=11,column=0,sticky='w',pady=5)

basis=tk.IntVar()
basis.set(100)
tk.Radiobutton(window,text='Slater-type orbitals',variable=basis,value=100,indicatoron=1,\
background=bg_button).grid(row=12,column=0,sticky='e')
tk.Radiobutton(window,text='Gaussian-type orbitals',variable=basis,value=200,indicatoron=1,\
background=bg_button).grid(row=12,column=1,sticky='w')
entry7=basis

for j,param in enumerate(params[8:-2]):
    i=j+8
    vars()['text'+str(i)],vars()['entry'+str(i)],vars()['input'+str(i)]= \
    new_param(window,'text'+str(i),'entry'+str(i),'input'+str(i),row=13+j, \
    dval=str(def_values[i]),text=param+': ')


txt4="Visualization options:"
text4=tk.Label(window,text=txt4,font=(font,fontsize),wraplength=length-10)
text4.configure(background=bg_color)
text4.grid(row=21,column=0,sticky='w',pady=5)

for j,param in enumerate(params[-2:]):
    i=j+17
    vars()['text'+str(i)],vars()['entry'+str(i)],vars()['input'+str(i)]= \
    new_param(window,'text'+str(i),'entry'+str(i),'input'+str(i),row=22+j, \
    dval=str(def_values[-2+j]),text=param+': ')


# adds a button to close the gui
close=tk.Button(window,text='Close',command=window.quit,
                bg=bg_button,fg=fg_button)
close.place(relx=0.5,rely=0.95,anchor='s')

# Start the GUI and check if user entered a proper number
# Convert input to float if it is a number
window.mainloop()

## Checks to see if user input is proper and sets everything to either an
## integer or float
p1=1
for i in range(0,N_params):
    try:
        if i in floats:
            if i!=2:
                vars()['p{0}'.format(i)]=float(vars()['entry{0}'.format(i)].get())
                print('{0} is set with value of {1}'.format(params[i],
                                                                vars()['p{0}'.format(i)]))
            elif i==2 and p1>1:
                vars()['p{0}'.format(i)]=float(vars()['entry{0}'.format(i)].get())
                print('{0} is set with value of {1}'.format(params[i],
                                                                  vars()['p{0}'.format(i)]))
        elif i in integers:
            vars()['p{0}'.format(i)]=int(vars()['entry{0}'.format(i)].get())
            print('{0} is set with value of {1}'.format(params[i],
                                                              vars()['p{0}'.format(i)]))
    except:
        if i in floats:
            if i!=2:
                vars()['p{0}'.format(i)]=def_values[i]
                print('{0} is set with value of {1}'.format(params[i],
                                                                vars()['p{0}'.format(i)]))
            elif i==2 and p1>1:
                vars()['p{0}'.format(i)]=def_values[i]
                print('{0} is set with value of {1}'.format(params[i],
                                                                  vars()['p{0}'.format(i)]))
        elif i in integers:
            vars()['p{0}'.format(i)]=def_values[i]
            print('{0} is set with value of {1}'.format(params[i],
                                                              vars()['p{0}'.format(i)]))


# Set the line for bond length equal to zero if there's only one atom
if p1==1:
    p2=0

# Making sure the input is a positive number, add as a check later (here or above)
#se default if not positive
# if d <= 0:
#     print('You did not enter a positive number, default will be utilized.')
#     d=

# Saving user input into text file
filename='init_params.txt'
file=open(filename,'w')

for i in range(0,N_params):
    file.write(str(vars()['p{0}'.format(i)])+'\n') #+','+str(p3)+',0')
