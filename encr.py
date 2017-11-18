from PIL import Image
import tkinter as tk
from tkinter import filedialog
import hashlib 
import binascii
import textwrap
import cv2
import numpy as np

''' 
GLOBAL Constants
'''
#Lorentz intial param
x0=10
y0=8/3
z0=28
#DNA-Encoding RULE #1 A = 00, T=01, G=10, C=11
dna={}
dna["00"]="A"
dna["01"]="T"
dna["10"]="G"
dna["11"]="C"
#DNA xor
dna["AA"]=dna["TT"]=dna["GG"]=dna["CC"]="A"
dna["AG"]=dna["GA"]=dna["TC"]=dna["CT"]="G"
dna["AC"]=dna["CA"]=dna["GT"]=dna["TG"]="C"
dna["AT"]=dna["TA"]=dna["CG"]=dna["GC"]="T"

def image_selector():                           #returns path to selected image
    path = "NULL"
    root = tk.Tk()
    root.withdraw()                             # we don't want a full GUI, so keep the root window from appearing
    path = filedialog.askopenfilename()         # show an "Open" dialog box and return the path to the selected file
    if path!="NULL":
        print("Image loaded!") 
    else:
        print("Error Image not loaded!")
    return path

def split_into_rgb_channels(image):
  red = image[:,:,2]
  green = image[:,:,1]
  blue = image[:,:,0]
  return red, green, blue

#Secure key generation
def securekey (iname):
    img = Image.open(iname)
    m, n = img.size
    print("pixels: {0}  width: {1} height: {2} ".format(m*n, m, n))
    pix = img.load()          
    plainimage = list()                         #_plainimage contains all the rgb values continuously
    for y in range(n):
        for x in range(m):
            for k in range(0,3):
                plainimage.append(pix[x,y][k])    
    key = hashlib.sha256()                      #key is made a hash.sha256 object  
    key.update(bytearray(plainimage))          #image data is fed to generate digest
    return key.hexdigest()

def update_lorentz (key):
    key_bin = bin(int(key, 16))[2:].zfill(256)  #covert hex key digest to binary
    k={}                                        #key dictionary
    key_32_parts=textwrap.wrap(key_bin, 8)      #slicing key into 8 parts
    num=1
    for i in key_32_parts:
        k["k{0}".format(num)]=i
        num = num + 1
    t1 = t2 = t3 = 0
    for i in range (1,12):
        t1=t1^int(k["k{0}".format(i)],2)
    for i in range (12,23):
        t2=t2^int(k["k{0}".format(i)],2)
    for i in range (23,33):
        t3=t3^int(k["k{0}".format(i)],2)   
    global x0 ,y0, z0
    x0=x0 + t1/256            
    y0=y0 + t2/256            
    z0=z0 + t3/256            
    print("{0} {1} {2}".format(x0,y0,z0))

def decompose_matrix(iname):
    image = cv2.imread(iname)
    blue,green,red = split_into_rgb_channels(image)
    for values, channel in zip((red, green, blue), (2,1,0)):
        img = np.zeros((values.shape[0], values.shape[1]), dtype = np.uint8)
        img[:,:] = (values)
        if channel == 0:
            B = np.asmatrix(img)
        elif channel == 1:
            G = np.asmatrix(img)
        else:
            R = np.asmatrix(img)
    B = np.unpackbits(B,axis=1)
    G = np.unpackbits(G,axis=1)
    R = np.unpackbits(R,axis=1)
    print(blue.shape)
    return B,G,R

def dna_encode(b,g,r,key):
    m,n = b.shape
    r_enc= np.chararray((m,int(n/2))).astype(str)
    g_enc= np.chararray((m,int(n/2))).astype(str)
    b_enc= np.chararray((m,int(n/2))).astype(str)
    for color,enc in zip((b,g,r),(b_enc,g_enc,r_enc)):
        idx=0
        for j in range(0,m):
            for i in range(0,n,2):
                if idx==(n/2):
                    idx=0
                enc[j,idx]=dna["{0}{1}".format(color[j,i],color[j,i+1])]
                idx+=1
    #encoded key matrix
    key_bin = bin(int(key, 16))[2:].zfill(256)
    Mk = np.zeros((m,n),dtype=np.uint8)
    x=0
    for j in range(0,m):
            for i in range(0,n):
                Mk[j,i]=key_bin[x%256]
                x+=1
    
    Mk_enc=np.chararray((m,int(n/2))).astype(str)
    idx=0
    for j in range(0,m):
        for i in range(0,n,2):
            if idx==(n/2):
                idx=0
            Mk_enc[j,idx]=dna["{0}{1}".format(Mk[j,i],Mk[j,i+1])]
            idx+=1
    return b_enc,g_enc,r_enc,Mk_enc

def xor_operation(b,g,r,mk):
    m,n = b.shape
    bx=np.chararray((m,n)).astype(str)
    gx=np.chararray((m,n)).astype(str)
    rx=np.chararray((m,n)).astype(str)
    
    for i in range(0,m):
        for j in range (0,n):
            bx[i,j] = dna["{0}{1}".format(b[i,j],mk[i,j])]
            gx[i,j] = dna["{0}{1}".format(g[i,j],mk[i,j])]
            rx[i,j] = dna["{0}{1}".format(r[i,j],mk[i,j])]
    return bx,gx,rx 

#program exec9
if __name__ == "__main__":
    file_path = image_selector()
    key = securekey(file_path)
    update_lorentz(key)
    blue,green,red=decompose_matrix(file_path)
    blue_e,green_e,red_e,Mk_e=dna_encode(blue,green,red,key)
    blue_final, green_final, red_final = xor_operation(blue_e,green_e,red_e,Mk_e)


    
