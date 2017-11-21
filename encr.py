from PIL import Image
import tkinter as tk
from tkinter import filedialog
import hashlib 
import binascii
import textwrap
import cv2
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# encoding=utf8  
import sys
from importlib import reload  
from bisect import bisect_left as bsearch

''' 
GLOBAL Constants
'''
# Lorenz paramters and initial conditions
a, b, c = 10, 2.667, 28
x0, y0, z0 = 0, 0, 0

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
# Maximum time point and total number of time points
tmax, N = 1, 10000

def lorenz(X, t, a, b, c):
    x, y, z = X
    x_dot = -a*(x - y)
    y_dot = c*x - y - x*z
    z_dot = -b*z + x*y
    return x_dot, y_dot, z_dot

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
    print("pixels: {0}  width: {2} height: {1} ".format(m*n, m, n))
    pix = img.load()          
    plainimage = list()                         #_plainimage contains all the rgb values continuously
    for y in range(n):
        for x in range(m):
            for k in range(0,3):
                plainimage.append(pix[x,y][k])    
    key = hashlib.sha256()                      #key is made a hash.sha256 object  
    key.update(bytearray(plainimage))          #image data is fed to generate digest
    return key.hexdigest() ,m ,n

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
    print(blue.shape,'\n',blue,'\n')
    return B,G,R

def dna_encode(b,g,r,key):
    
    b = np.unpackbits(b,axis=1)
    g = np.unpackbits(g,axis=1)
    r = np.unpackbits(r,axis=1)
    m,n = b.shape
    r_enc= np.chararray((m,int(n/2)))
    g_enc= np.chararray((m,int(n/2)))
    b_enc= np.chararray((m,int(n/2)))
    
    print (b,m,n,'\n','\n',b.shape,b_enc.shape)
    for color,enc in zip((b,g,r),(b_enc,g_enc,r_enc)):
        idx=0
        for j in range(0,m):
            for i in range(0,n,2):
                #print("({0},{1}) ".format(idx,i))
                enc[j,idx]=dna["{0}{1}".format(color[j,i],color[j,i+1])]
                idx+=1
                if (i==n-2):
                    idx=0
                    break
    
    #encoded key matrix
    key_bin = bin(int(key, 16))[2:].zfill(256)
    Mk = np.zeros((m,n),dtype=np.uint8)
    x=0
    for j in range(0,m):
            for i in range(0,n):
                Mk[j,i]=key_bin[x%256]
                x+=1
    
    Mk_enc=np.chararray((m,int(n/2)))
    idx=0
    for j in range(0,m):
        for i in range(0,n,2):
            if idx==(n/2):
                idx=0
            Mk_enc[j,idx]=dna["{0}{1}".format(Mk[j,i],Mk[j,i+1])]
            idx+=1
    b_enc=b_enc.astype(str)
    g_enc=g_enc.astype(str)
    r_enc=r_enc.astype(str)
    Mk_enc=Mk_enc.astype(str)
    return b_enc,g_enc,r_enc,Mk_enc

def xor_operation(b,g,r,mk):
    m,n = b.shape
    bx=np.chararray((m,n))
    gx=np.chararray((m,n))
    rx=np.chararray((m,n))
    
    for i in range(0,m):
        for j in range (0,n):
            bx[i,j] = dna["{0}{1}".format(b[i,j],mk[i,j])]
            gx[i,j] = dna["{0}{1}".format(g[i,j],mk[i,j])]
            rx[i,j] = dna["{0}{1}".format(r[i,j],mk[i,j])]
    return bx,gx,rx 

def gen_chaos_seq(m,n):
    global x0,y0,z0,a,b,c,N
    N=m*n*4
    x= np.array((m,n*4))
    y= np.array((m,n*4))
    z= np.array((m,n*4))
    t = np.linspace(0, tmax, N)
    f = odeint(lorenz, (x0, y0, z0), t, args=(a, b, c))
    x, y, z = f.T
    print(len(x),m*n*4)
    x=x[:(N)]
    y=y[:(N)]
    z=z[:(N)]
    print(x,y,z)
    return x,y,z

def plot(x,y,z):
    # Plot the Lorenz attractor using a Matplotlib 3D projection
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # Make the line multi-coloured by plotting it in segments of length s which
    # change in colour across the whole time series.
    s = 10
    c = np.linspace(0,1,n)
    for i in range(0,n-s,s):
        ax.plot(x[i:i+s+1], y[i:i+s+1], z[i:i+s+1], color=(1,c[i],0), alpha=0.4)
    # Remove all the axis clutter, leaving just the curve.
    ax.set_axis_off()
    plt.show()

def sequence_indexing(x,y,z):
    n=len(x)
    fx=np.empty((n),dtype=np.uint8)
    fy=np.empty((n),dtype=np.uint8)
    fz=np.empty((n),dtype=np.uint8)
    for func, var in zip((fx,fy,fz),(x,y,z)):
        seq=sorted(var)
        for k1 in range(0,n):
            t = var[k1]
            k2 = bsearch(seq, t, 0, len(seq))
            func[k1]=int(k2)            
    print(fx,fy,fz,"\n",len(fx),len(fy),len(fz))
    return fx,fy,fz
        
def scramble(fx,fy,fz,b,r,g):
    p,q=b.shape
    size = p*q
     
    bx=b.reshape(size).astype(str)
    gx=g.reshape(size).astype(str)
    rx=r.reshape(size).astype(str)

    print("scramble:oooo ",size,len(fx))
    bx_s=np.chararray((size))
    gx_s=np.chararray((size))
    rx_s=np.chararray((size))

    for i in range(size):
            idx = fz[i]
            bx_s[i] = bx[idx]
    for i in range(size):
            idx = fy[i]
            gx_s[i] = gx[idx]
    for i in range(size):
            idx = fx[i]
            rx_s[i] = rx[idx]     
    bx_s=bx_s.astype(str)
    gx_s=gx_s.astype(str)
    rx_s=rx_s.astype(str)
    
    print(bx,"\n scramble: ",bx_s)        
    b_s=np.chararray((p,q))
    g_s=np.chararray((p,q))
    r_s=np.chararray((p,q))

    b_s=bx_s.reshape(p,q)
    g_s=gx_s.reshape(p,q)
    r_s=rx_s.reshape(p,q)

    print(b.astype(str),"\n","\n",g.astype(str),"\n","\n",r.astype(str),"\n","-----------------")
    print(b_s,"\n","\n",g_s,"\n","\n",r_s,"\n","-----------------")
    return b_s,g_s,r_s

#program exec9
if (__name__ == "__main__"):
    reload(sys)  
    file_path = image_selector()
    key,m,n = securekey(file_path)
    update_lorentz(key)
    blue,green,red=decompose_matrix(file_path)
    blue_e,green_e,red_e,Mk_e=dna_encode(blue,green,red,key)

    blue_final, green_final, red_final = xor_operation(blue_e,green_e,red_e,Mk_e)
    x,y,z=gen_chaos_seq(m,n)
    
    #plot(x,y,z)
    fx,fy,fz=sequence_indexing(x,y,z)
    blue_scrambled,green_scrambled,red_scrambled = scramble(fx,fy,fz,blue_final,red_final,green_final)
    
