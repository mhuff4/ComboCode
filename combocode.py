import numpy as np
from pyhdf.SD import SD, SDC
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
from matplotlib.pyplot import plot, draw, show
from matplotlib.colors import LinearSegmentedColormap


def line_select_callback(eclick, erelease):
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))


def toggle_selector(event):
    print(' Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print(' RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print(' RectangleSelector activated.')
        toggle_selector.RS.set_active(True)


def make_plot():
    filename = 'sop_s25386.hdf'
    hdf = SD(filename, SDC.READ)
    data = hdf.select('Streak_array').get().astype('float')
    rotated = np.rot90(data, 2)
    subtract = rotated[1] - rotated[0]
    fig, current_ax = plt.subplots()
    current_ax.imshow(subtract)

    toggle_selector.RS = RectangleSelector(current_ax, line_select_callback, drawtype='box', useblit=True, button=[1, 3], minspanx=5, minspany=5, interactive=True, spancoords='pixels')
    plt.connect('key_press_event', toggle_selector)
    draw()


make_plot()
show()

startx, endx, starty, endy = input('Zoom in coords?').split()
print(startx, endx, starty, endy)

filename = 'asbo1_s80114.hdf'
hdf = SD(filename, SDC.READ)
data = hdf.select('Streak_array').get().astype('float')
rotated = np.rot90(data, 2)
subtract = rotated[1] - rotated[0]
S=subtract[int(startx):int(endx), int(starty):int(endy)]

#Create and plot the 1-D fourier transform of S (s) at a fixed time t
s = np.fft.fft(S, axis=0)
s_shift=np.fft.fftshift(s)
plt.figure(2)
plt.imshow(np.real(s_shift))
plt.title("One-D fourier transform at fixed time")


#Sum the fourier transform horizontally and plot it (s-shift)
s_shift_list=[]
for i in range(len(s_shift)):
    sum=0
    for j in range(len(s_shift[i])):
        sum+=np.real(s_shift[i,j])
    s_shift_list.append(sum)
plt.figure(3)
plt.plot(s_shift_list)

#set the frequency shift to be half the length of the box size you chose:
sz = np.rint((float(endx) - float(startx))/2.0)
#Set the negative frequencies to zero, then transform back and plot the new image
d_shift=np.zeros_like(s)
d_shift[0:int(sz),:]=s_shift[0:int(sz),:]
d=np.fft.ifftshift(d_shift)
plt.figure(4)
plt.title("One-D FM with negative frequences=0")
plt.imshow(np.real(d_shift))
D=np.fft.ifft(d,axis=0)
plt.figure(5)
plt.title('Filtered Image')
plt.imshow(np.real(D))


#Wrapped phase
W=np.angle(D)
plt.figure(6)
plt.title("Wrapped Phase")
plt.imshow(W)

#Unwrap the phase.
ratio=D[:,1:]/D[:,:-1] #the angle of dividing complex numbers is the angle difference
dphi = np.zeros_like(W)
dphi[:,1:] = np.angle(ratio)
U = -np.cumsum(dphi, 1)
plt.figure(7)
plt.title('Unwrapped fringes')
plt.imshow(U/(2*np.pi))

plt.show()