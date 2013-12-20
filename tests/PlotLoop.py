# $Id$

import matplotlib.pyplot as plt
import sys

fig = None
subplot = None
cid = None
callback = None

def Init():

	global fig, subplot
	fig = plt.figure(1)
	subplot = fig.add_subplot(111)

def RealCallback(event):

	ok = callback()
	if not ok:
		sys.exit(0)

def RegisterLoopCallback(loop):

	global cid, callback
	callback = loop
	cid = fig.canvas.mpl_connect("button_press_event", RealCallback)

def GetPlot():

	return subplot

def Clear():

	subplot.clear()

def Draw():

	fig.canvas.draw()

def Start():

	plt.show()

