import matplotlib.pyplot as plt

plt.rc('axes', grid=True)
plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)

#textsize = 9
#left, width = 0.1, 0.8
#rect1 = [left, 0.7, width, 0.2]
#rect2 = [left, 0.3, width, 0.4]
#rect3 = [left, 0.1, width, 0.2]


fig = plt.figure(facecolor='white')
axescolor  = '#f6f6f6'  # the axies background color

ax1 = fig.add_subplot(111)

#ax1 = fig.add_axes(rect1, axisbg=axescolor)  #left, bottom, width, height
#ax2 = fig.add_axes(rect2, axisbg=axescolor, sharex=ax1)
#ax2t = ax2.twinx()
#ax3  = fig.add_axes(rect3, axisbg=axescolor, sharex=ax1)

### plot the relative strength indicator
yUB =  (5, 4.5, 4.0, 2.0, 2.0)
yLB =  (1.0, 1.2, 1.4, 1.97, 1.99)
iter = (0.123, 0.324, 0.426,
        0.529, 0.6291)
yUB1 = (4.5, 4.4, 4.3, 2.4, 2.3, 1.3, 1.0)
yLB1 = (0.0, 0.1, -0.1, 0.9, 0.8, 0.8, 1.0)
iter1 = (0.01, 0.04, 0.08, 0.09, 0.10, 0.11, 0.123)

yUB2 = (4.4, 4.3, 2.4, 2.3, 1.2, 1.2, 1.2)
yLB2 = (0.1, 0.1, -0.1, 0.5, 0.4, 0.6, 1.2)
iter2 = (0.123, 0.124, 0.126, 0.22, 0.23, 0.31, 0.324)

fillcolor = 'darkgoldenrod'

ax1.plot(iter, yUB, color='b')
ax1.plot(iter, yLB, color='g')

ax1.plot(iter1, yUB1, color='r')
ax1.plot(iter1, yLB1, color='k')

ax1.plot(iter2, yUB2, color='r')
ax1.plot(iter2, yLB2, color='k')
#ax1.axhline(70, color=fillcolor)
#ax1.axhline(30, color=fillcolor)
#ax1.fill_between(r.date, rsi, 70, where=(rsi>=70), facecolor=fillcolor, edgecolor=fillcolor)
#ax1.fill_between(r.date, rsi, 30, where=(rsi<=30), facecolor=fillcolor, edgecolor=fillcolor)
#ax1.text(0.6, 0.9, '>70 = overbought', va='top', transform=ax1.transAxes, fontsize=textsize)
#ax1.text(0.6, 0.1, '<30 = oversold', transform=ax1.transAxes, fontsize=textsize)
#ax1.set_ylim(0, 6.0)
#ax1.set_yticks([2.0,7.0])
#ax1.text(0.025, 0.95, 'RSI (14)', va='top', transform=ax1.transAxes, fontsize=textsize)
#ax1.set_title('%s daily'%ticker)

# rotates and right aligns the x labels, and moves the bottom of the
# axes up to make room for them
fig.autofmt_xdate()


plt.show()
