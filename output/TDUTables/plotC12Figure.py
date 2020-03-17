import matplotlib.pyplot as plt


c123M = [0.349, 0.370, 0.379, 0.377, 0.370, 0.360, 0.358, 0.343, 0.328, 0.314,
         0.303, 0.293, 0.303, 0.286, 0.295, 0.284, 0.277, 0.280, 0.272, 0.294,
         0.271, 0.273, 0.291, 0.285, 0.289, 0.278]
o163M = [0.015, 0.018, 0.022, 0.022, 0.023, 0.022, 0.023, 0.028, 0.028, 0.026,
         0.024, 0.023, 0.022, 0.021, 0.020, 0.019, 0.018, 0.018, 0.018, 0.017,
         0.022, 0.019, 0.022, 0.024, 0.025, 0.024]

c124M = [0.200, 0.262, 0.300, 0.309, 0.322, 0.325, 0.331, 0.313, 0.304, 0.293,
         0.286, 0.294, 0.283, 0.283, 0.283, 0.289, 0.293, 0.300, 0.304, 0.361,
         0.352, 0.320, 0.342]
o164M = [0.003, 0.006, 0.010, 0.017, 0.019, 0.02, 0.013, 0.023, 0.023, 0.023,
         0.022, 0.022, 0.021, 0.021, 0.022, 0.024, 0.024, 0.024, 0.027, 0.027,
         0.026, 0.028, 0.030]

c125M = [0.301, 0.295, 0.295, 0.299, 0.299, 0.295, 0.295, 0.294, 0.297, 0.307,
         0.349, 0.386, 0.380, 0.346, 0.307, 0.301, 0.302, 0.317, 0.318, 0.312,
         0.337, 0.320]
o165M = [0.011, 0.014, 0.018, 0.016, 0.017, 0.016, 0.015, 0.018, 0.020, 0.017,
         0.020, 0.022, 0.019, 0.017, 0.024, 0.021, 0.015, 0.016, 0.018, 0.016,
         0.016, 0.020]

fig, axArr = plt.subplots(2, sharex = True)

axArr[0].plot(range(1, len(c123M) + 1), c123M, 'bo-', lw = 2,
        markersize = 10, label = "3 M$_\odot$")
axArr[1].plot(range(1, len(o163M) + 1), o163M, 'bo-', lw = 2,
        markersize = 10)

axArr[0].plot(range(1, len(c124M) + 1), c124M, 'g^-', lw = 2,
        markersize = 10, label = "4 M$_\odot$")
axArr[1].plot(range(1, len(o164M) + 1), o164M, 'g^-', lw = 2,
        markersize = 10)

axArr[0].plot(range(1, len(c125M) + 1), c125M, 'ys-', lw = 2,
        markersize = 10, label = "5 M$_\odot$")
axArr[1].plot(range(1, len(o165M) + 1), o165M, 'ys-', lw = 2,
        markersize = 10)

axArr[0].legend(loc = 0, numpoints = 1, prop = {'size': 12})
axArr[0].set_ylabel("$^{12}$C mass fraction")
axArr[0].minorticks_on()

axArr[1].set_ylabel("$^{16}$O mass fraction")
axArr[1].set_xlabel("TDU number")
axArr[1].minorticks_on()

fig.subplots_adjust(hspace = 0)
for ax in axArr:
    ax.xaxis.label.set_fontsize(12)
    ax.yaxis.label.set_fontsize(12)

fig.canvas.draw()


#for ii in range(len(axArr)):
    #axLabs = axArr[ii].get_yticklabels()
    #axLabs[0] = axLabs[-1] = ""
    #if ii == 0:
        #axLabs[1] = ""
    #elif ii == 1:
        #axLabs[-2] = ""
    
    #axArr[ii].set_yticklabels(axLabs)

plt.savefig('c12O16MassFractions.pdf')
plt.show()
